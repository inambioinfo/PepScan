# PepScan Pipeline
# Song Cao, Ding Lab
# scao@genome.wustl.edu 

my $version = 1.2;
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m"; 
(my $usage = <<OUT) =~ s/\t+//g;
This script will run the piptide identification pipeline on LSF cluster.
Pipeline version: $version
$yellow		Usage: perl $0 <run_folder> <run_dir_mgf> <run_dir_out> <type_mass> <step_number> $normal
<run_folder> = full path of the folder holding files of mzml and vcf file
<run_folder> = full path of the folder holding the mgf file
<run_dir_out> = full path of the folder which outputs the protein expresion and peptides variants 
<type_mass> = 0, 1 for phosphoproteome and proteome, respectively. 
<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)
[1] Run Quilts for the annotation of snvs and in-house script for the annotation of indels
[2] Create fasta DB for MSGF+
[3] Run MSGF+ for one sample
[4] Run MSGF+ for multiple samples
[5] Convert mzid to tsv file
[6] Extract reporter ion intensity
[7] Get protein expression
[8] Get the number of SAAVs, independent on mass type
$normal
OUT

die $usage unless @ARGV == 5;
my ( $run_dir,$run_dir_mgf, $run_dir_out, $type_mass,$step_number ) = @ARGV;
if ($run_dir =~/(.+)\/$/) {
	$run_dir = $1;
}

die $usage unless ((($step_number >=0)&&($step_number <= 8)));
#####################################################################################
# values need to be modified to adapt to local environment
my $email = "scao\@genome.wustl.edu";
#####################################################################################
# everything else below should be automated
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-1];

my $run_script_path = `dirname $0`;
#my $run_script_path = `pwd`;
chomp $run_script_path;
$run_script_path = "/usr/bin/perl ".$run_script_path."/";

#my $file_number_of_msgf=24;

#store job files here

if (! -d $HOME."/shfiles") {
	`mkdir $HOME"/shfiles"`;
}

my $job_files_dir = $HOME."/shfiles";

#store SGE output and error files here
if (! -d $HOME."/LSF_DIR") {
	`mkdir $HOME"/LSF_DIR"`;
}

my $lsf_file_dir = $HOME."/LSF_DIR";
my $mod_msgf; 

if($type_mass==0) { $mod_msgf="/gscmnt/gc2108/info/medseq/proteogenomics/analyses/pipeline_step_1_peptides_detection/cptac/breast/MSGF+_modifications/MSGFDB_CompRef_Phospho_iTRAQ_20ppm.txt"; }
if($type_mass==1) { $mod_msgf="/gscmnt/gc2108/info/medseq/proteogenomics/analyses/pipeline_step_1_peptides_detection/cptac/breast/MSGF+_modifications/MSGFDB_CompRef_iTRAQ_20ppm.txt"; }
# obtain script path

my $run_script_path = `dirname $0`;
chomp $run_script_path;
$run_script_path = "/usr/bin/perl ".$run_script_path."/";
my $MSGFPLUS="/gscmnt/gc2108/info/medseq/proteogenomics/analyses/pipeline_step_1_peptides_detection/MSGF+/bin/MSGFPlus.jar";
my $QUILTS1="./quilts-WU/v1.0/run_quilts_tophat.pl"; 
my $reporter="reporter_intensity.pl"; 
my $IonSR="mgf_select_reporter_ions_only.pl";
my $Rindel="protein_seq_for_indel_using_bed.pl";
my $PBED="/gscuser/scao/gc2108/info/medseq/proteogenomics/databases/refseq_human_20130727_proteome/proteome.bed";

my $current_job_file = "no";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";
my $msgf_fasta; 
my $job_name="";
my $run_dir_msgf; 
my $run_dir_quilts;
#my $run_dir_mgf; 
#my $run_dir_out="/gscmnt/gc3027/info/medseq/proteomics/cptac/breast/proteindata"; 

$run_dir_quilts=$run_dir."/quilts"; 

$run_dir_indel=$run_dir."/indel";

if(!(-d $run_dir_quilts)) { print "no quilts directory\n"; exit(2);  }
if(!(-d $run_dir_indel))  { print "no indel directory\n"; exit(2); }

if($type_mass==0) { $run_dir_msgf=$run_dir."/phosphoproteome"; }
else { 
     if($type_mass==1) { $run_dir_msgf=$run_dir."/proteome"; }
     else { die $usage; }
     }

if($type_mass==0) { $run_dir_mgf=$run_dir_mgf."/phosphoproteome"; }
else {
     if($type_mass==1) { $run_dir_mgf=$run_dir_mgf."/proteome"; }
     else { die $usage; }
     }

opendir(DH, $run_dir_quilts) or die "Cannot open dir quilts $run_dir: $!\n";
my @sample_dir_list_quilts = readdir DH;
close DH;

opendir(DH, $run_dir_indel) or die "Cannot open dir indel $run_dir: $!\n";
my @sample_dir_list_indel = readdir DH;
close DH;

opendir(DH, $run_dir_msgf) or die "Cannot open dir msgf $run_dir: $!\n";
my @sample_dir_list_msgf = readdir DH;
close DH;


opendir(DH, $run_dir_mgf) or die "Cannot open dir mgf $run_dir_mgf: $!\n";
my @sample_dir_list_mgf = readdir DH;
close DH;

my $HOME = $ENV{HOME};

if ($step_number < 16) {

      if($step_number==0 || $step_number == 1)
        {

                for (my $i=0;$i<@sample_dir_list_quilts;$i++) 
				{
                	$sample_name = $sample_dir_list_quilts[$i];
                	if (!($sample_name =~ /\./)) 
						{
                     		$sample_full_path = $run_dir_quilts."/".$sample_name;
                        	print $sample_full_path,"\n";                   
                        	if (-d $sample_full_path) 
                        	 { 
                                   $job_name=$sample_name; 
                                   $sample_full_name=$sample_full_path; 
                                   print $red, "\nSubmitting quiltsjobs for the sample ",$sample_full_name, "...",$normal, "\n";
                                   &submit_job_array_quilts(1,$sample_full_name);}
								   #&submit_job_array_indel(1,$sample_full_name);
                   		}
         		}

				for (my $i=0;$i<@sample_dir_list_indel;$i++) 
				{
                	$sample_name = $sample_dir_list_indel[$i];
                	if (!($sample_name =~ /\./))
                    	{
                        	$sample_full_path = $run_dir_indel."/".$sample_name;
                        	print $sample_full_path,"\n";
                        	if (-d $sample_full_path)
                             	{
                                   $job_name=$sample_name;
                                   $sample_full_name=$sample_full_path."/".$sample_name.".vcf";
                                   print $red, "\nSubmitting indel jobs for the sample ",$sample_full_name, "...",$normal, "\n";
                                   &submit_job_array_indel(1,$sample_full_name);
                                 }
                         }
                 }
			
		}

       if($step_number==0 || $step_number == 2)
         {
                for (my $i=0;$i<@sample_dir_list_msgf;$i++) 
				{
                $sample_name = $sample_dir_list_msgf[$i];
                if (!($sample_name =~ /\./)) {
                        $sample_full_path = $run_dir_msgf."/".$sample_name;
                        print $sample_full_path,"\n";
			$job_name=$sample_name; 
                        if (-d $sample_full_path)
                        {
                           	$msgf_fasta=$run_dir_msgf."/".$sample_name."/".$sample_name.".fasta";
				print $gray, "\nSubmitting jobs for creating fasta db for msgf the sample ",$sample_full_path, "...",$normal, "\n";		
                                        &submit_job_array_db(1);
                                                                }
                                        }
		 }
          }
	
	if($step_number==0 || $step_number == 3)
	{
		for (my $i=0;$i<@sample_dir_list_msgf;$i++) {
		$sample_name = $sample_dir_list_msgf[$i];
		if (!($sample_name =~ /\./)) {
			$sample_full_path = $run_dir_msgf."/".$sample_name;
 			print $sample_full_path,"\n"; 			
			if (-d $sample_full_path) 
			{ 
		           opendir(DH,$sample_full_path); 	
		           my @sample_mzmlgz = readdir DH; 
			   {
			    $msgf_fasta=$run_dir_msgf."/".$sample_name."/".$sample_name.".fasta"; 

       			   foreach my $gz_file_name (@sample_mzmlgz)
				     {
                                     	if (($gz_file_name =~ /\.mzML\.gz$/ || $gz_file_name =~ /\.mzML$/)) 
					 {
     				    	  $gz_file_name=~s/\.gz$//g;
				    	  $job_name=$gz_file_name; 
                                    	  $sample_full_name=$run_dir_msgf."/".$sample_name."/".$gz_file_name; 
					#$msgf_fasta=$run_dir_msgf."/".$sample_name."/".$sample_name.".fasta"; 
				    	  print $cyan, "\nSubmitting jobs for initial msgf run: sample ",$sample_full_name, "...",$normal, "\n";
					  &submit_job_single_init_msgf(1,$sample_full_name);
					  last;  
                                          }
					}
				      }
			    }
			  }
			  close DH; 
			}
	    }


        if($step_number==0 || $step_number ==4)
        {
                for (my $i=0;$i<@sample_dir_list_msgf;$i++) {
                $sample_name = $sample_dir_list_msgf[$i];
                if (!($sample_name =~ /\./)) {
                        $sample_full_path = $run_dir_msgf."/".$sample_name;
                        print $sample_full_path,"\n";
                        if (-d $sample_full_path)
                        {
                           opendir(DH,$sample_full_path);
                           my @sample_mzmlgz = readdir DH;
                           {
                            $msgf_fasta=$run_dir_msgf."/".$sample_name."/".$sample_name.".fasta";
                           foreach my $gz_file_name (@sample_mzmlgz)
                                     {
                                        if (($gz_file_name =~ /\.mzML\.gz$/ || $gz_file_name =~ /\.mzML$/))
                                        {
                                        $gz_file_name=~s/\.gz$//g;
                                        $job_name=$gz_file_name;
                                        $sample_full_name=$run_dir_msgf."/".$sample_name."/".$gz_file_name;
                                        #$msgf_fasta=$run_dir_msgf."/".$sample_name."/".$sample_name.".fasta"; 
                                        print $cyan, "\nSubmitting jobs for running msgf: sample ",$sample_full_name, "...",$normal, "\n";
                                        &submit_job_array_msgf(1,$sample_full_name);
                                           }

                                        }
                                      }
                            }
                          }
                          close DH;
                        }
            }
 


        if($step_number==0 || $step_number ==5)
        {
                for (my $i=0;$i<@sample_dir_list_msgf;$i++) {
                $sample_name = $sample_dir_list_msgf[$i];
                if (!($sample_name =~ /\./)) {
                        $sample_full_path = $run_dir_msgf."/".$sample_name;
                        print $sample_full_path,"\n";
                        if (-d $sample_full_path)
                        {
                           opendir(DH,$sample_full_path);
                           my @sample_mzmlgz = readdir DH;
                           {
                            $msgf_fasta=$run_dir_msgf."/".$sample_name."/".$sample_name.".fasta";
                           foreach my $gz_file_name (@sample_mzmlgz)
                                     {
                                        if (($gz_file_name =~ /\.gz$/ || $gz_file_name =~ /\.mzML$/))
                                        {
                                        $gz_file_name=~s/\.gz$//g;
                                        $job_name=$gz_file_name;
                                        $sample_full_name=$run_dir_msgf."/".$sample_name."/".$gz_file_name;
                                        #$msgf_fasta=$run_dir_msgf."/".$sample_name."/".$sample_name.".fasta"; 
                                        print $cyan, "\nSubmitting jobs for getting tsv files: sample ",$sample_full_name, "...",$normal, "\n";
                                        &submit_job_array_tsv(1,$sample_full_name);
                                           }

                                        }
                                      }
                            }
                          }
                          close DH;
                        }
            }

	 if($step_number==0 || $step_number ==6)
         {
		for($i=0;$i<@sample_dir_list_mgf;$i++) {

     		$sample_name = $sample_dir_list_mgf[$i];

     		if (!($sample_name =~ /\./))
       		{
        	$sample_full_path = $run_dir_mgf."/".$sample_name;
        	print $sample_full_path,"\n";
        	if (-d $sample_full_path)
         	{
          	$job_name=$sample_name;
          #$sample_full_name=$sample_full_path;
          	print $yellow, "\nSubmitting mgf jobs for the sample ",$sample_full_path, "...",$normal, "\n";
          	&submit_job_array_mgf(1,$sample_full_path);
         # <STDIN>;
         	}
        	}
      		}
	  }

       if($step_number==0 || $step_number ==7)
	 	{
	   		if($type_mass==1)
	   		{
           		print $purple, "\nSubmitting jobs for getting protein expression ","...",$normal, "\n";
	   			&submit_job_protein_expression(1); }
	  # else 
	   #{ 
	    # print $purple, "\nSubmitting job..","phosphoprotein..",$normal, "\n";
            # &submit_job_phosphoprotein_expression(1);
	   #	}
         }

	 if($step_number==0 || $step_number ==8)
         {
	   print $purple, "\nSubmitting jobs for getting the number of SAAVS ","...",$normal, "\n";
           &submit_job_variant(1);
         }


}
	

sub submit_job_array_quilts{
        my ($step_by_step,$sid) = @_;
        if ($step_by_step) {
                $hold_job_file = "no";
        }else{
                $hold_job_file = $current_job_file;
        }
        $current_job_file = "j1_".$job_name."_quilts.sh";
        open (QUILTS, ">$job_files_dir/$current_job_file") or die $!;
        print QUILTS "#!/bin/bash\n\n";
        print QUILTS "#BSUB -n 1\n";
        print QUILTS "#BSUB -R \"rusage[mem=12000]\"","\n";
        print QUILTS "#BSUB -M 12000000\n";
        print QUILTS "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
        print QUILTS "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
        print QUILTS "#BSUB -J $current_job_file\n";
        print QUILTS "SAMPLE_DIR=".$sample_full_path,"\n";
        print QUILTS "PROTEIN_DIR=".$sample_full_path."/protein","\n"; 
        #print MSGF "#\$ -t 1-$file_number_of_msgf:1","\n"; #must be a decimal number, the value must be determined when this job file is generated. cannot be a variable
        print QUILTS 'if [ ! -d $PROTEIN_DIR ]', "\n";
        print QUILTS "then\n";
        print QUILTS "perl $QUILTS1 \${SAMPLE_DIR}","\n";
        print QUILTS "fi\n";
        close QUILTS;
        $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
        system ($bsub_com);
}


sub submit_job_array_indel{
        my ($step_by_step,$sid) = @_;
        if ($step_by_step) {
                $hold_job_file = "no";
        }else{
                $hold_job_file = $current_job_file;
        }
        $current_job_file = "j1indel_".$job_name."_indel.sh";
        open (INDEL, ">$job_files_dir/$current_job_file") or die $!;
        print INDEL "#!/bin/bash\n\n";
        print INDEL "#BSUB -n 1\n";
        print INDEL "#BSUB -R \"rusage[mem=12000]\"","\n";
        print INDEL "#BSUB -M 12000000\n";
        print INDEL "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
        print INDEL "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
        print INDEL "#BSUB -J $current_job_file\n";
        #print MSGF "#\$ -t 1-$file_number_of_msgf:1","\n"; #must be a decimal number, the value must be determined when this job file is generated. cannot be a variable
        print INDEL "perl $Rindel $PBED $HDIR $sid","\n";
        close INDEL;
        $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
        system ($bsub_com);
}

sub submit_job_array_db{
        my ($step_by_step) = @_;
        if ($step_by_step) {
                $hold_job_file = "no";
        }else{
                $hold_job_file = $current_job_file;
        }
        $current_job_file = "j2_".$job_name."_db.sh";
        open (FASTADB, ">$job_files_dir/$current_job_file") or die $!;
        print FASTADB "#!/bin/bash\n\n";
        print FASTADB "#BSUB -n 1\n";
        print FASTADB "#BSUB -R \"rusage[mem=12000]\"","\n";
        print FASTADB "#BSUB -M 12000000\n";
        print FASTADB "#BSUB -J $current_job_file\n";
        print FASTADB "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
        print FASTADB "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
	print FASTADB "rm -rf $sample_full_path/$sample_name.revCat*\n"; 
        print FASTADB "         ".$run_script_path."db_snpindel.pl $sample_full_path\n";
        close FASTADB;
        $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
        system ($bsub_com);
    }      

sub submit_job_single_init_msgf{
        my ($step_by_step,$sid) = @_;
        if ($step_by_step) {
                $hold_job_file = "no";
        }else{
                $hold_job_file = $current_job_file;
        }
        $current_job_file = "j3_".$job_name."_init_msgf.sh";
        open (MSGF1, ">$job_files_dir/$current_job_file") or die $!;
        print MSGF1 "#!/bin/bash\n\n";
        print MSGF1 "#BSUB -n 1\n";
        print MSGF1 "#BSUB -R \"rusage[mem=12000]\"","\n";
        print MSGF1 "#BSUB -M 12000000\n";
        print MSGF1 "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
        print MSGF1 "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
        print MSGF1 "#BSUB -J $current_job_file\n";
        print MSGF1 "SAMPLE_DIR=".$sample_full_path,"\n";
        print MSGF1 "msgfOUT=$sid.mzid\n";#full path
        print MSGF1 "msgftsv=$sid.mzid.tsv\n";
        print MSGF1 "msgfzip=$sid.gz\n";
        #print MSGF1 "msgffasta=$sid.fasta\n";
        print MSGF1 "msgfIN=$sid\n";
        print MSGF1 'if [ ! -f $msgfIN ]',"\n";
        print MSGF1 "then\n";
        print MSGF1 'gunzip $msgfzip',"\n";
        print MSGF1 "fi\n";
        print MSGF1 'if [ -s $msgfIN ]',"\n";
        print MSGF1 "then\n";
# Note(rmashl): for TGI systems, thread count should match BSUB -n value
        print MSGF1 "            java -Xmx10G  -jar $MSGFPLUS -thread 1 -t 20ppm  -ti \"-1,2\" -ntt 2  -tda 1 -n 5 -minLength 6 -maxLength 50 -minCharge 1 -maxCharge 6 -mod $mod_msgf -d $msgf_fasta -s \${msgfIN} -o \${msgfOUT}","\n";
        print MSGF1 '            tail -5 ${msgfOUT}|grep MzIdentML',"\n";
        print MSGF1 '            CHECK1=$?',"\n";
        print MSGF1 '            while [ ${CHECK1} -eq 1 ]',"\n";
        print MSGF1 "            do\n";
        print MSGF1 "                    java -Xmx10G  -jar $MSGFPLUS -thread 1 -t 20ppm  -ti \"-1,2\" -ntt 2  -tda 1 -n 5 -minLength 6 -maxLength 50 -minCharge 1 -maxCharge 6 -mod $mod_msgf -d $msgf_fasta -s  \${msgfIN} -o  \${msgfOUT}","\n";
        print MSGF1 '                    tail -5 ${msgfOUT}|grep MzIdentML',"\n";
        print MSGF1 '                    CHECK1=$?',"\n";
        print MSGF1 "            done\n";
        print MSGF1 "fi";
        close MSGF1;
        $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
        system ($bsub_com);
}
 
sub submit_job_array_msgf{
	my ($step_by_step,$sid) = @_;
	if ($step_by_step) {
		$hold_job_file = "no";
	}else{
		$hold_job_file = $current_job_file;
	}
	$current_job_file = "j4_".$job_name."_msgf.sh";
	open (MSGF, ">$job_files_dir/$current_job_file") or die $!;
	print MSGF "#!/bin/bash\n\n"; 
	print MSGF "#BSUB -n 1\n"; 
 	print MSGF "#BSUB -R \"rusage[mem=12000]\"","\n"; 
	print MSGF "#BSUB -M 12000000\n"; 
	print MSGF "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n"; 
	print MSGF "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n"; 
	print MSGF "#BSUB -J $current_job_file\n"; 
	print MSGF "SAMPLE_DIR=".$sample_full_path,"\n";
	#print MSGF "#\$ -t 1-$file_number_of_msgf:1","\n"; #must be a decimal number, the value must be determined when this job file is generated. cannot be a variable
	print MSGF "msgfOUT=$sid.mzid\n";#full path
    print MSGF "msgftsv=$sid.mzid.tsv\n";
	print MSGF "msgfzip=$sid.gz\n";
	#print MSGF "msgffasta=$sid.fasta\n"; 
    print MSGF "msgfIN=$sid\n";
    print MSGF 'if [ ! -f $msgfIN ]',"\n";
	print MSGF "then\n";
	print MSGF 'gunzip $msgfzip',"\n";
    print MSGF "fi\n";
	print MSGF 'if [ -s $msgfIN ]',"\n"; 
	print MSGF "then\n";
	print MSGF "		java -Xmx10G  -jar $MSGFPLUS  -thread 1 -t 20ppm  -ti \"-1,2\" -ntt 2  -tda 1 -n 5 -minLength 6 -maxLength 50 -minCharge 1 -maxCharge 6 -mod $mod_msgf -d $msgf_fasta -s \${msgfIN} -o \${msgfOUT}","\n";
	print MSGF '		tail -5 ${msgfOUT}|grep MzIdentML',"\n";
	print MSGF '		CHECK1=$?',"\n";
	print MSGF '		while [ ${CHECK1} -eq 1 ]',"\n";
	print MSGF "		do\n";
	print MSGF "			java -Xmx10G  -jar $MSGFPLUS -thread 1 -t 20ppm  -ti \"-1,2\" -ntt 2  -tda 1 -n 5 -minLength 6 -maxLength 50 -minCharge 1 -maxCharge 6 -mod $mod_msgf -d $msgf_fasta -s  \${msgfIN} -o  \${msgfOUT}","\n";
	print MSGF '			tail -5 ${msgfOUT}|grep MzIdentML',"\n";
	print MSGF '			CHECK1=$?',"\n";
	print MSGF "		done\n";
        #print MSGF "java -Xmx10G -XX:PermSize=32M -XX:MaxPermSize=256M -cp $MSGFPLUS  edu.ucsd.msjava.ui.MzIDToTsv -i  \${msgfOUT} -showDecoy 1 -o \${msgftsv}","\n"; 
	print MSGF "fi";
	close MSGF;
	$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
	system ($bsub_com);
     }

sub submit_job_array_tsv{
        my ($step_by_step,$sid) = @_;
        if ($step_by_step) {
                $hold_job_file = "no";
        }else{
                $hold_job_file = $current_job_file;
        }
        $current_job_file = "j5_".$job_name."_tsv.sh";
        open (MSGF2, ">$job_files_dir/$current_job_file") or die $!;
        print MSGF2 "#!/bin/bash\n\n";
        print MSGF2 "#BSUB -n 1\n";
        print MSGF2 "#BSUB -R \"rusage[mem=12000]\"","\n";
        print MSGF2 "#BSUB -M 12000000\n";
        print MSGF2 "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
        print MSGF2 "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
        print MSGF2 "#BSUB -J $current_job_file\n";
        print MSGF2 "SAMPLE_DIR=".$sample_full_path,"\n";
        #print MSGF "#\$ -t 1-$file_number_of_msgf:1","\n"; #must be a decimal number, the value must be determined when this job file is generated. cannot be a variable
        print MSGF2 "msgfOUT=$sid.mzid\n";#full path
        print MSGF2 "msgftsv=$sid.mzid.tsv\n";
        print MSGF2 "msgfzip=$sid.gz\n";
        #print MSGF2 "msgffasta=$sid.fasta\n";
        print MSGF2 "msgfIN=$sid\n";
        print MSGF2 'if [ -f $msgftsv ]',"\n";
        print MSGF2 "then\n";
        print MSGF2 'rm $msgftsv',"\n";
        print MSGF2 "fi\n";
        print MSGF2 'if [ -s $msgfOUT ]',"\n";
        print MSGF2 "then\n";
	#print MSGF2 "java -Xmx3G -XX:PermSize=32M -XX:MaxPermSize=1000M -cp $MSGFPLUS  edu.ucsd.msjava.ui.MzIDToTsv -i  \${msgfOUT} -showDecoy 1 -o \${msgftsv}","\n";
    	print MSGF2 'while [ ! -s $msgftsv ]',"\n";
	print MSGF2 "do\n";
        print MSGF2 "java -Xmx10G -XX:PermSize=32M -XX:MaxPermSize=256M -cp $MSGFPLUS  edu.ucsd.msjava.ui.MzIDToTsv -i  \${msgfOUT} -showDecoy 1 -o \${msgftsv}","\n";
	print MSGF2 "done\n";
        print MSGF2 "fi";
        close MSGF2;
        $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
        system ($bsub_com);
     }

sub submit_job_array_mgf{
        my ($step_by_step,$SAMPLE_DIR) = @_;
        if ($step_by_step) {
                $hold_job_file = "no";
        }else{
                $hold_job_file = $current_job_file;
        }
        $current_job_file = "j6_".$job_name."_mgf.sh";
        open (MGF, ">$job_files_dir/$current_job_file") or die $!;
        print MGF "#!/bin/bash\n\n";
        print MGF "#BSUB -n 1\n";
        print MGF "#BSUB -R \"rusage[mem=1200]\"","\n";
        print MGF "#BSUB -M 1200000\n";
        print MGF "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
        print MGF "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
        print MGF "#BSUB -J $current_job_file\n";
        #print MGF "SAMPLE_DIR=".$sample_full_path,"\n";
        #print MGF "PROTEIN_DIR=".$sample_full_path."/protein","\n";
        #print MSGF "#\$ -t 1-$file_number_of_msgf:1","\n"; #must be a decimal number, the value must be determined when this job file is generated. cannot be a variable
        print MGF "if [ -d $SAMPLE_DIR ]", "\n";
        print MGF "then\n";
        print MGF "perl $IonSR $SAMPLE_DIR","\n";
        print MGF "fi\n";
        close MGF;
        $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
        system ($bsub_com);
}

sub submit_job_protein_expression{
        my ($step_by_step) = @_;
        if ($step_by_step) {
                $hold_job_file = "no";
        }else{
                $hold_job_file = $current_job_file;
        }
        $current_job_file = "j7_"."proteinexpression.sh";
        open (PRO, ">$job_files_dir/$current_job_file") or die $!;
        print PRO "#!/bin/bash\n\n";
        print PRO "#BSUB -n 1\n";
        print PRO "#BSUB -R \"rusage[mem=12000]\"","\n";
        print PRO "#BSUB -M 12000000\n";
        print PRO "#BSUB -J $current_job_file\n";
        print PRO "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
        print PRO "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
        print PRO "         ".$run_script_path."protein_expression_refseq.pl $run_dir_msgf $run_dir_mgf $run_dir_out\n";
        close PRO;
        $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
        system ($bsub_com);
    }

#sub submit_job_phosphoprotein_expression{
 #       my ($step_by_step) = @_;
  #      if ($step_by_step) {
   #             $hold_job_file = "no";
    #    }else{
     #           $hold_job_file = $current_job_file;
      #  }
      #  $current_job_file = "j7_"."phosphoproteinexpression.sh";
       # open (PRO, ">$job_files_dir/$current_job_file") or die $!;
       # print PRO "#!/bin/bash\n\n";
       # print PRO "#BSUB -n 1\n";
       # print PRO "#BSUB -R \"rusage[mem=12000]\"","\n";
       # print PRO "#BSUB -M 12000000\n";
       # print PRO "#BSUB -J $current_job_file\n";
        #print PRO "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
        #print PRO "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
        #print PRO "         ".$run_script_path."phosphoprotein_expression_v1.pl $run_dir_msgf $run_dir_mgf $run_dir_out\n";
        #close PRO;
        #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
        #system ($bsub_com);
    #}

sub submit_job_variant{
        my ($step_by_step) = @_;
        if ($step_by_step) {
                $hold_job_file = "no";
        }else{
                $hold_job_file = $current_job_file;
        }
        $current_job_file = "j8_"."variant.sh";
        open (VARI, ">$job_files_dir/$current_job_file") or die $!;
        print VARI "#!/bin/bash\n\n";
        print VARI "#BSUB -n 1\n";
        print VARI "#BSUB -R \"rusage[mem=12000]\"","\n";
        print VARI "#BSUB -M 12000000\n";
        print VARI "#BSUB -J $current_job_file\n";
        print VARI "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
        print VARI "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
        print VARI "         ".$run_script_path."peptide_table_msgf.pl $run_dir 0 $run_dir_out\n";
	print VARI "         ".$run_script_path."peptide_table_msgf.pl $run_dir 1 $run_dir_out\n";
	print VARI "cat $run_dir_out/Phosphoproteome_msgf.txt $run_dir_out/Proteome_msgf.txt > $run_dir_out/DB_MSGF.txt\n"; 
	print VARI "         ".$run_script_path."get_confirmed_snps_v1.pl $run_dir $run_dir_out/DB_MSGF.txt $run_dir_out\n";
        close VARI;
        $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
        system ($bsub_com);
    }

