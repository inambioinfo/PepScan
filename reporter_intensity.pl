
#script for getting reporter intensity
# Song Cao, Ding Lab
# scao@genome.wustl.edu 

my $usage='perl reporter_intensity.pl $run_dir
'; 

die $usage unless @ARGV == 1;
my ( $run_dir) = @ARGV;
if ($run_dir =~/(.+)\/$/) {
        $run_dir = $1;
}

opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
 
my @sample_dir_list_mgf = readdir DH;
close DH;

my $HOME = $ENV{HOME};

if (! -d $HOME."/shfiles") {
        `mkdir $HOME"/shfiles"`;
}

my $job_files_dir = $HOME."/shfiles";

#store SGE output and error files here
if (! -d $HOME."/LSF_DIR") {
        `mkdir $HOME"/LSF_DIR"`;
}

my $lsf_file_dir = $HOME."/LSF_DIR";

my $IonSR="/gscuser/scao/scripts/peptidesDect/mgf_select_reporter_ions_only.pl";
my $sample_name; 
my $sample_full_path; 
my $job_name;
my $hold_job_file; 
my $current_job_file; 
 
for($i=0;$i<@sample_dir_list_mgf;$i++) {

     $sample_name = $sample_dir_list_mgf[$i];

     if (!($sample_name =~ /\./))
       {
 	$sample_full_path = $run_dir."/".$sample_name;
  	print $sample_full_path,"\n";
        if (-d $sample_full_path)
     	 {
          $job_name=$sample_name;
          #$sample_full_name=$sample_full_path;
          print $yellow, "\nSubmitting jobs for the sample ",$sample_full_path, "...",$normal, "\n";
          &submit_job_array_mgf(1,$sample_full_path);
	 # <STDIN>;
         }
      	}
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
