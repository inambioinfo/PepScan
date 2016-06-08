
#script for getting protein expression per sample
# Song Cao, Ding Lab
# scao@genome.wustl.edu 

my $usage='perl protein_expression_v1.pl $msgf_dir $mgf_dir $out_dir
';

die $usage unless @ARGV == 3;
my ($msgf_dir,$mgf_dir,$out_dir) = @ARGV;
if ($msgf_dir =~/(.+)\/$/) {
        $msgf_dir = $1;
}
if ($mgf_dir =~/(.+)\/$/) {
        $mgf_dir = $1;
}
if ($out_dir =~/(.+)\/$/) {
        $out_dir = $1;
}

if(!(-d $out_dir)) { `mkdir $out_dir`; }

my %protein_expression=(); 
my %intens1=(); 
my %intens2=(); 
my %intens3=(); 
my %intens4=(); 
my %pe_control=();
my %ensp_2_gn=(); 
my $proteindb="/gscmnt/gc2108/info/medseq/proteogenomics/analyses/pipeline_step_1_peptides_detection/Quilts/databases/refseq_human_20130727_proteome/proteome.fasta"; 
my $proteinbed="/gscmnt/gc2108/info/medseq/proteogenomics/analyses/pipeline_step_1_peptides_detection/Quilts/databases/refseq_human_20130727_proteome/proteome.bed"; 
open(IN,"<$proteinbed");

while(<IN>)
{
  chomp;
  my $line=$_;
  my @ss=split("\t",$line); 
		     my @genes=split(/-/,$ss[3]);
		     $genes[0]=~s/\>//g; 
		     $ensp_2_gn{$genes[0]}=$genes[1]; 		      
		     #print $genes[0],"\t",$genes[1],"\n"; 
}
close IN; 

opendir(DH1, $msgf_dir) or die "Cannot open dir $msgf_dir: $!\n";
my @sample_dir_list_msgf = readdir DH1;
close DH1;

opendir(DH2, $mgf_dir) or die "Cannot open dir $mgf_dir: $!\n";
my @sample_dir_list_mgf = readdir DH2;
close DH2;

my $filen=(split("/",$msgf_dir))[-1]; 

print $filen,"\n"; 

for($i=0;$i<@sample_dir_list_mgf;$i++) {

   	my $sample_dir=$sample_dir_list_mgf[$i]; 
   	my $file_reporter=$mgf_dir."/".$sample_dir."/merged-d20-iTRAQ4-1000-2/merged_selected.txt"; 
   	print $file_reporter,"\n"; 
	my $runid;	
   	if($sample_dir=~/(\d+)TCGA/) { $runid=$1; }
	#print $runid,"\n";
 	#my $runid=
  	open(IN,"<$file_reporter"); 
   	while(<IN>)
   	{
     	chomp;
     	my $line=$_;  
		if($line=~/^filename/) { next; }
		else 
	  	{
	    	my @mgf_a=split("\t",$line); 
	    	my $fn=$mgf_a[0];
	    	my @fns=split(/\_/,$fn); 
	    	$fn=$fns[1]; 
	    	$fn=join("_",$fn,$fns[2],$fns[3]);
			my $repid = $fns[-1];
			$repid=~s/\.mgf//g; 
	    	my $scan=$mgf_a[1]; 
	    	my $is1=$mgf_a[4]; 
	    	my $is2=$mgf_a[5]; 
	    	my $is3=$mgf_a[6]; 
	    	my $is4=$mgf_a[7];  
			#print $repid,"\n";
			my $scan=$mgf_a[1]."R".$runid."P".$repid;
			#print $scan,"\n";
			#<STDIN>;
	    	$intens1{$fn}{$scan}=$is1*0.929+$is2*0.002; 
	    	$intens2{$fn}{$scan}=0.059*$is1+0.923*$is2+0.030*$is3+0.001*$is4; 
	    	$intens3{$fn}{$scan}=0.002*$is1+0.056*$is2+0.924*$is3+0.004*$is4; 
	    	$intens4{$fn}{$scan}=0.001*$is2+0.045*$is3+0.924*$is4; 
	  } 	
   }	    
   #last; 
 }

my $fcont=$out_dir."/RefProtE_9_11_14".$filen; 
open(IC,">$fcont"); 

for (my $i=0;$i<@sample_dir_list_msgf;$i++) {
   my  $sample_name = $sample_dir_list_msgf[$i];
   if (!($sample_name =~ /\./)) {
   my $sample_full_path = $msgf_dir."/".$sample_name;
   print $sample_full_path,"\n";
   if (-d $sample_full_path)
     {
      	opendir(DH,$sample_full_path);
      	my @sample_mzmlgz = readdir DH;
      	undef %pe_control;  
		my $runid; 
	  	if($sample_name=~/(\d+)TCGA/) { $runid=$1; }
      	foreach my $tsv_file_name (@sample_mzmlgz)
      	{
      	if (($tsv_file_name =~ /mzML\.mzid\.tsv$/))
           {
        	my $tsv_full_name=$msgf_dir."/".$sample_name."/".$tsv_file_name;
			my @file_info = split( "\_" , $tsv_full_name );
			my $repid = $file_info[-1];
        	$repid =~ s/\.mzML\.mzid\.tsv//;
        	print $tsv_full_name,"\t",$runid,"\t",$repid,"\n";
        	&read_tsv_file($tsv_full_name,$runid,$repid);
           }
      	}
      close DH;
      print IC $sample_name; 
      foreach my $gene (sort {$pe_control{$b} <=> $pe_control{$a}} keys %pe_control)
      {
        print IC "\t", $gene,",",$pe_control{$gene}; 
      }
      print IC "\n";  
     }
    }
   }

my $outfile=$out_dir."/SamplesProtE_9_11_14".$filen;
open(OUT,">$outfile");
#print OUT $sample; 
foreach my $sample (sort keys %protein_expression)
  {
   print OUT $sample; 
   foreach my $gene (sort {$protein_expression{$sample}{$b} <=> $protein_expression{$sample}{$a}}keys %{$protein_expression{$sample}})
        {
 	 print OUT "\t", $gene,",",$protein_expression{$sample}{$gene}; 
        }
   print OUT "\n"; 
 }
close OUT; 

sub read_tsv_file() 
 {
   my ($stsv,$runid,$repid)=@_;
   open(IN,"<$stsv");
   my $tsv_name=(split(/\//,$stsv))[-1];
   my %scans=();
   #my %pe_control=(); 
   while(<IN>)
   {
     chomp;
     my $line=$_;
     my $nn;
     if($line=~/scan=(\d+)/) {
     $nn=$1."R".$runID."P".$repid;
     #print $nn,"\n";
     if(!defined $scans{$nn})
     	{
          if(!($line=~/XXX_/))
          	{
            	my @lines=split(/\t/,$line);
            	my @proteins=split(/;/,$lines[9]);
            	my @genes=split(/-/,$proteins[0]);
				my $ensp=substr($genes[0],0,15);
				if($ensp=~/NP/)
				{ 
				if (defined $ensp_2_gn{$ensp}) 
				{
	    		my $gene=$ensp_2_gn{$ensp};
        		my $evalue=$lines[12];
				my $qvalue=$lines[14];
				my $fn=$lines[0];
				$fn=~s/\.mzML//g;
        		if($qvalue<=0.01 && $gene ne "" && !($gene=~/pre/) && !($gene=~/\(/) && !($gene=~/\)/)) {	 
	        		my @ss=split(/\_/,$fn);
					$fn=$ss[1];
                	$fn=join("_",$fn,$ss[2],$ss[3]);
					if(defined $intens1{$fn}{$nn})
					{ 
					my $is1=$intens1{$fn}{$nn}; 
					my $is2=$intens2{$fn}{$nn}; 
					my $is3=$intens3{$fn}{$nn}; 
					my $is4=$intens4{$fn}{$nn}; 
					my $normis=$is1+$is2+$is3+$is4;
					$protein_expression{$ss[1]}{$gene}+=$is1/$normis; 
					$protein_expression{$ss[2]}{$gene}+=$is2/$normis; 
					$protein_expression{$ss[3]}{$gene}+=$is3/$normis;    
					$pe_control{$gene}+=$is4/$normis;
		#print $pe_control{$gene},"\n";	
					}      										
           		}
        		}
			}
        	}
		}
         $scans{$nn}=1;
      }
    }
  }   
