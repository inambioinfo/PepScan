#!/usr/local/bin/perl
use strict;

my $error=0;
my $dir="";
my $mz_error="";
my $type="";
my $min_intensity="";
my $min_reporters="";

my $usage = '
mgf intensity report
Please give the directory which contains mgf files. 
';

die $usage unless scalar @ARGV >= 1;

if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $dir="."; } 
if ($ARGV[1]=~/\w/) { $mz_error=$ARGV[1];} else { $mz_error=20; } 
if ($ARGV[2]=~/\w/) { $type=$ARGV[2];} else { $type="iTRAQ4"; } 
if ($ARGV[3]=~/\w/) { $min_intensity=$ARGV[3];} else { $min_intensity=1000; } 
if ($ARGV[4]=~/\w/) { $min_reporters=$ARGV[4];} else { $min_reporters=2; } 

my $com1="rm $dir/merged-d20-iTRAQ4-1000-2 -rf";
my $com2="rm $dir/not-selected-d20-iTRAQ4-1000-2 -rf";
my $com3="rm $dir/selected-d20-iTRAQ4-1000-2 -rf";

print $com1,"\n";
print $com2,"\n";
print $com3,"\n";
system($com1);
system($com2);
system($com3);
#<STDIN>;

my @reporters=();

if ($type=~/^iTRAQ8$/)
{
	#@reporters=(113.1078,114.1112,115.1083,116.1116,117.1150,118.1120,119.1154,121.1221);
	@reporters=(113.1078,114.1112,115.1082,116.1116,117.1150,118.1120,119.1153,121.1220);
}
else
{
	if ($type=~/^iTRAQ4$/ or $type=~/^iTRAQ/)
	{
		@reporters=(114.1112,115.1083,116.1116,117.1150);
	}
	else
	{
		if ($type=~/^TMT6$/ or $type=~/^TMT/)
		{
			#@reporters=(126.127725,127.124760,128.134433,129.131468,130.141141,131.138176);
			@reporters=(126.127725,127.131079,128.134433,129.137787,130.141141,131.138176);
		}
		else
		{
			print "Error\n";
			$error=1;
		}
	}
}

if ($error==0)
{
	my $count_all_spectra=0;
	my $count_spectra=0;
	my $max_max=0;
	my %stat=();
	mkdir(qq!$dir/selected-d$mz_error-$type-$min_intensity-$min_reporters!);
	mkdir(qq!$dir/selected-d$mz_error-$type-$min_intensity-$min_reporters/reporter_cal!);
	mkdir(qq!$dir/not-selected-d$mz_error-$type-$min_intensity-$min_reporters!);
	mkdir(qq!$dir/merged-d$mz_error-$type-$min_intensity-$min_reporters!);

	if (opendir(dir,"$dir"))
	{
		my @allfiles=readdir dir;
		closedir dir;
		#open (OUT_MERGED,">$dir/merged-d$mz_error-$type-$min_intensity-$min_reporters/merged_selected.mgf");
		#open (OUT_CAL_MERGED,">$dir/merged-d$mz_error-$type-$min_intensity-$min_reporters/merged_selected_cal.mgf");
		open (OUT_MERGED_TABLE,">$dir/merged-d$mz_error-$type-$min_intensity-$min_reporters/merged_selected.txt");
		print OUT_MERGED_TABLE qq!filename\tscan\trt\tMS1_intensity!;
		foreach my $reporter (@reporters)
		{
			my $reporter_=int($reporter);
			print OUT_MERGED_TABLE qq!\t$type-$reporter_!;
		}
		print OUT_MERGED_TABLE qq!\t$type-sum\n!;
		open (OUT_MERGED_NOT,">$dir/merged-d$mz_error-$type-$min_intensity-$min_reporters/merged_not_selected.txt");
		open (OUT_COUNT,">$dir/merged-d$mz_error-$type-$min_intensity-$min_reporters/count.txt");
		print OUT_COUNT "filename\tselected_spectra\tspectra\n";
		open (OUT_DM_INT,">$dir/merged-d$mz_error-$type-$min_intensity-$min_reporters/dm_intensity.txt");
		print OUT_DM_INT "reporter\tdm\tintensity\n";
		open (OUT_DM_INT_,">$dir/merged-d$mz_error-$type-$min_intensity-$min_reporters/dm_intensity_corrected.txt");
		print OUT_DM_INT_ "reporter\tdm\tintensity\n";
		foreach my $filename (@allfiles)
		{
			if ($filename=~/\.mgf$/i || $filename=~/\.mgf\.gz$/i)
			{

 				if($filename=~/\.mgf\.gz$/i) { `gunzip $filename`; }

				my $count_all_spectra_=0;
				my $count_spectra_=0;
				my $max_max_=0;
				my $filename_=$filename;
				$filename_=~s/\.mgf$//i;
				if(open (IN, "$dir/$filename"))
				{
					#open (OUT_NOT,">$dir/not-selected-d$mz_error-$type-$min_intensity-$min_reporters/$filename");
					#open (OUT_CAL,">$dir/selected-d$mz_error-$type-$min_intensity-$min_reporters/reporter_cal/$filename");
					if(open (OUT,">$dir/selected-d$mz_error-$type-$min_intensity-$min_reporters/$filename"))
					{
						my $pepmass=0;
						my $ms1_intensity="";
						my $title="";
						my $scans="";
						my $charge="";
						my $rt="";
						my @mz=();
						my @intensity=();
						my $header="";
						my $footer="";
						my $started_reading_fragments=0;
						my $done_reading_fragments=0;
						my $points=0;
						my $line="";
						my @top20_title=();
						my @top20_mz=();
						my @top20_int=();
						my @all_mz=();
						my @all_int=();
						my @all_count=();
						my $all_count_max=0;
						while($line=<IN>)
						{
							if ($line=~/^TITLE=(.*)$/)
							{
								$title=$1;
								#if ($title=~/scan=(\d+)/)
								#{
								 #   $scans=$1; 
								    #print $scans,"\n"; 
								#}
								#else
								#{
								#	if ($title=~/^Scan ([0-9\.]+), Time=([0-9\.\-\+edED]+)/)
								#	{
								#		$scans=$1;
								#	
										$rt=$2;
								#	}
								#}
								$scans=(split(/\./,$title))[-2]; 
								#print $scans,"\n"; 
								#<STDIN>;
							}
							if ($line=~/^PEPMASS=([0-9\.\-\+edED]+)\s?([0-9\.\-\+edED]*)\s*$/)
							{
								$pepmass=$1;
								$ms1_intensity=$2;
								#print "pepmass=",$pepmass,"ms1_intensity=",$ms1_intensity,"\n"; 
							}
							if ($line=~/^CHARGE=([0-9\.\-\+]+)\s*$/)
							{
								$charge=$1;
								#print "charge=",$charge,"\n"; 
							}
							if ($line=~/^RTINSECONDS=([0-9\.\-\+]+)\s*$/)
							{
								$rt=$1;
								#print "rt=",$rt,"\n"; 
							}
							if ($line=~/^SCANS=([0-9\.\-\+]+)\s*$/)
							{
								$scans=$1;
								
							}
							if ($line=~/^([0-9\.\+\-edED]+)\s([0-9\.\+\-edED]+)/)
							{
								$started_reading_fragments=1;
								$mz[$points]=$1;
								$intensity[$points]=$2;
								#print $points,"\t",$mz[$points], "\t", $intensity[$points],"\n"; 
								$points++;
							}
							else
							{
								if ($started_reading_fragments==1)
								{
									$done_reading_fragments=1;
								}
							}
							if ($started_reading_fragments==0)
							{
								$header.=$line;
							}
							else
							{
								if ($done_reading_fragments==1)
								{
									$footer.=$line;
									#----------------------------------------------------------#
									my $max=0;
									my $sum=0;
									my @sum=();
									my @reporters_found=();
									my %reporter_dm=();
									my %reporter_intensity=();
									my $reporter_max=$reporters[0];
									my $reporter_count=0;
									foreach my $reporter (@reporters)
									{
										my $max_=0;
										my $dm_at_max=0;
										for(my $i=0;$i<$points;$i++)
										{
											if (abs($reporter-$mz[$i])<$mz_error*$reporter/1e+6)
											{
												$sum+=$intensity[$i];
												$sum[$reporter_count]+=$intensity[$i];
												if ($max<$intensity[$i]) { $max=$intensity[$i]; }
												if ($max_<$intensity[$i]) 
												{ 
													$max_=$intensity[$i]; 
													$dm_at_max=($mz[$i]-$reporter)/($reporter/1e+6);
												}
												if ($min_intensity<$intensity[$i])
												{
													$reporters_found[$reporter_count]=1;
												}
											}
										}
										if (0<$max_)
										{
											#print OUT_DM_INT "$reporter\t$dm_at_max\t$max_\n";
											$reporter_dm{"$reporter"}="$dm_at_max";
											$reporter_intensity{"$reporter"}="$max_";
											if ($reporter_intensity{"$reporter_max"}<=$max_)
											{
												$reporter_max=$reporter;
											}
											if (-$mz_error<$dm_at_max)
											{
												my $k=int(($dm_at_max+$mz_error)/(2*$mz_error)*200);
												$stat{"dm#$reporter#$k"}++;
											}
											my $log2_int=int(log($max_)/log(2));
											$stat{"$filename#$reporter#$log2_int"}++;
											$stat{"$reporter#$log2_int"}++;
										}
										$reporter_count++;
									}
									if (0<$max)
									{
										foreach my $reporter (@reporters)
										{
											if ($reporter_intensity{"$reporter"}=~/\w/)
											{
												my $dm=$reporter_dm{"$reporter"}-$reporter_dm{"$reporter_max"};
												#if(-$mz_error<$dm and $dm<$mz_error)
												#{
												#	print OUT_DM_INT_ qq!$reporter\t$dm\t$reporter_intensity{"$reporter"}\n!;
												#}
											}
										}
									}
									my $reporters_found=0;
									my $reporter_count=0;
									foreach my $reporter (@reporters)
									{
										if ($reporters_found[$reporter_count]==1)
										{
											$reporters_found++;
										}
										$reporter_count++
									}
									$stat{"$filename#$reporters_found"}++;
									$stat{"$reporters_found"}++;
									if ($reporters_found>=$min_reporters)
									{
										if ($max_max<$max) { $max_max=$max; }
										#print OUT $header;
										#for(my $i=0;$i<$points;$i++)
										#{
										#	print OUT "$mz[$i] $intensity[$i]\n";
										#}
										#print OUT $footer;
										#print OUT_CAL $header;
										#for(my $i=0;$i<$points;$i++)
										#{
										#	my $mz=$mz[$i]*(1-$reporter_dm{"$reporter_max"}/1e+6);
										#	print OUT_CAL "$mz $intensity[$i]\n";
										#}
										#print OUT_CAL $footer;
										$header=~s/TITLE=/TITLE=$filename, /m;
										#print OUT_MERGED $header;
										#for(my $i=0;$i<$points;$i++)
										#{
										#	print OUT_MERGED "$mz[$i] $intensity[$i]\n";
										#}
										#print OUT_MERGED $footer;
										#print OUT_CAL_MERGED $header;
										#for(my $i=0;$i<$points;$i++)
										#{
										#	my $mz=$mz[$i]*(1-$reporter_dm{"$reporter_max"}/1e+6);
										#	print OUT_CAL_MERGED "$mz $intensity[$i]\n";
										#}
										#print OUT_CAL_MERGED $footer;
										print OUT_MERGED_TABLE qq!$filename\t$scans\t$rt\t$ms1_intensity!;

										for(my $k=0;$k<$reporter_count;$k++)
										{
											my $sum_=$sum[$k]/(1.0*$sum);
											print OUT_MERGED_TABLE qq!\t$sum_!;
										}
										print OUT_MERGED_TABLE qq!\t$sum\n!;
										$count_spectra_++;
										$count_spectra++;
									}
									else
									{
										#print OUT_NOT $header;
										#for(my $i=0;$i<$points;$i++)
										#{
										#	print OUT_NOT "$mz[$i] $intensity[$i]\n";
										#}
										#print OUT_NOT $footer;
										$header=~s/TITLE=/TITLE=$filename, /m;
										#print OUT_MERGED_NOT $header;
										#for(my $i=0;$i<$points;$i++)
										#{
										#	print OUT_MERGED_NOT "$mz[$i] $intensity[$i]\n";
										#}
										#print OUT_MERGED_NOT $footer;
									}
									#----------------------------------------#
									
									$pepmass="";
									$title="";
									$charge="";
									@mz=();
									@intensity=();
									$header="";
									$footer="";
									$started_reading_fragments=0;
									$done_reading_fragments=0;
									$points=0;
									$count_all_spectra_++;
									$count_all_spectra++;
								}
							}
						}
						close(OUT);
					}
					close(OUT_NOT);
					close(OUT_CAL);
					close(IN);
					
					#if(open (OUT,">$dir/selected-d$mz_error-$type-$min_intensity-$min_reporters/$filename_-intensity-distribution.stat"))
					#{
					#	print OUT "intensity";
					#	foreach my $reporter (@reporters)
					#	{
					#		print OUT "\t$reporter";
					#	}
					#	print OUT "\n";
					#	for(my $i=0;$i<=100;$i++)
					#	{
					#		print OUT $i;
					#		foreach my $reporter (@reporters)
					#		{
					#			if ($stat{"$filename#$reporter#$i"}!~/\w/) { $stat{"$filename#$reporter#$i"}=0; }
					#			print OUT qq!\t$stat{"$filename#$reporter#$i"}!;
					#		}
					#		print OUT "\n";
					#	}
					#	close(OUT);
					#}		
					#if(open (OUT,">$dir/selected-d$mz_error-$type-$min_intensity-$min_reporters/$filename_-reporter-count.stat"))
					#{
					#	print OUT "number_of_reporters\tcount\n";
					#	for(my $i=0;$i<=10;$i++)
					#	{
					#		if ($stat{"$filename#$i"}!~/\w/) { $stat{"$filename#$i"}=0; }
					#		print OUT qq!$i\t$stat{"$filename#$i"}\n!;
					#	}
					#	close(OUT);
					#}
					print "$count_spectra_ spectra selected ($mz_error, $type, $min_intensity, $min_reporters) out of $count_all_spectra_ in \"$filename\".\n";
					#print OUT_COUNT "$filename\t$count_spectra_\t$count_all_spectra_\n";
				} 
				else
				{
					print "Could not open \"$dir/$filename\".\n";
					$error=1;
				}	
			}
		}
		close(OUT_COUNT);
		close(OUT_MERGED);
		close(OUT_CAL_MERGED);
		close(OUT_MERGED_TABLE);
		close(OUT_MERGED_NOT);
		close(OUT_DM_INT);
		close(OUT_DM_INT_);
#	        exit; 
	
		
		#if(open (OUT,">$dir/merged-d$mz_error-$type-$min_intensity-$min_reporters/intensity-distribution.stat"))
		#{
		#	print OUT "intensity";
		#	foreach my $reporter (@reporters)
		#	{
		#		print OUT "\t$reporter";
		#	}
		#	print OUT "\n";
		#	for(my $i=0;$i<=100;$i++)
		#	{
		#		print OUT $i;
		#		foreach my $reporter (@reporters)
		#		{
		#			if ($stat{"$reporter#$i"}!~/\w/) { $stat{"$reporter#$i"}=0; }
		#			print OUT qq!\t$stat{"$reporter#$i"}!;
		#		}
		#		print OUT "\n";
		#	}
		#	close(OUT);
		#}		
		#if(open (OUT,">$dir/merged-d$mz_error-$type-$min_intensity-$min_reporters/reporter-count.stat"))
		#{
		#	print OUT "number_of_reporters\tcount\n";
		#	for(my $i=0;$i<=10;$i++)
		#	{
		#		if ($stat{"$i"}!~/\w/) { $stat{"$i"}=0; }
		#		print OUT qq!$i\t$stat{"$i"}\n!;
		#	}
		#	close(OUT);
		#}			
		#if(open (OUT,">$dir/merged-d$mz_error-$type-$min_intensity-$min_reporters/dm-distribution.stat"))
		#{
		#	print OUT "dm";
		#	foreach my $reporter (@reporters)
		#	{
		#		print OUT "\t$reporter";
		#	}
		#	print OUT "\n";
		#	for(my $i=0;$i<=200;$i++)
		#	{
				
		#		my $dm=$i/200*(2*$mz_error)-$mz_error;
		#		print OUT $dm;
		#		foreach my $reporter (@reporters)
		#		{
		#			if ($stat{"dm#$reporter#$i"}!~/\w/) { $stat{"dm#$reporter#$i"}=0; }
		#			print OUT qq!\t$stat{"dm#$reporter#$i"}!;
		#		}
		#		print OUT "\n";
		#	}
		#	close(OUT);
		#}				
		print "\n$count_spectra spectra selected ($mz_error, $type, $min_intensity, $min_reporters) out of $count_all_spectra in \"merged\".\n";
	}
}
