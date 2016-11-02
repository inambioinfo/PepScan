PepScan version 0.1

PepScan pipeline is a fully automated and modular software package designed for protein quantification and the identification of peptides from mass spectra data. It works on LSF job scheduler. 

Usage: perl PepScan.pl data_dir output_dir step_number


Example: perl PepScan.pl /gscuser/scao/gc2521/dinglab/cptac/test_pepscan /gscuser/scao/gc2521/dinglab/cptac/test_pepscan/results 3

step_number: Run this pipeline step by step. (Running all steps if step number is 0)

[1] Use Quilts for the annotation of snvs and in-house script to do annotation of indels. 

[2] Create fasta DB for running MSGF+.

[3] Run MSGF+ for one sample.

[4] Run MSGF+ for multiple samples.

[5] Convert mzid to tsv file.

[6] Extract reporter ion intensity.

[7] Estimate protein expression.

[8] Report SNVs supported by mass spectra data.  
 
If you have any questions and suggestions for PepScan pipeline, please contact scao@wustl.edu. 
