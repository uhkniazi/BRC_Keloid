# Title
Identifying what goes wrong during wound healing to cause excessive scars known as keloids.

**DB ID** Project.id=1, Data.id=2  

[Evernote Link - Private](https://www.evernote.com/shard/s288/nl/38698211/0a8bc5dc-e07d-40e4-9077-725c94e97bcd?title=00%20Keloid%20RNA-Seq)
## Comments
300 fastq files processed for quality checks. All seems fine, some trimming on the 3' side will be required on the first library.  

## Scripts
1. rna_seq_qa.R
  * Process a batch of FASTQ files and produces html output reports. Adds results to Project.File table.  
2. rna_seq_qa_2.R 
  * Similar to previous script but creates one large object with all ~ 260 FASTQ files from the first experiment, and creates one pdf document with read quality plots. **NOTE** memory and time intensive.  
3. merge_fastq.R
  * Merge a batch of FASTQ files according to sample ids and paired end direction.  
4. write_trimmomatic_script.R
  * Writes a script file in the AutoScripts folder that can be used on the HPC to run trimmomatic on the fastq files. **NOTE** create the correct directories for input and output in the folder where the script will be executed.  
5. write_trimmomatic_array_job.R
  * similar to previous script but creates an array job script and a parameter file in the AutoScripts folder.  
6. rna_seq_qa_3.R  
  * Similar to previous scripts 1 and 2, however uses the updated version of CFastqQuality class. Saves the object information in the Projects.MetaFile table in the database and the object in the appropriate folder. Outputs the relevant figures. **NOTE** memory and time intensive.  
7. write_hisat2_array_job.R  
  * array job for hpc to perform alignment using hisat2.  
8. write_samtools_array_job.R  
  * array job for cleaning up sam file after alignment and producing bam file.  
9. bam_files_qa.R  
  * uses the CBamScaffold class to look at the 22 chromosomes of each bam file (depending on the sequencing run) and produces summary plots for the alignment coverage quality.  
10. counts_from_bams.R  
  * uses all the bam files to and Txdb object to count over laps over the genome and save the results.  
11. clustering_rna_seq_counts.R  
  * loads the sample information (covariates etc) and counts matrix from the database. Produces PCA plots for the count matrix, coloured by covariate of choice, using various transformations of the data e.g. scaling, log, rlog, normalised. Plots the average gene expression density distribution. Merges the technical replicates, normalises and repeats the plots.  
  

  
  

