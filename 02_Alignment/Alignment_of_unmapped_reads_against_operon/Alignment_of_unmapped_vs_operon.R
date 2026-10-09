#This is the script for aligning unmapped read files against operon sequence. Make sure you have unmapped read 1 and unmapped read2 files for each sample after running .sh script and the sequence file for operon (W105_operon.str) before doing this.
#For this alignment, the original sequence file (W105_operon.str) has been manually converted to .fasta file (operon.fasta) because it is easier to build index of .fasta file.

#Note: This script was run in R 4.5.1 and the package Rsubread was of version 2.22.1.


#Note: if Rsubread is not installed, make sure BiocManager package is installed in your system and try to install Rsubread from BiocManager.
#To do this, uncomment the lines below and run it:
# if (require("BiocManager", quietly=TRUE)==FALSE) {install.packages("BiocManager")}
# BiocManager::install("Rsubread")

#Load Rsubread library
library(Rsubread) 

# Build an index for the operon (index for the operon sequence will be made using basename "operon_index" so it will be used for the parameter "index" inside align() function for all replicates)
buildindex(basename = "operon_index", reference = "operon.fasta")


#### Alignment of unmapped files (Evolved-replicate-1) ####
Ev_1_r1 <-"path/to/your/unmapped read file1 for Ev1/unmapped_R1_Ev1.fastq"
Ev_1_r2 <-"path/to/your/unmapped read file2 for Ev1/unmapped_R2_Ev1.fastq"

#Do the alignment (Note: change nthreads according to the number of available CPU cores)
align1 <- align(index="operon_index", readfile1= Ev_1_r1, readfile2= Ev_1_r2, type= "rna", output_file = "unmapped_vs_operon_Ev1.BAM", nthreads= 14)

#### Alignment of unmapped files (Evolved-replicate-2) #### 
Ev_2_r1 <-"path/to/your/unmapped read file1 for Ev2/unmapped_R1_Ev2.fastq"
Ev_2_r2 <-"path/to/your/unmapped read file2 for Ev2/unmapped_R2_Ev2.fastq"

#Do the alignment
align2 <- align(index="operon_index", readfile1= Ev_2_r1, readfile2= Ev_2_r2, type= "rna", output_file = "unmapped_vs_operon_Ev2.BAM", nthreads= 14)

#### Alignment of unmapped files (Evolved-replicate-3) #### 
Ev_3_r1 <-"path/to/your/unmapped read file1 for Ev3/unmapped_R1_Ev3.fastq"
Ev_3_r2 <-"path/to/your/unmapped read file2 for Ev3/unmapped_R2_Ev3.fastq"

#Do the alignment
align3 <- align(index="operon_index", readfile1= Ev_3_r1, readfile2= Ev_3_r2, type= "rna", output_file = "unmapped_vs_operon_Ev3.BAM", nthreads= 14)

#### Alignment of unmapped files (Parent-replicate-1) #### 
Pr_1_r1 <-"path/to/your/unmapped read file1 for Pr1/unmapped_R1_Pr1.fastq"
Pr_1_r2 <-"path/to/your/unmapped read file2 for Pr1/unmapped_R2_Pr1.fastq"

#Do the alignment
align4 <- align(index="operon_index", readfile1= Pr_1_r1, readfile2= Pr_1_r2, type= "rna", output_file = "unmapped_vs_operon_Pr1.BAM", nthreads= 14)

#### Alignment of unmapped files (Parent-replicate-2) #### 
Pr_2_r1 <-"path/to/your/unmapped read file1 for Pr2/unmapped_R1_Pr2.fastq"
Pr_2_r2 <-"path/to/your/unmapped read file2 for Pr2/unmapped_R2_Pr2.fastq"

#Do the alignment
align5 <- align(index="operon_index", readfile1= Pr_2_r1, readfile2= Pr_2_r2, type= "rna", output_file = "unmapped_vs_operon_Pr2.BAM", nthreads= 14)

#### Alignment of unmapped files (Parent-replicate-3) #### 
Pr_3_r1 <-"path/to/your/unmapped read file1 for Pr3/unmapped_R1_Pr3.fastq"
Pr_3_r2 <-"path/to/your/unmapped read file2 for Pr3/unmapped_R2_Pr3.fastq"

#Do the alignment
align6 <- align(index="operon_index", readfile1= Pr_3_r1, readfile2= Pr_3_r2, type= "rna", output_file = "unmapped_vs_operon_Pr3.BAM", nthreads= 14)
