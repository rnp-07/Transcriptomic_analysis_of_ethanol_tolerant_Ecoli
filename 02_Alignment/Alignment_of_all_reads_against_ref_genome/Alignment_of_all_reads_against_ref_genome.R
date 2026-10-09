#This is the script for aligning trimmed reads (paired-end reads) of each of the 6 replicates (3 parental and 3 evolved) against reference genome of E coli W (GCF_000258145.1_ASM25814v1)
#Note: This script was run in R 4.5.1 and the package Rsubread was of version 2.22.1.


#Note: if Rsubread is not installed, make sure BiocManager package is installed in your system and try to install Rsubread from BiocManager.
#To do this, uncomment the lines below and run it:
# if (require("BiocManager", quietly=TRUE)==FALSE) {install.packages("BiocManager")}
# BiocManager::install("Rsubread")


#Load Rsubread library
library(Rsubread)

#Build index for reference genome (ensure the reference genome file "GCF_000258145.1_ASM25814v1_genomic.fna" is in working directory)
#Note: reference genome index will be made using basename "ref_gen_index" so it will be used for the parameter "index" inside align() function for all replicates

buildindex("ref_gen_index", "GCF_000258145.1_ASM25814v1_genomic.fna")



#### Alignment of trimmed files (Parent-replicate-1) ####

#Indicate the file path location for trimmed read1 and read2 files for Parent-replicate-1
Pr_1_r1 <-"path/to/your/trimmed/readfile/PARENT_REP1_trimmed_1_val_1.fq.gz"
Pr_1_r2 <-"path/to/your/trimmed/readfile/PARENT_REP1_trimmed_2_val_2.fq.gz"

#Do the alignment (Note: change nthreads according to the number of available CPU cores)
align1 <- align(index="ref_gen_index", readfile1= Pr_1_r1, readfile2= Pr_1_r2, type= "rna", output_file = "Pr_1_alignment.BAM", nthreads= 14)


#### Alignment of trimmed files (Parent-replicate-2) ####

#Indicate the file path locations for trimmed read1 and read2 files for Parent-replicate-2 as done for Parent-replicate-1
Pr_2_r1 <-"path/to/your/trimmed/readfile/PARENT_REP2_trimmed_1_val_1.fq.gz"
Pr_2_r2 <-"path/to/your/trimmed/readfile/PARENT_REP2_trimmed_2_val_2.fq.gz"

#Do the alignment
align2 <- align(index="ref_gen_index", readfile1= Pr_2_r1, readfile2= Pr_2_r2, type= "rna", output_file = "Pr_2_alignment.BAM", nthreads= 14)

#### Alignment of trimmed files (Parent-replicate-3) ####

#Indicate the file path locations for trimmed read1 and read2 files for Parent-replicate-3
Pr_3_r1 <- "path/to/your/trimmed/readfile/PARENT_REP3_trimmed_1_val_1.fq.gz"
Pr_3_r2 <- "path/to/your/trimmed/readfile/PARENT_REP3_trimmed_2_val_2.fq.gz"

#Do the alignment
align3 <- align(index="ref_gen_index", readfile1= Pr_3_r1, readfile2= Pr_3_r2, type= "rna", output_file = "Pr_3_alignment.BAM", nthreads= 14)

#### Alignment of trimmed files (Evolved-Replicate-1) ####

#Indicate the file path locations for trimmed read1 and read2 files for Evolved-replicate-1
Ev_1_r1 <- "path/to/your/trimmed/readfile/EVOLVED_REP1_trimmed_1_val_1.fq.gz"
Ev_1_r2 <- "path/to/your/trimmed/readfile/EVOLVED_REP1_trimmed_2_val_2.fq.gz"

#Do the alignment
align4 <- align(index="ref_gen_index", readfile1= Ev_1_r1, readfile2= Ev_1_r2, type="rna", output_file="Ev_1_alignment.BAM", nthreads=14)

#### Alignment of trimmed files (Evolved-Replicate-2) ####

#Indicate the file path locations for trimmed read1 and read2 files for Evolved-replicate-2
Ev_2_r1 <- "path/to/your/trimmed/readfile/EVOLVED_REP2_trimmed_1_val_1.fq.gz"
Ev_2_r2 <- "path/to/your/trimmed/readfile/EVOLVED_REP2_trimmed_2_val_2.fq.gz"

#Do the alignment
align5 <- align(index="ref_gen_index", readfile1= Ev_2_r1, readfile2= Ev_2_r2, type="rna", output_file="Ev_2_alignment.BAM", nthreads=14)

#### Alignment of trimmed files (Evolved-Replicate-3) ####

#Indicate the file path locations for trimmed read1 and read2 files for Evolved-replicate-3
Ev_3_r1 <- "path/to/your/trimmed/readfile/EVOLVED_REP3_trimmed_1_val_1.fq.gz"
Ev_3_r2 <- "path/to/your/trimmed/readfile/EVOLVED_REP3_trimmed_2_val_2.fq.gz"

#Do the alignment
align6 <- align(index="ref_gen_index", readfile1= Ev_3_r1, readfile2= Ev_3_r2, type="rna", output_file="Ev_3_alignment.BAM", nthreads=14)