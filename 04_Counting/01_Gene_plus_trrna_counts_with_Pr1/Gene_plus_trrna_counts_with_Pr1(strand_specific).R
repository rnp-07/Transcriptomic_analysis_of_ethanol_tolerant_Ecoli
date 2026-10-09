#This is a script for finding the counts of all protein coding genes plus genes for non-coding rna biotypes (trna,rrna, etc.)- so of 4795 genes

#Note: This script was run in R 4.5.1 and the package Rsubread was of version 2.22.1.
#If Rsubread is not installed, make sure BiocManager package is installed in your system and try to install Rsubread from BiocManager.
#To do this, uncomment the lines below and run it:
# if (require("BiocManager", quietly=TRUE)==FALSE) {install.packages("BiocManager")}
# BiocManager::install("Rsubread")

#Load Rsubread library
library(Rsubread)

#### featureCounts analysis of all 6 alignment files ####

bamfiles= c("Pr_1_alignment.BAM", "Pr_2_alignment.BAM", "Pr_3_alignment.BAM", "Ev_1_alignment.BAM", "Ev_2_alignment.BAM", "Ev_3_alignment.BAM") #Give the names of the alignment files here and make sure they are in working directory

gff <- "/path/to/gff file/genomic.gff" #Give the file location to the gff file here (genomic.gff is added here)

fc_new<- featureCounts(files= bamfiles, annot.ext= gff, isGTFAnnotationFile= TRUE, GTF.featureType= "gene", GTF.attrType="ID", isPairedEnd= TRUE, countReadPairs=TRUE,strandSpecific=c(0,2,2,2,2,2), minMQS= 0, nthreads=14) #Pr_1(Parent-replicate-1) was found to be unstranded while rest of the samples were found to be reverse-stranded from RSeQC analysis so O (unstranded) and 2 (reverse-stranded) have been provided for parameter strandSpecific 


#Extracting counts and summary data from the list (fc_new) and saving them as .txt files

fc_new_gene_counts <- fc_new$counts
write.table(fc_new_gene_counts, file= "Gene_plus_trrna_counts_with_Pr1_ss.txt", sep= "\t", quote= FALSE, row.names=TRUE)

#If you want summary stats, you can run the following two lines
fc_new_gene_summ <- fc_new$stat
write.table(fc_new_gene_summ, file= "Gene_plus_trrna_summ_with_Pr1_ss.txt", sep= "\t", quote= FALSE, row.names= TRUE)