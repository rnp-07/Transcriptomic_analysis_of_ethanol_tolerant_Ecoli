#This is the script for quantification of operon genes (pdc and adhB)

#Note: This script was run in R 4.5.1 and the package Rsubread was of version 2.22.1.
#If Rsubread is not installed, make sure BiocManager package is installed in your system and try to install Rsubread from BiocManager.
#To do this, uncomment the lines below and run it:
# if (require("BiocManager", quietly=TRUE)==FALSE) {install.packages("BiocManager")}
# BiocManager::install("Rsubread")

#Load Rsubread library
library(Rsubread)

#### featureCounts analysis of all 6 alignment files ####

#Make sure the alignment files (.BAM) which were obtained after aligning unmapped reads against operon in the previous step are in your working directory
bam_files <- c("unmapped_vs_operon_Pr1.BAM", "unmapped_vs_operon_Pr2.BAM", "unmapped_vs_operon_Pr3.BAM", "unmapped_vs_operon_Ev1.BAM", "unmapped_vs_operon_Ev2.BAM", "unmapped_vs_operon_Ev3.BAM") 

#Make sure operon.gff is in your working directory. The gff provided here is a minimal, custom annotation file which includes only the essential coordinate, feature type and strand orientation required by featureCounts to quantify the adhB and pdc genes in this operon.
#Also change nthreads according to the number of available cores in your CPU

operon_counts <- featureCounts(bam_files, annot.ext= "operon.gff", isGTFAnnotationFile = TRUE, GTF.featureType = "gene", GTF.attrType= "ID", isPairedEnd= TRUE, countReadPairs=TRUE, strandSpecific=c(0,2,2,2,2,2), minMQS=0, nthreads=14 ) 

#Save the counts as a .txt file
adhb_pdc_counts<- operon_counts$counts
write.table(adhb_pdc_counts, file = "adhb_pdc_counts(unmapped)_ss.txt", sep = "\t", quote= FALSE, row.names = TRUE)

