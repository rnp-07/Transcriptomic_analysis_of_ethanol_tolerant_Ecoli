#This is the script for finding gene length of all the genes mentioned in gff file
#The tool used was "biomartr" v 1.0.7 in R v 4.5.1.

#If you haven't installed biomartr, please run the following lines:
if (!require("biomartr", quietly=TRUE) {
install.packages("biomartr")
}

library("biomartr")
gff<- read_gff("genomic.gff") #Make sure this gff file is in working directory or provide its location here

#Extract start coordinate, end coordinate and attribute columns of only those rows where "type" is either "gene" or "pseudogene"
gff_subset<- gff[gff$type=="gene"|gff$type=="pseudogene", c(4,5,9)]
gff_subset[1:3, ]

#Extract IDs of all those genes and pseudogenes
gff_subset_ids<- sub(".*ID=([^;]+).*", "\\1", gff_subset$attribute)
gff_subset_ids[1:5]

#Put those ids as the first column, and length as the second column in the new dataframe
Ids_with_length<- data.frame(IDs=gff_subset_ids, Length= gff_subset$end-gff_subset$start+1)

