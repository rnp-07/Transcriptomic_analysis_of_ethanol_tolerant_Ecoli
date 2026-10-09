#This is a script for combining counts of protein-coding genes and pseudogenes and removing counts of genes for non-coding rna biotypes (trna,rrna, etc.) 

library(biomartr)
gff<-read_gff("genomic.gff") #Replace "genomic.gff" with another name if the reference genome annotation file's name is different and make sure it is in working directory

#Look for the rows with "exon" in the "type" column because non-coding rna biotypes have "exon" as their child metafeature and try to get their Parent IDs (starting with "rna-...") from "attribute" column

exon_parents <- sub(".*Parent=([^;]+).*", "\\1", gff$attribute[gff$type == "exon"])
rna_ids <- unique(exon_parents)

#Extract the 'ID' from every single row in the gff file
all_ids <- sub(".*ID=([^;]+).*", "\\1", gff$attribute)

#Take out only those rows of the "attribute" column where all_ids match rna_ids
rna_rows <- gff$attribute[all_ids %in% rna_ids]

# Now, extract the Parent IDs (the one starting with "gene-...") from those rna_rows
nc_gene_ids <- sub(".*Parent=([^;]+).*", "\\1", rna_rows)
trrna_list <- unique(nc_gene_ids)


#Import two dataframes: one containing counts of genes plus t/rrna genes and another containing pseudogenes counts
gene_trrna_counts<- read.table("Gene_plus_trrna_counts_with_Pr1_ss.txt", header=TRUE, sep="")
pseudogene_counts<- read.table("Pseudogene_counts_with_Pr1_ss.txt", sep="", header=TRUE)

#Combine two dataframes
genes_plus_pseudogenes_plus_trrna<- rbind(gene_trrna_counts, pseudogene_counts)
dim(genes_plus_pseudogenes_plus_trrna)
genes_plus_pseudogenes_plus_trrna[4903:4905, ] #Just checking if these three entries are of last three entries of pseudogenes counts

#Take out list of all gene ids (i.e. gene ids of all protein coding genes+pseudogenes+t/rrna)
all_gene_ids<- rownames(genes_plus_pseudogenes_plus_trrna)
length(all_gene_ids)

#Try to check each element of all_gene_ids in the list of trrna_list
true_false<- all_gene_ids %in% trrna_list
not_true_false<- !true_false #Flip it because we want TRUE in only those positions where there are genes plus pseudogenes not t/rrna

#Now extract the counts using this logical vector (not_true_false)
genes_plus_pseudogenes<- genes_plus_pseudogenes_plus_trrna[not_true_false, ]
dim(genes_plus_pseudogenes) #Should be 4786 rows and 6 columns

#Save it
write.table(genes_plus_pseudogenes, "genes_plus_pseudogenes_counts_wo_trrna_with_Pr1_ss.txt", quote=FALSE, sep="\t", row.names=TRUE)
