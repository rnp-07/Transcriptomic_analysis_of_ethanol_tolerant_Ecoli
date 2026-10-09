#This is the script for doing enrichKEGG analysis

#Set your working directory by uncommenting the line below and putting the location in the placeholder
#setwd("<place the location here>")

library("clusterProfiler")
library("org.EcK12.eg.db")

#Get background list of genes (with gene names only)- this is the universe
back_list<- read.table("Background_list.txt", header=TRUE, sep= "\t") 
nrow(back_list)
background <- back_list$Gene_Name
length(background)

#Map gene symbols of our strain to entrezid of K12 cause symbols are not recognized as keyType in enrichKEGG function
back_k12<-bitr(geneID=background, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.EcK12.eg.db)

head(back_k12)
nrow(back_k12)

#Get the list of significant genes (with gene names only)- this is present in each folder named after threshold of logfc (e.g. "no_logfc_threshold", "abs_logfc_greater_than_log2", "abs_logfc_greater_than_log3")

sig_list<- read.table("Significant_genes_list - Copy.txt", header=TRUE, sep= "\t") 
sig_gene_names <- sig_list$Gene_Name
length(sig_gene_names)

#Map gene symbols to entrezid like before
sig_k12<-bitr(geneID=sig_gene_names, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.EcK12.eg.db)

head(sig_k12)
nrow(sig_k12)


#Do enrichKEGG analysis
#Remember keyType should be "ncbi-geneid" not "kegg" as entrezid is recognized as "ncbi-geneid" here. 
keg <- enrichKEGG(gene= sig_k12$ENTREZID, universe=back_k12$ENTREZID, organism= 'eco', keyType='ncbi-geneid', pvalueCutoff=0.05, pAdjustMethod = "BH", minGSSize = 3)

head(keg)
write.csv(keg, "kegg_k12.csv", row.names=FALSE)

#Barplot
bp<- barplot(keg)
ggplot2::ggsave("Barplot.png", bp, width=12, height=12)

#Dot plot
dp <- dotplot(keg)
ggplot2::ggsave("Dotplot.png", dp, width = 10, height = 12)

#Plot emapplot to see which pathways are related
library(enrichplot)
keg2 <- pairwise_termsim(keg) 
emapplot(keg2, label_format = 10)

#Pathview Plots for each enriched pathway:
library("pathview")

#Load the internal dataset
data(gene.idtype.list)

# View the full list of supported Gene IDs
print(gene.idtype.list)

#Since symbols are supported, load the full table containing gene symbols and their corresponding logfc values obtained from DESeq2 analysis
name_logfc<- read.table("gene_names_with_logfc.txt", sep="\t", header=T)
name_logfc[1:5, ]

gene_data_vector2<- name_logfc$log2FoldChange
names(gene_data_vector2)<- name_logfc$Gene_Name
gene_data_vector2[1:10]

# Now run Pathview

pathway_id <- keg@result$ID[2] #This is just for 2nd enriched pathway- if you want for other pathway, simply change the number 
pathview(gene.data=gene_data_vector2, pathway.id=pathway_id, species="eco", gene.idtype="symbol", limit=list(gene=2,cpd=1))