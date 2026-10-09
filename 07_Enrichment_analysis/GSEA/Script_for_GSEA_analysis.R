#This is the script for doing GSEA (both GSEA-GO and GSEA-KEGG) analysis
#Before doing it, do DESeq2 analysis because it will be easy to import the list of all genes with logfc values required for this analysis.


#Set your working directory to the location where there are required files for doing DESeq2 analysis
#setwd("<place the correct location here after uncommenting this line>")

library("DESeq2")

#Load the count matrix containing counts of protein coding genes+pseudogenes+adhB+pdc
cts<- read.table("genes_plus_pseudogenes_counts_with_adhb_pdc.txt", sep="\t", row.names=1, header=T)
cts[1:5, ]
cts<- cts[, -1] #Remove parent-replicate 1

#Number of relevant genes
nrow(cts)

#Making a data frame where it will have one column named "condition" and rownames will be sample names
conditions <- c("Parent", "Parent", "Evolved", "Evolved", "Evolved")
sample <- colnames(cts)
coldata<- data.frame(condition = conditions)
rownames(coldata) <- sample
coldata

#Changing datatype to factor
coldata$condition <- factor(coldata$condition)

#Apparently I have to do this otherwise R will treat as "Evolved" as reference since it is alphabetically earlier
coldata$condition <- relevel(coldata$condition, ref = "Parent")

#I have to create a DESeq object containing all these information
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

#I have to do Pre-filtering step where I need to keep only those genes whose counts is >= 10
keep <- rowSums(counts(dds) >= 10) >= 1
dds <- dds[keep, ]
length(rownames(dds))

#Then I have to run DESeq analysis using a function DESeq()
dds <- DESeq(dds)

# See the results
res<- results(dds)
summary(res)

resOrdered <- res[order(res$padj), ]
head(resOrdered)

#Get gene symbols for the gene ids
gene_ids_deseq <- rownames(resOrdered)
id_gene_names<- read.csv("ID_gene_names.csv", header=TRUE, sep=",") #This is a list of all gene ids with their corresponding gene names
id_gene_names[1:5,]

#Get the indices of rownames(gene_ids_deseq) in Gene_ID column
ind<- match(gene_ids_deseq, id_gene_names$Gene_ID)
ind[1:10]

#use the indices to extract the names
deseq_id_with_names<- id_gene_names[ind, ]
deseq_id_with_names[1:10, ]
nrow(deseq_id_with_names)

logfc_with_names<- resOrdered[,2]
logfc_with_names[c(100, 284) ]
names(logfc_with_names)<- deseq_id_with_names[,2]
names(logfc_with_names)[c(100,284)]<- c("adhB", "pdc")

#Check if NA is present in the names
sum(is.na(names(logfc_with_names)))

#To see if any name is duplicated
name_counts<- table(names(logfc_with_names))
name_counts[name_counts>1]

#Make an average of duplicated symbols
averaged_logfc <- tapply(logfc_with_names, names(logfc_with_names), mean, na.rm = TRUE)
averaged_logfc[1:10]

#Sort it in decreasing order
logfc_ordered<- averaged_logfc[rev(order(averaged_logfc))]
logfc_ordered[1:8]
final_gene_list <- as.numeric(logfc_ordered)
names(final_gene_list) <- names(logfc_ordered)

#Load the required libraries
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.EcK12.eg.db)

# Run Gene Ontology (GO) GSEA
 gse_go <- gseGO(geneList     = final_gene_list,
                OrgDb        = org.EcK12.eg.db,
                keyType      = "SYMBOL",
                ont          = "CC",     
                minGSSize    = 3,        
                maxGSSize    = 800,       
                pvalueCutoff = 0.05,
                verbose      = TRUE,
                pAdjustMethod = "BH")
# Check the results table
head(as.data.frame(gse_go))

write.csv(as.data.frame(gse_go), "GSEA_GO_CC_padj0.05.csv")

gseaplot2(gse_go, geneSetID = "GO:0046930", title = "Enrichment Plot: Pore Complex", base_size= 16, pvalue_table = TRUE)

#For GSEA KEGG analysis

#First convert the symbols into entrezids because symbols are not recognized
sym_to_ent<- bitr(names(final_gene_list), fromType="SYMBOL", toType="ENTREZID", OrgDb= org.EcK12.eg.db)
sym_to_ent[1:5, ]

kegg_list<- final_gene_list[names(final_gene_list) %in% sym_to_ent$SYMBOL]
names(kegg_list)<- sym_to_ent$ENTREZID[match(names(kegg_list), sym_to_ent$SYMBOL)] 

gse_kegg<- gseKEGG(geneList     = kegg_list,
                organism        = "eco",
                keyType      = "ncbi-geneid",     
                minGSSize    = 3,        
                maxGSSize    = 800,       
                pvalueCutoff = 0.05,
                verbose      = TRUE,
                pAdjustMethod = "BH")
head(as.data.frame(gse_kegg))

write.csv(as.data.frame(gse_kegg), "GSEA_KEGG_0.05.csv")
gseaplot2(gse_kegg, geneSetID="eco02040", title="Enrichment Plot: Flagellar Assembly", base_size=16, pvalue_table=TRUE)
gseaplot2(gse_kegg, geneSetID="eco00440", title="Enrichment Plot: Phosphonate & Phosphinate Metabolism", base_size = 16, pvalue_table=TRUE)
