#This is the script for doing DESeq2 analysis for all samples excluding Parent-replicate-1
#Note: If you haven't installed DESeq2, please uncomment 4 lines below (from if() to }) and run them (in my case, the analysis was done using DESeq2 v 1.48.1)
# if (require("BiocManager")==FALSE) {
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#}

#Set your working directory by uncommenting the line below and putting the location in the placeholder
#setwd("<paste your working directory here>")
library("DESeq2")

#Load your count matrix
cts<- read.table("genes_plus_pseudogenes_counts_with_adhb_pdc.txt", sep="\t", row.names=1, header=T)
cts <- cts[, -1] # Remove Parent-replicate-1
cts[1:5, 1:5]

#Making a data frame where it will have one column named "condition" and rownames will be sample names
conditions <- c("Parent", "Parent", "Evolved", "Evolved", "Evolved")
sample <- colnames(cts)
coldata<- data.frame(condition = conditions)
rownames(coldata) <- sample
coldata

#Changing datatype to factor
coldata$condition <- factor(coldata$condition)

#Set Parent as the reference level because if you don't do this, "Evolved" will be treated as reference level since it is alphabetically earlier and logfc will be calculated as Parent vs Evolved instead of Evolved vs Parent
coldata$condition <- relevel(coldata$condition, ref = "Parent")

#Create a DESeq DataSet object containing all these information
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

#Do the pre-filtering step where you need to keep only those genes whose counts is >= 10 in at least one sample
keep <- rowSums(counts(dds) >= 10) >= 1
dds <- dds[keep, ]
length(rownames(dds))

#Run DESeq analysis using a function DESeq()
dds <- DESeq(dds)

# See the results
res<- results(dds)
summary(res)

resOrdered <- res[order(res$padj), ]
head(resOrdered)

#Save the results
write.csv(as.data.frame(resOrdered), file = "DESeq_results_wo_Pr1_ss.csv")
# NOTE ON EXCEL WORKBOOK PRESENT IN THE REPOSITORY INSTEAD OF THIS CSV :
# The "DESeq_results_wo_Pr1_ss.xlsx" file provided in the repository was created 
# manually from this csv "DESeq_results_wo_Pr1_ss.csv" generated above. The first sheet "All_genes" basically contains stats for all genes analyzed by DESeq2 while the next sheets
# "sig_genes", "upregulated", "downregulated", "abslogfc_greater_than_log2" and "abslogfc_greater_than_log3" are created after applying padj<0.05, 
# logfc>0 & padj<0.05, logfc<0 & padj<0.05, |logfc|>log2 & padj<0.05 and |logfc|>log3 & padj<0.05 filters respectively, for easier viewing and downstream analysis.

#Do this to write gene symbols for each gene ids
gene_ids_deseq <- rownames(resOrdered)
id_gene_names<- read.csv("ID_gene_names.csv", header=TRUE, sep=",") #This is the file containing gene names for each gene Ids extracted from gff file
id_gene_names[1:5,]

# Add entries for adhB and pdc manually
pdc<- c("gene_pdc", "pdc")
adhB <- c("gene_adhB", "adhB")
id_gene_names<- rbind(id_gene_names, adhB, pdc)
id_gene_names[(nrow(id_gene_names)-1):nrow(id_gene_names),]

#Get the indices of rownames(gene_ids_deseq) in Gene_ID column
ind<- match(gene_ids_deseq, id_gene_names$Gene_ID)
ind[1:10]

#Use the indices to extract the names
deseq_id_with_names<- id_gene_names[ind, ]
deseq_id_with_names[1:10, ]
deseq_id_with_names[deseq_id_with_names$Gene_ID=="gene_adhB"|deseq_id_with_names$Gene_ID=="gene_pdc", ]

#Save it
write.table(deseq_id_with_names, "deseq_id_with_names_wo_Pr1.txt", quote=F, row.names=FALSE, sep="\t")

#See the normalized counts
n_counts<- counts(dds, normalized=TRUE)
head(n_counts)

#Put gene names along with gene ids in n_counts
n_counts_with_names<- merge(n_counts, deseq_id_with_names, by.x="row.names", by.y="Gene_ID")
head(n_counts_with_names)
rownames(n_counts_with_names)<- n_counts_with_names$Row.names
n_counts_with_names<- n_counts_with_names[, c(7, 2:6)] 
write.table(n_counts_with_names, "n_counts_wo_Pr1.txt", sep="\t", quote=F, row.names=TRUE)

# Now, do "Variance Stabilizing Transformation" of original data so that it is useful for visualization
vsd <- vst(dds, blind = FALSE)

# MA plot:
png("MA_plot_wo_Pr1.png", width = 8, height = 6, units = "in", res = 300)
plotMA(res, ylim = c(-3, 3))
dev.off()

# Sample-to-sample distance heatmap:

library("pheatmap")

# Calculate the distances between all of your samples
# We use t(assay(vsd)) to transpose the matrix, 
# because we are calculating distances between samples (columns), not genes (rows) since dist() can calculate distance between only rows not columns
trans<- t(assay(vsd))
sampleDists <- dist(trans)
sampleDists

# We need to convert the distance object into a full matrix
class(sampleDists)
sampleDistMatrix <- as.matrix(sampleDists)

# Create heatmap
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists, filename= "sample-to-sample-distance-heatmap_wo_Pr1.png")

#PCA-plot: Load libraries first
library("ggplot2")
library("ggrepel")

# Add a new column of labels in your vsd object
colnames(vsd)
vsd$labels = colnames(vsd) #Use these sample names for creating extra column named "labels"

pca<- plotPCA(vsd, intgroup = c("condition", "labels"), returnData=TRUE) #returnData=TRUE will save the PCA data in the variable pca 
attr(pca, "percentVar") #You will get percent variance of each axis, PC1 and PC2 from attr() function
variance_per <- round(attr(pca, "percentVar")*100)

#Use PCA data to draw a plot using ggplot
ggplot(pca, aes(PC1, PC2, color = condition, label = labels))+
geom_point(size=3)+
geom_text_repel(show.legend=FALSE)+ #show.legend=FALSE will not print text alongside dots in the legend
xlab(paste("PC1: ", variance_per[1], "% variance", sep=""))+
ylab(paste("PC2: ", variance_per[2], "% variance", sep=""))+
theme_bw()
ggsave("PCA_plot_wo_Pr1.png", width=8, height=6, dpi=300)

#Volcano Plot:
library(EnhancedVolcano)

# Attach gene symbols using the 'id_gene_names' dataframe built earlier
res$Symbol <- id_gene_names$Gene_Name[match(rownames(res), id_gene_names$Gene_ID)]

# Generate the Volcano plot and assign it to a variable 'p'
p <- EnhancedVolcano(res, lab = res$Symbol, x = 'log2FoldChange', y = 'padj', title = 'Volcano plot (Parent vs Evolved)',
pCutoff = 0.05, FCcutoff = 1.0, pointSize = 1.2, colAlpha = 0.4, max.overlaps = Inf,               
legendPosition = 'top', legendIconSize = 4.0)

print(p)

# Save the plots
ggsave("Volcano_Plot_wo_Pr_1.png", plot = p, width = 10, height = 8, dpi = 300)
ggsave("Volcano_Plot_wo_Pr_1.tiff", plot = p, width = 10, height = 8, dpi = 300, compression = "lzw")