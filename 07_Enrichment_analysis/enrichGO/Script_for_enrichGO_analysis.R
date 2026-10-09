#This is the script for doing enrichGO analysis and generating the required plots including cnet plot and heatmap plot.

#Set your working directory by uncommenting the line below and putting the location in the placeholder
#setwd("<put the location here>")

#Get background list (or universe) of genes (with gene names only)- This is the list of all genes analyzed by DESeq2
background<- read.table("Background_list.txt", header=TRUE)
nrow(background)
colnames(background)
wfl<- grepl("^WFL", background$Gene_Name)
background_gnames_only<- background[!wfl, ] #Just extract the ones which have got gene symbols (not the ones starting with WFL...)
length(background_gnames_only)

#Get significant genes list (with gene names only)- this is present in each folder named after threshold of logfc (e.g. "no_logfc_threshold", "abs_logfc_greater_than_log2", "abs_logfc_greater_than_log3")
sig<- read.table("Significant_genes_list - Copy.txt", header=TRUE) #Move the file to the working directory before doing this
nrow(sig)
colnames(sig)
wfl<- grepl("^WFL", sig$Gene_Name)
sig_gnames_only<- sig[!wfl, ] #Just extract the ones which have got gene symbols (not the ones starting with WFL...) 
length(sig_gnames_only)

#Check whether background and sig genes are vector or not
class(background_gnames_only)
class(sig_gnames_only)

#Do enrichGO analysis
library("clusterProfiler")
library("org.EcK12.eg.db")

ego<-enrichGO(gene=sig_gnames_only, universe=background_gnames_only, OrgDb=org.EcK12.eg.db, keyType="SYMBOL", ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2, minGSSize=3) #Change ont to others if you want to see enriched "CC" terms or "MF" terms

head(ego)
dim(ego)

ego_df <- as.data.frame(ego)
write.csv(ego_df, "enrichGO_results.csv", row.names = FALSE)

#Barplot
bp<- barplot(ego, showCategory=31)
library(ggplot2)
ggsave("Barplot.png", bp, width=12, height=12)

#Dot plot
dp <- dotplot(ego, showCategory=31)
ggsave("Dotplot.png", dp, width = 10, height = 12)

#Plot emapplot to see which pathways are related
library(enrichplot)
ego2 <- pairwise_termsim(ego) 
emapplot(ego2, label_format = 10)


#For plotting cnet plot
name_logfc<- read.table("gene_names_with_logfc.txt", sep="\t", header=T)
name_logfc[1:5, ]
foldchanges<- name_logfc$log2FoldChange
names(foldchanges)<- name_logfc$Gene_Name

foldchanges[1:5]
options(ggrepel.max.overlaps = Inf)
p<- cnetplot(ego, foldChange= foldchanges, showCategory=nrow(ego))
p <- p + scale_color_gradientn(colors=c("red","white","green"),limits=c(-2,2),oob=scales::squish,name="Log2 FC")
ggsave("cnetplot.pdf", width = 30, height = 20)


# For plotting heatmap
p <- heatplot(ego, foldChange = foldchanges, showCategory = 31)

# Add black borders to the boxes so that individual tiles look distinct
p$layers[[1]]$aes_params$colour <- "black"
p$layers[[1]]$aes_params$size   <- 0.3

# Colors
p<-p+scale_fill_gradient2(name="logfc",low="green",mid="white",high="red",midpoint=0,limits=c(-3,3),oob=scales::squish)

# Add a vertical Line arising from each gene name upto the top
p <- p + theme_minimal() + theme(panel.grid.major.x = element_line(color = "grey50", linetype = "solid", size = 0.3), panel.grid.major.y = element_blank(),

# Ensure gene names align perfectly with the lines
axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
axis.text.y = element_text(size = 10, color = "black"),
    
# Remove the grey background so lines pop out
panel.background = element_rect(fill = "white", color = NA))

# Save with enough width to see the lines clearly
ggsave("heatmap_with_logfc.png", p, width = 35, height = 12)