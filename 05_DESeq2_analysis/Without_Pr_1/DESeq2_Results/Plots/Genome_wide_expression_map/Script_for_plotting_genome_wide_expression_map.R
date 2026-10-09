#This is the script for making genome wide expression map
#First run DESeq2 on all samples excluding Parent-replicate-1

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

#Load a file named "all_genes_information_from_gff.csv" where there are several columns of information (Gene_IDs, Gene_Names, seqid, type, start and strand) for all genes+pseudogenes extracted from gff file

gff_inf<- read.csv("all_genes_information_from_gff.csv", sep=",", header=T)
gff_inf[1:5, ]
gff_clean<- gff_inf[, c(1,2,3,5)]
gff_clean[1:5, ]

#Make a dataframe of DESeq2 results and add a column "Gene_IDs" in it
deseq_df<- as.data.frame(resOrdered)
deseq_df<- data.frame(Gene_IDs= rownames(deseq_df), deseq_df)
deseq_df[1:5, ]

#Now, merge two dataframes; gff_clean and deseq_df by "Gene_IDs" and keep all the rows of gff_clean because there were some genes which were not analyzed by DESeq2 which we want to show in the plot
plot_df<- merge(gff_clean, deseq_df, by="Gene_IDs", all.x=TRUE)
nrow(plot_df)
plot_df[1:10, ]

#Just extract the required columns for the plot
req_plot <- plot_df[, c(1,3,4,6,10)]
req_plot[1:10, ]


#Now make the plot using req_plot dataframe
#Load the required libraries
library(ggplot2)
library(dplyr)

# load the t/rrna gene list (the one with just gene IDs of t/rrna genes)
trna_table <- read.table("trna_rrna_gene_list.txt", header=TRUE, sep="\t")
trna_ids <- trna_table[,1]


# Create two new columns: 1) "Plot_LogFC" that will have same logfc values as in log2FoldChange column but NA will be changed to 0 because we want those genes to appear on the plot and 
# 2) "Category" where the gene will be classified as either "t/rRNAs+Filtered out" (Filtered out means those genes which were filtered out during Pre-filtering step during DESeq2 analysis i.e. genes having less than 10 raw counts in all samples) or "Upregulated" or "Downregulated" or "Non-Significant"

req_plot_colored <- req_plot %>%
mutate(Plot_LogFC = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
Category = case_when(Gene_IDs %in% trna_ids | is.na(padj) ~ "t/rRNAs+Filtered out", padj < 0.05 & log2FoldChange > 0 ~ "Upregulated", padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",TRUE ~ "Non-Significant"))

# Arrange the data based on categories mentioned in "Category" column so that "Non-Significant" comes first, "t/rRNAs+Filtered out" comes second, "Downregulated" comes third and "Upregulated" comes fourth
# This makes the points of "Upregulated" to come on the foreground while "Non-Significant" to be at the background

req_plot_colored <- req_plot_colored %>% arrange(factor(Category, levels = c("Non-Significant", "t/rRNAs+Filtered out", "Downregulated", "Upregulated")))

# Pick 4 colors for 4 categories
p<- ggplot(req_plot_colored, aes(x = start, y = Plot_LogFC)) + geom_point(aes(color = Category, alpha = Category), size = 1.5) +
scale_color_manual(values = c("Upregulated" = "forestgreen", "Downregulated" = "firebrick", "Non-Significant" = "grey70", "t/rRNAs+Filtered out" = "gold")) +

scale_alpha_manual(values = c("Upregulated" = 1, "Downregulated" = 1, "Non-Significant" = 0.6, "t/rRNAs+Filtered out" = 0.6)) +
facet_wrap(~seqid, ncol = 1, scales = "free_x") + theme_bw() + labs(title = "Genome-Wide Expression Map", x = "Genomic Coordinate (bp)", y = "Log2 Fold Change") +
geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3)

ggsave("genome_wide_expression_map.png", plot=p, width=10, height=8, dpi=300)