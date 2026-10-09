#This is the script for making bar plot for housekeeping genes' normalized counts

library(ggplot2)
library(tidyr)
library(dplyr)

#Load the counts file
counts_df <- read.delim("housekeeping_counts.txt", header = TRUE, sep = "\t")
print(colnames(counts_df))

#Reshape the data from "wide" to "long" so that there will be three columns, 1) there will be gene names, 2)  there will be normalized counts of each gene in all samples and 3) there will be sample names.

counts_long <- counts_df %>%
pivot_longer(cols = -Gene, names_to = "Samples",values_to = "Normalized_Counts")

#Define the exact order for the samples
sample_order <- c("Parent.replicate.1", "Parent.replicate.2", "Parent.replicate.3", "Evolved.replicate.1", "Evolved.replicate.2", "Evolved.replicate.3")

# Change the type of data inside "Samples" column into factor datatype so that the first heirarchical position is of "Parent.replicate.1", second  is of "Parent.replicate.2" and so on (this is important for which bar will appear first, which bar will appear second, etc.)

counts_long_ordered <- counts_long %>%
mutate(Samples = factor(Samples, levels = sample_order))

#Create the plot
p <- ggplot(counts_long_ordered, aes(x = Samples, y = Normalized_Counts, fill = Samples)) +
  geom_col() + 
  facet_wrap(~ Gene, scales = "free_y") +
  labs(title = "Normalized counts of 6 housekeeping genes across samples", x = "Samples", y = "Normalized Counts") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), axis.text.x = element_blank(), axis.ticks.x = element_blank())

#Print the plot and save it
print(p)
ggsave("housekeeping_genes_plot.png", plot = p, device = "png", width = 10, height = 8, dpi = 300)