#This is the script for plotting the normalized counts distribution curve

#Import the table containinig normalized counts of all genes in all samples including Parent-replicate-1
dist<- read.table("normalized_counts_w_Pr1.txt", sep="\t", header=TRUE, row.names=1 )
dist[1:5, ]

library(ggplot2)
library(reshape2)
library(scales)

#Melt the data (Wide to Long)
long_data <- melt(dist, id.vars = c("Gene_Name"))
colnames(long_data) <- c("Gene_Name", "Sample", "Count")

#Assign the plot to a variable 'p'
p <- ggplot() +
# Density Plot
geom_density(data = long_data, aes(x = Count + 1, fill = Sample), alpha = 0.3) +

# Scale
scale_x_log10(labels = scales::label_comma()) + 
  
# Themes
theme_minimal() +
theme(plot.title = element_text(size = 26, face = "bold"),
axis.title = element_text(size = 18, face = "bold"),
axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
axis.text.y = element_text(size = 18),
strip.text = element_text(size = 20, face = "bold"),
legend.title = element_text(size = 18, face = "bold"),
legend.text = element_text(size = 16),
legend.position = "bottom") +
        
# Labels and Facets
labs(title = "Gene counts distribution", x = "Normalized Counts + 1 (Log Scale)", y = "Density") +
facet_wrap(~Sample)

#Print the plot
print(p)

#Save it
ggsave("n_counts_distr_curve.png", plot = p, width = 14, height = 10, dpi = 300)