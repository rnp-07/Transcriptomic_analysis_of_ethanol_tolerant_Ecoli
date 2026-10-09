#This is the script for converting the entrezids into gene symbols which are outputted in enrichKEGG results (in csv files, you will see entrezids instead of symbols because entrezids were used while doing enrichKEGG analysis. This script is useful for converting those entrezids back to symbols).

#Set your working directory by uncommenting the line below and putting the location in the placeholder
#setwd("<write the working directory here>")

library("clusterProfiler")
library("org.EcK12.eg.db")

# Write those entrezids which need to be converted into gene symbols in one .txt file named as "List_to_be_converted.txt" and import it
ids<- read.table("List_to_be_converted.txt", header=FALSE, sep="\t")
gene_ids<- as.character(ids$V1)
gene_ids
class(gene_ids)


# Extract the entire SYMBOL map to a standard Data Frame
all_symbols_df <- toTable(org.EcK12.egSYMBOL)

# View the table(It will have 'gene_id' and 'symbol' columns)
head(all_symbols_df)

gene_ids_df <- data.frame(gene_id = gene_ids)
gene_ids_df

#Now get the symbols by merging the data frame with the all_symbols_df
mapped_data <- merge(gene_ids_df, all_symbols_df, by="gene_id")

#If you want to put all those gene symbols in one list and separate them by ","
final_output <- paste(mapped_data$symbol, collapse = ", ")
final_output

#Save it to a txt file
write.table(final_output, "mapped_gene_symb.txt", row.names=FALSE, quote=F, col.names=FALSE)