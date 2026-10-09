library(biomartr)
gff<-read_gff("genomic.gff")

#Look for the rows with "exon" in the "type" column because non-coding rna biotypes have "exon" as their child metafeature and try to get Parent IDs (starting with "rna-...") from them
exon_parents <- sub(".*Parent=([^;]+).*", "\\1", gff$attribute[gff$type == "exon"])
rna_ids <- unique(exon_parents)

#Extract the 'ID' from every single row in the GFF
all_ids <- sub(".*ID=([^;]+).*", "\\1", gff$attribute)

#Take out only those rows of the attribute column where all_ids match rna_ids
rna_rows <- gff$attribute[all_ids %in% rna_ids]

# Now, extract the Parent IDs (the one starting with "gene-...") from those rna_rows
nc_gene_ids <- sub(".*Parent=([^;]+).*", "\\1", rna_rows)
nc_gene_ids <- unique(nc_gene_ids)
