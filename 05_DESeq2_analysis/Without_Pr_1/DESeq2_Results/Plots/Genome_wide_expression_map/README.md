# Genome-Wide Expression Map

This directory contains the scripts and required datasets to generate a genome-wide map of differential expression (including the chromosome and plasmids) based on the DESeq2 results (excluding the Parent-replicate-1 sample).

### Inputs
* **`genes_plus_pseudogenes_counts_with_adhb_pdc.txt`**: The main count matrix containing the raw counts of all protein-coding genes and pseudogenes (including *adhB* and *pdc*), which is required to rerun the DESeq2 analysis before plotting.
* **`all_genes_information_from_gff.csv`**: A reference file containing the genomic coordinates and metadata for all genes, extracted directly from the GFF file. This file contains 6 columns, namely `Gene_IDs`, `Gene_Names`, `seqid` (sequence ID), `type` (feature type), `start` (start coordinate), and `strand` (strand orientation).

    This file is required to accurately map the log2 fold-changes to their actual physical locations on the genome.
* **`trna_rrna_gene_list.txt`**: A simple list containing only the Gene IDs of all non-coding RNA biotypes (tRNA, rRNA, etc.). This is used by the script to categorize these specific biotypes on the expression map.

### Script
* **`Script_for_plotting_genome_wide_expression_map.R`**: The complete R script required to generate the final genome-wide expression map.

### Output
* **`genome_wide_expression_map.png`**: A faceted dot plot visualizing the log2 fold-change of every feature mapped across the entire genome. Genes are color-coded into four specific categories for visual clarity:
  1. **Upregulated:** $\log_2\text{FC} > 0$ and $p_\text{adj} < 0.05$ 
  2. **Downregulated:** $\log_2\text{FC} < 0$ and $p_\text{adj} < 0.05$ 
  3. **Non-Significant:** Analyzed genes that did not meet the $p_\text{adj} < 0.05$ threshold.
  4. **t/rRNAs + Filtered Out:** Non-coding RNA biotypes (tRNA, rRNA, etc.), plus those genes (or pseudogenes) that were dropped during the DESeq2 pre-filtering step (the ones that had less than 10 raw counts across all samples).

