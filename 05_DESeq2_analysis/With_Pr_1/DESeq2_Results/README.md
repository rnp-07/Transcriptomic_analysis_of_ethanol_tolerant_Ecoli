# DESeq2 Analysis (Including Parent-Replicate-1)

This folder contains files related to the DESeq2 analysis of all samples, including Parent-replicate-1.

### Inputs
* **`genes_plus_pseudogenes_counts_with_adhb_pdc.txt`**: The primary input file containing raw counts of all protein-coding genes and pseudogenes obtained in the previous step, with *adhB* and *pdc* counts appended at the end.
* **`ID_gene_names.csv`**: This is a list of gene ids and their respective gene symbols which were extracted from the GFF file using similar methodology as applied for RNA biotype IDs during the counting step. This is used to add gene symbols after getting DESeq2 results because DESeq2 results will just give gene ids as gene identifier.

### Outputs & Subdirectories
* **`DESeq_results_with_Pr1.csv`**: Contains the final DESeq2 differential expression results for all analyzed samples.
* **`normalized_counts_w_Pr1.txt`**: Contains the normalized counts calculated by DESeq2 for each gene across all 6 samples, mapped alongside their gene symbols.
* **`Plots/`**: A folder containing various visualization plots (MA plot, PCA plot, and sample-to-sample distance heatmap) obtained after DESeq2 analysis.
* **`Housekeeping_genes/`**: This folder contains an R script (`SCRIPT_for_bar_plot_of_housekeeping_genes_counts.R`) which imports the normalized counts for 6 specific housekeeping genes (*gyrB*, *rpoD*, *gapA*, *arcA*, *rpoB*, *dnaN*) from a count matrix file (`housekeeping_counts.txt`). The script then generates a bar plot (`housekeeping_genes_plot.png`) that allows to see if there was a drastic drop/spike in gene expression in one of these stable housekeeping genes across the samples. Sudden change in one or more of these genes' expression indicates potential RNA degradation or contamination (technical failure) rather than a true biological difference. 
* **`normalized_counts_distribution_curve/`**: Contains the R script (`Script_for_n_counts_distributn_curve.R`) which imports the complete normalized counts matrix (`normalized_counts_w_Pr1.txt`) generated after DESeq2 analysis and generates a density plot (`n_counts_distr_curve.png`). This plot helps to see the relative proportion of genes at a given normalized counts level (counts+1 precisely) where the X-axis is logarithmically transformed (meaning normalized counts+1 increase by 10 fold with each grid spacing in the X-axis i.e. 1, 10, 100, 1000, etc.).
