# Over-representation Analysis (enrichGO)

This directory contains the script, input files and the outputs of over-representation analysis (ORA) using the Gene Ontology: Biological Process (GO-BP) database using `enrichGO` function from `clusterProfiler` package. The analysis was performed independently across three different fold-change thresholds: i) no logfc threshold ii) $|\log_2\text{fc}|>log_22$ and iii) $|\log_2\text{fc}|>log_23$


### Input Files
* **`Background_list.txt`**: The universe/background list of all 4,641 genes analyzed by DESeq2. 

* **`gene_names_with_logfc.txt`**: A reference file containing the log2 fold-change values for all 4,641 genes. This was just taken out from DESeq2 results. This file was strictly required to map fold-change color gradients onto the `cnetplot` and `heatmap plot` visualizations.

### Script
* **`Script_for_enrichGO_analysis.R`**: The R script used to filter the gene lists, execute the `enrichGO` function, and generate the resulting CSV data and pathway plots.

### Subdirectories (threshold related)
There are 3 subdirectories here related to 3 thresholds which were taken for this analysis: "`no_logfc_threshold/` (all significant genes i.e. $p_\text{adj} < 0.05$)", "`abs_logfc_greater_than_log2/`" (i.e. significant genes satisfying $|\log_2\text{fc}|>log2$ threshold) and "`abs_logfc_greater_than_log3/`"(i.e. significant genes satisfying $|\log_2\text{fc}|>log3$ threshold). For each threshold, the list of complying genes is written in the file `Significant_genes_list - Copy.txt` (present in each folder). The file must be used from the correct folder depending on which threshold you are using for the analysis. The results of the analysis (`.csv` file and visualization plots such as barplot, dotplot, cnetplot, heatmap, etc.) are also present in each folder for each threshold.

*Note on gene symbols:
For all `enrichGO` analyses, uncharacterized genes which didn't have a formal gene symbol (the ones which were denoted by the `WFL` prefix) were excluded prior to analysis. As a result, the universal background list of 4,641 genes analyzed by DESeq2 was reduced down to a final list of **3,955** annotated genes.*
*The total number of significant genes satisfying each of the three thresholds was similarly reduced:*
* **`no_logfc_threshold`:** 962 annotated genes (filtered down from 1,142 total significant genes).
* **`abs_logfc_greater_than_log2`:** 485 annotated genes (filtered down from 648).
  
* **`abs_logfc_greater_than_log3`:**   234 annotated genes (filtered down from 379).

