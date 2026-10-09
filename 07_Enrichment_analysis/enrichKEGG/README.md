# KEGG Pathway Enrichment Analysis

This directory contains the script, input files and the outputs of over-representation analysis (ORA) using the KEGG database using `enrichKEGG` function from the `clusterProfiler` package. Similar to GO analysis, this analysis was also performed independently across three different fold-change thresholds: i) no logfc threshold ii) $|\log_2\text{fc}|>log_22$ and iii) $|\log_2\text{fc}|>log_23$


### Input Files
* **`Background_list.txt`**: The universal background list of all 4,641 genes analyzed by DESeq2.
* **`gene_names_with_logfc.txt`**: A txt file containing the log2 fold-change values for all analyzed genes. This is required by the `pathview` package to put color gradients based on fold-change onto the KEGG pathway maps.
### Note on ID mapping
The `enrichKEGG` function specifically requires Entrez IDs rather than standard Gene Symbols. Therefore, all gene lists were mapped to *E. coli* K-12 Entrez IDs using the `org.EcK12.eg.db` database prior to analysis. Because some uncharacterized or strain-specific genes do not have a direct K-12 Entrez counterpart, the total number of analyzed genes drops during mapping:
* **Background List:** 4,641 total genes $\rightarrow$ 3,686 successfully mapped Entrez IDs.
* **No logFC Threshold:** 1,142 significant genes $\rightarrow$ 906 successfully mapped Entrez IDs.
* **$|\log_2\text{fc}|>log_22$:** 648 significant genes $\rightarrow$ 462 successfully mapped Entrez IDs.
* **$|\log_2\text{fc}|>log_23$:** 379 significant genes $\rightarrow$ 221 successfully mapped Entrez IDs.


### Script
* **`Script_for_enrichKEGG_analysis.R`**: The main R script used to map identifiers, execute the KEGG enrichment analysis, and generate standard plots (Barplot, Dotplot, Emapplot) and Pathview maps.

    *NOTE: This script outputs the results as a .csv file where the identifier used are EntrezIDs rather than Gene Symbols. So, in oder to convert those EntrezIDs back to Gene Symbols, you can run `Script_for_converting_entrezids_into_symbols.R`.*
* **`Script_for_converting_entrezids_into_symbols.R`**: As mentioned above, this script was used to convert EntrezIDs back into Gene Symbols for the final result tables.

### Threshold Subdirectories
Similar to GO analysis, there are 3 subdirectories here as well which are related to 3 thresholds  taken for KEGG pathway enrichment analysis: "`no_logfc_threshold/` (all significant genes i.e. $p_\text{adj} < 0.05$)", "`abs_logfc_greater_than_log2/`" (i.e. significant genes satisfying $|\log_2\text{fc}|>log2$ threshold) and "`abs_logfc_greater_than_log3/`"(i.e. significant genes satisfying $|\log_2\text{fc}|>log3$ threshold). For each threshold, the list of complying genes is written in the file `Significant_genes_list - Copy.txt` (present in each folder). The file must be used from the correct folder depending on which threshold you are using for the analysis. The results of the analysis are also present in each folder. 

*(**Important:** Don't get confused with `.xlsx` file and `.csv` file present in these folders as both contain the same enrichment analysis output stats for each analysis except the fact that `.xlsx` file also includes Gene Symbols for easier review.*
