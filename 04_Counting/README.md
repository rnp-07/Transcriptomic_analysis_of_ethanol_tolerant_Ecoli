# Counting (Read Quantification)

This directory contains the R scripts and output count matrices for quantifying the aligned RNA-Seq reads using `featureCounts` (from the `Rsubread` package). The quantification process is divided into four distinct steps.

### Directory Structure & Workflow

* **`01_Gene_plus_trrna_counts_with_Pr1/`**
This folder contains an R script (`Gene_plus_trrna_counts_with_Pr1(strand_specific).R`) and the resulting count matrix (`Gene_plus_trrna_counts_with_Pr1_ss.txt` file) for all features annotated as "gene" in the reference GFF file (`genomic.gff`) and grouped by their unique "ID" attribute. This captures both protein-coding genes and non-coding RNA biotypes (e.g., rRNA, tRNA).

* **`02_Pseudogene_counts_with_Pr1/`**
This folder contains an R script (`Pseudogene_counts_with_Pr1_ss.R`) and the resulting count matrix (`Pseudogene_counts_with_Pr1_ss.txt` file) for features annotated as "pseudogene" in the reference GFF file (`genomic.gff`) which are grouped by their unique "ID" attribute.

* **`03_adhb_pdc_counts_from_unmapped/`**
Contains a custom GFF (`operon.gff`), an R script (`adhb_pdc_counts(unmapped)_ss.R`), and the resulting count matrix (`adhb_pdc_counts(unmapped)_ss.txt` file) for the chromosomally integrated ethanologenic genes (*pdc* and *adhB*). Because these synthetic genes were absent from the wild-type reference genome, read pairs hadn't mapped onto them during the initial whole genome-alignment. So, the unmapped read pairs which were extracted from the original BAM files using SAMtools and aligned to the synthetic operon sequence (detailed in the directory `02_Alignment/Alignment_of_unmapped_reads_against_operon/`) are quantified in this folder.

* **`04_Counts_of_genes_and_pseudogenes(wo_trrna)_with_Pr1/`**
Contains an R script `trna_rrna_gene_list.R` that extracts the IDs of non-coding RNA biotypes from the GFF file and another R script `genes_plus_pseudogenes_counts_wo_trrna_with_Pr1_ss.R` that merges the two previously generated count matrices (from folders 1 and 2), and filter out the non-coding RNA counts (rRNA, tRNA).

### Note on Data Availability & Large Files
Due to GitHub's file size limitations, the intermediate `.BAM` files required to generate these counts are not hosted in this repository. All final count matrices (`.txt`) are provided within their respective sub-directories to allow immediate downstream differential expression analysis.