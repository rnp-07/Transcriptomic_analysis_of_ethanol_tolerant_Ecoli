# TPM Calculation

This directory contains an R script `Finding_gene_length_from_gff.R` that parses the reference annotation (`genomic.gff`) to calculate the exact length of each gene and pseudogene and an excel file `tpm_and_percentiles.xlsx` that consists of 4 sheets:

  * **`tpm_of_all_genes`**: Contains the calculated TPM for all protein-coding genes and pseudogenes (excluding tRNA/rRNA biotypes), with the chromosomally integrated *pdc* and *adhB* genes appended.
  * **`average_tpm_and_percentiles`**: Contains genes ordered by their average TPM percentiles in descending order, calculated separately for both the Parent and Evolved strains.
  * **`pdc_and_adhB`**: Highlights native genes whose expression levels fall within a close margin (±1%) of the average TPM counts of the integrated *pdc* and *adhB* genes in both the Parent and Evolved strains.
  * **`Mutated_genes`**: Contains gene expression data specifically related to the mutated loci (*yeeJ* and *waaV* loci).