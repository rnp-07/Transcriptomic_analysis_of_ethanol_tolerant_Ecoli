## Alignment of Trimmed Reads Against the Reference Genome

This directory contains the R script and summary data used to align the pre-processed, trimmed paired-end reads against the *E. coli* W reference genome.

### Files

* **`Alignment_of_all_reads_against_ref_genome.R`**: R script responsible for building the reference genome index and aligning all the trimmed paired-end read files against the reference sequence (`GCF_000258145.1_ASM25814v1_genomic.fna`).
* **`BAM_Summary.csv`**: A compiled spreadsheet containing the alignment metrics and statistics extracted from the individual `.BAM.summary` files outputted for each biological replicate after the alignment.

### Outputs

**1. Reference Index**:
Before mapping the reads, the R script must generate an index of the reference genome. This process creates 6 internal files in the working directory (using the prefix specified in the R script):
`.00.b.array` | `.00.b.tab` | `.lowinf` | `.reads` | `.files` | `.log`

**2. Alignment Results**
Once the alignment process is completed, 3 resultant files are generated for each individual sample:
`.BAM` | `.BAM.indel.vcf` | `.BAM.summary`

### Note on Data Availability & Large Files
Due to GitHub's file size limitations, the raw sequencing files, generated reference index files, and the resultant alignment files (e.g., `.BAM`, `.vcf`) are not hosted in this repository. These files are fully reproducible and will be automatically generated on your local machine when you execute the provided scripts.