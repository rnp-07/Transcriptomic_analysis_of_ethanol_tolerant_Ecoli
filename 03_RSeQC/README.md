# Post-Alignment Quality Control (RSeQC)
This directory contains the script, reference annotations, and output reports used to evaluate the integrity of the aligned RNA-Seq libraries using RSeQC (v5.0.4).

### Prerequisites & Environment Setup

The specific tools that were used in this workflow were:
* **RSeQC suite**: v5.0.4
* **gffread**: v0.12.7
* **ucsc-gtftogenepred**: v482
* **ucsc-genepredtobed**: v482
* **samtools**: v1.22.1

To install these exact dependencies automatically, set up an isolated environment:

```bash
# Create the environment and install exact versions of RSeQC, SAMtools, gffread, and UCSC conversion tools
conda create -n rseqc_env -c bioconda -c conda-forge rseqc=5.0.4 gffread=0.12.7 ucsc-gtftogenepred=482 ucsc-genepredtobed=482 samtools=1.22.1 -y

# Activate the environment
conda activate rseqc_env
```
### Directory Structure & Files

**Main Script & Annotations:**
* **`Script_for_running_RSeQC_scripts.sh`**: The main script that converts the annotation file (`.gff` into `.genePred` and then to `.bed`), sorts and indexes the BAM files, and runs the RSeQC suite.
* **`genomic.gff`**: The original reference genome annotation file.
* **`gtf_genepred_bed/`**: A folder containing the reference genome annotation file in `.gtf`, `.genePred` and `.bed` formats.

**Output Directories:**
The results of the 6 RSeQC scripts (`bam_stat.py`, `geneBody_coverage.py`, `infer_experiment.py`, `inner_distance.py`, `read_duplication.py`, and `read_distribution.py`) are present in folders named after them:
* `bam_stat_output/`
* `gene_body_coverage/`
* `infer_experiment_output/`
* `inner_distance_output/`
* `read_distribution_output/`
* `read_duplication_output/`

### Note on Data Availability & Large Files
Due to GitHub's file size limitations, neither the raw nor the sorted BAM files (`.sorted.bam` and `.sorted.bam.bai`) are hosted in this repository.