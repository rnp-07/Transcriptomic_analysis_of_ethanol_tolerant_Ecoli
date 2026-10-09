## Alignment of Unmapped Reads Against the Synthetic Operon

This directory contains the scripts and reference files required to extract unmapped reads from the reference genome alignment and align them against the synthetic engineered operon.

### Files

* **`BASH_script_for_getting_unmapped_reads.sh`**: A bash script that extracts unmapped paired-end read files from the BAM files generated in the previous whole-genome alignment step.
* **`W105_operon.str`**: The original sequence file for the synthetic operon.
* **`operon.fasta`**: The operon sequence manually converted into FASTA format to facilitate index building. (The word "operon" was manually added after the `>` character on the first line of .fasta file).
* **`Alignment_of_unmapped_vs_operon.R`**: R script that builds the index from the `.fasta` file and aligns the unmapped read pairs against the synthetic operon sequence.

### Outputs

When you run the R script to generate the FASTA index, 6 internal files will be generated in the working directory with the prefix specified in the script:
`.00.b.array` | `.00.b.tab` | `.lowinf` | `.reads` | `.files` | `.log`

Once the alignment is complete, 3 resultant files are generated for each sample:
`.BAM` | `.BAM.indel.vcf` | `.BAM.summary`

**Next Step:**
The aligned reads are then quantified. The script and results for this counting step are located in the `04_Counting/03_adhb_pdc_counts_from_unmapped/` directory.

### Note on Data Availability & Large Files
Due to GitHub's file size limitations, the raw sequencing files, generated reference index files, and the resultant alignment files (e.g., `.BAM`, `.vcf`) are not hosted in this repository. These files are fully reproducible and will be automatically generated on your local machine when you execute the provided R and Bash scripts.