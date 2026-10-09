## Read Alignment

This main directory is divided into two sub-directories, isolating the whole-genome alignment from the synthetic operon alignment.

### Directory Structure

* **`Alignment_of_all_reads_against_ref_genome/`**:
Contains the R script and summary results for aligning the trimmed paired-end reads directly against the *E. coli* W reference genome.

* **`Alignment_of_unmapped_reads_against_operon/`**:
Contains the bash script required to extract unmapped paired-end readfiles from the whole-genome alignment files, along with the R script to align those specific unmapped read pairs against the synthetic operon sequence.

*Note: Each of these folders contains its own dedicated `README.md` with granular instructions, exact file dependencies, and expected outputs for that specific step.*