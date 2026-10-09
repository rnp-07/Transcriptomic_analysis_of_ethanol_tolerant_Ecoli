#!/bin/bash
#The following line was the script that was used for getting unmapped read pairs (reads and their mate both being unmapped) and saving them separately as R1 and R2 read files for Evolved-replicate-1. Similar command was used by changing the input and output filenames for rest of the samples.

samtools fastq -f 12 -1 unmapped_R1_Ev1.fastq -2 unmapped_R2_Ev1.fastq Ev_1_alignment.bam



