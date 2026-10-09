#!/bin/bash
#This script needs to be run to get reference genome annotation file in bed format (original is in gff format), to get bam files sorted and indexed which are necessary for some of the RSeQC scripts and to run specific RSeQC scripts


#Change directory to the location where gff file is present by uncommenting the following line and putting the file location in the placeholder
#cd <path_to_gff_file>

#Convert gff file to genePred
gffread genomic.gff -T -o genomic.gtf
gtfToGenePred genomic.gtf genomic.genePred

#Convert the resulting genePred file into a standard BED12 format
genePredToBed genomic.genePred genomic.bed

#Then, change directory to the location where .bam files are present by uncommenting the following line and putting the file location in the placeholder
#cd <path_to_your_bam_files>

#Make an array containing prefixes of the bam files
sample=("Pr_1_alignment" "Pr_2_alignment" "Pr_3_alignment" "Ev_1_alignment" "Ev_2_alignment" "Ev_3_alignment")

#First run the loop to make sorted and indexed bam files

for i in "${sample[@]}"; do
#The following line does sorting using 4 CPU cores. If you want to change, you can change 4 to your desired number depending on your computer's specifications
samtools sort -@ 4 -o "${i}_sorted.bam" "${i}.bam" 
samtools index "${i}_sorted.bam"
done

#Second, run the loop to run 5 RSeQC scripts for all those sorted bam files

for j in "${sample[@]}"; do 
bam_stat.py -i "${j}_sorted.bam" > "${j}_bam_stat.txt" 
infer_experiment.py -i "${j}_sorted.bam" -r genomic.bed -q 0 > "${j}_infer_experiment.txt"  
read_distribution.py -i "${j}_sorted.bam" -r genomic.bed > "${j}_read_distribn.txt"
inner_distance.py -i "${j}_sorted.bam" -r genomic.bed -o "${j%_alignment}_inner_distance"
read_duplication.py -i "${j}_sorted.bam" -q 0 -o "${j}"
done

#Third, run geneBody_coverage.py once for all samples combined
geneBody_coverage.py -i Pr_1_alignment_sorted.bam,Pr_2_alignment_sorted.bam,Pr_3_alignment_sorted.bam,Ev_1_alignment_sorted.bam,Ev_2_alignment_sorted.bam,Ev_3_alignment_sorted.bam -r genomic.bed -o gene_body_coverage
