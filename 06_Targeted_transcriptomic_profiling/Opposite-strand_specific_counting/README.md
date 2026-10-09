# Opposite Strand-Specific Counting:

### Background
The initial read quantification after the whole-genome alignment was performed using reverse-stranded settings (`strandSpecific=2` for all samples except the unstranded Parent-replicate-1) based on RSeQC results (see `infer_experiment.py` results). However, prior genomic analysis had identified a mutation in the *waaV* gene of the evolved strain, raising a specific question: **Did this mutation disrupt transcriptional termination?**

If you see the orientation, *waaV* sits on the negative strand, while its upstream neighbor, *rfaC*, is on the positive strand. If *waaV* fails to terminate properly, its elongated transcripts will overlap *rfaC*. Because the library preparation was inferred to be reverse-stranded from RSeQC script, read1s map in the opposite direction of their source genes. Therefore, any read-through transcripts originating from *waaV* and overlapping with *rfaC* would have read1s map onto *rfaC* in the *forward* direction. 

To detect this specific read-through event, reads were re-quantified using a forward-stranded parameter (`strandSpecific=c(1,1,1,1,1)`) for all 5 replicates (`Pr_2`, `Pr_3`, `Ev_1`, `Ev_2` and `Ev_3`). 

### Files
This directory contains data evaluating the raw read counts for each gene under these two opposing strandedness assumptions.

* **`counts_considering_opp_strandedness.xlsx`**: The primary data workbook containing side-by-side raw read counts for all genes across the 5 analyzed samples (`Pr_2`, `Pr_3`, `Ev_1`, `Ev_2`, `Ev_3`) under both strandedness assumptions (forward and reverse).

  Alongside the raw counts, the sheet also includes "Ratio (Forward/Reverse)" and "Log2(Forward/Reverse)" to quantitatively assess the difference in read counts between the two strandedness assumptions.
  
  The highly relevant genes belonging to the *rfa* operon (including *waaV* and *rfaC*) are manually highlighted in yellow within the spreadsheet for quick assessment.