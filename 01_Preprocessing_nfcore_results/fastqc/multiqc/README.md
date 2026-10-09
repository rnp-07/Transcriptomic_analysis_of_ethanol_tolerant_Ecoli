In this folder named "multiqc", you will find one html file ("multiqc_report.html") that has compiled all FastQC results into one report. 

The names of the samples or readfiles here are in such a way that the first part of the name tells you the biological condition, so you'll see either "PARENT" for the parental strain or "EVOLVED" for the evolved one. It is followed by replicate number (so REP1, REP2 or REP3 for three replicates).
And then, it is either followed by "_raw" (the original readfiles which were not trimmed) or "_trimmed" (the readfiles obtained after TrimGalore's trimming). 
Then, it is followed by "_1" or "_2" which means _1 is for forward read while _2 is for reverse read. For forward read (_1) and reverse read(_2) of trimmed files only, it is also followed by "_val_1" and "_val_2" respectively which is just TrimGalore's way of naming the trimmed files. 

So, for ex- if you see EVOLVED_REP1_raw_1, it means it is a raw read file of EVOLVED replicate-1 and it is a forward readfile. Similarly, if you see EVOLVED_REP3_trimmed_2_val_2, it means it is a trimmed file for EVOLVED replicate-3 and it is a reverse readfile.

Thus, for raw readfiles for 3 replicates of each of the condition (PARENT vs EVOLVED), there should be 12 raw readfiles and for trimmed readfiles for those 12 raw readfiles, there should be another 12 trimmed readfiles.  