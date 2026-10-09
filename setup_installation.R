# Check and install BiocManager if it is missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# These are the packages required for this analysis in R
req_pkgs <- c("Rsubread", "biomartr", "DESeq2", "ggplot2", "dplyr", "pheatmap", "ggrepel", "tidyr", "reshape2", "scales", "EnhancedVolcano", "clusterProfiler", "org.EcK12.eg.db", "enrichplot", "pathview")


#Save the list of installed packages in a vector
inst<- installed.packages()[, "Package"]

#Find missing packages
missing_pkgs <- req_pkgs[!(req_pkgs %in% inst)]

# Install only the missing packages using BiocManager
if(length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  BiocManager::install(missing_pkgs, update = FALSE)
} else {
  message("All required packages for the pipeline are already installed!")
}