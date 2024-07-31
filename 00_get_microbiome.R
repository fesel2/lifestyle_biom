#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(mia)
  library(curatedMetagenomicData)
})

folder.data <- "/home/lukas/Documents/Master/biomML/smoking_metagenomics/data/"

# 1. prepare the data sets of interest from curated metagenomic data
filtered <- sampleMetadata %>%
  filter(study_condition == "control")

get_build_phyloseq <- function(study) {
  tryCatch({
    # Debug print
    print(paste("Processing study:", study))
    
    # Subset for study
    samples <- filtered[filtered$study_name == study,]
    print(paste("Number of samples for", study, ":", nrow(samples)))
    
    # Get the associated samples, note: the absolute counts are requested
    Study.tse <- curatedMetagenomicData::returnSamples(
      samples,
      "relative_abundance",
      counts = TRUE,
      rownames = "long")
    print(paste("Sample extraction successful for", study))
    
    # Convert to phyloseq object
    Study.phyl <- mia::makePhyloseqFromTreeSummarizedExperiment(
      Study.tse, 
      assay.type = "relative_abundance")
    print(paste("Conversion to phyloseq object successful for", study))
    
    # Save it to disk for further analysis
    f <- paste0(folder.data, study, ".rds")
    saveRDS(Study.phyl, file = f)
    print(paste("Saved phyloseq object for", study, "to", f))
    
    # Clean up
    rm(list = c("Study.tse", "Study.phyl"))
  }, error = function(e) {
    message(paste("Error processing study", study, ":", e$message))
  })
}

# Studies of interest
studs = c("SchirmerM_2016", "YachidaS_2019",
          "KeohaneDM_2020")

for (stud in studs) {
  get_build_phyloseq(stud)
}
