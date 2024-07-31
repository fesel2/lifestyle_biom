#!/usr/bin/env Rscript

# IMPORT MODULES AND DATA ......................................................

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(phyloseq)
  library(ggpubr)
  library(ggplot2)
  library(pROC)
  library(ggrepel)
  library(SIAMCAT)
  library(patchwork)
})

wd <- "/home/lukas/Documents/Master/biomML/smoking_microbiome/"

# import the prepared helper functions
source(paste0(wd, "03_helperfunctions.R"))
# prepare an outputfolder
folder <- paste0(wd, "results/")

# read the data
phyl <- readRDS(paste0(wd, "data/SchirmerM_2016.rds"))
message("IMPORT: Data successfully read in")

# FILTERING ....................................................................

# plot QC distributions
quality_checks(phyl@otu_table, folder)

# save the unfiltered phyloseq object before further processing
phyl.raw <- phyl

# only keep informative taxa which are found in at least 20% of the samples
phyl.shape1 <- dim(phyl@otu_table)
phyl <- filter_taxa(phyl, function(x) sum(x > 0) > (0.2*length(x)), TRUE)
phyl.shape2 <- dim(phyl@otu_table)

message("FILTERING: ",phyl.shape1[1] - phyl.shape2[1], 
        " Taxa are not found in at least 20% of the samples (n=",
        round(0.2*dim(phyl@otu_table)[2]),
        ") and discarded for futher processing.")

# Transform counts into relative abundances and save as phyl.rel
phyl.rel <- transform_sample_counts(phyl, function(x) x/sum(x))

# Filtering for taxa which cover at least 0.00001 of the total relative
# abundance as mean
phyl.rel <- filter_taxa(phyl.rel, function(x) mean(x) > 1e-5, TRUE)

phyl.shape3 <- dim(phyl.rel@otu_table)

message("FILTERING: ", phyl.shape2[1] - phyl.shape3[1],
        " taxa cover less than 0.001% on average across the samples 
        and therefore filtered out.")

# Filter the samples: keeo only current and no smoker

# generate a new column including the ex smoker column
n <- phyl.rel@sam_data[, c("smoker", "ever_smoker")]
n$smoking_status <- ifelse(n$smoker == "yes", "current smoker", "no smoker")
n$smoking_status[n$smoker == "no" & n$ever_smoker == "yes"] <- "former smoker"

phyl.rel@sam_data$smoking_status <- n$smoking_status
phyl@sam_data$smoking_status <- n$smoking_status

# drop the former smokers
smoke_only <- phyl.rel@sam_data$smoking_status == "current smoker" | phyl.rel@sam_data$smoking_status == "no smoker"
phyl.rel.smoke_only <- prune_samples(smoke_only, phyl.rel)
phyl.rel <- phyl.rel.smoke_only

phyl.smoke_only <- prune_samples(smoke_only, phyl)
phyl<- phyl.smoke_only


phyl.shape4 <- dim(phyl.rel.smoke_only@otu_table)

message("FILTERING: ", phyl.shape3[2] - phyl.shape4[2],
        " samples are either NAs or former smoker and discarded")

# print Final dimensions
message("FILTERING: A total of ", phyl.shape4[1], " taxa are remaining and ",
        phyl.shape4[2], " samples are remaining")
feature <- "smoking_status"

class_balance <- table(phyl.rel@sam_data[, feature])
print(class_balance)

# extract the phenotype of interest for simplified processing
phen <- as.data.frame(phyl.rel@sam_data[, feature])
names(phen) <- "feature"

# COMPOSITION ANALYIS ..........................................................

# use prepared functions to extract phylum data and built a stacked barplot
cT <- create_compositional_table(phyl.rel, group = "phylum")
cT_lf <- select_comp_features(cT, feature)

# sort samples by Firmicutes amount
cT$firm_bac_ratio <- cT$Firmicutes / cT$Bacteroidota
sample_order <- cT[order(cT$firm_bac_ratio, decreasing = FALSE),]$sample
cT_lf$sample <- factor(cT_lf$sample, levels = sample_order)

# create the figure 
fig <- create_figure(cT_lf, folder)

# Create the heatmap data
hm_data <- as.data.frame(phyl.rel@sam_data[sample_order, feature])
hm <- create_heatmap(hm_data, feature)

# combine and save the plots
combined_plot <- hm +  fig + plot_layout(heights = c(1, 10))
ggsave(paste0(folder, "bacterial_composition_hm.png"), 
       plot = combined_plot, width = 10, height = 4)

message("COMPOSITIONAL: stacked barplot with feature heatmap successfully 
        created")

# Firm/Bac Ratios
tmp <- phyl.rel@sam_data[, feature]
tmp$sample <- rownames(tmp)
cT <- left_join(cT, tmp, by = "sample")

# statistical test
cats <- unique(cT[[feature]])
cat1 <- cT[cT[[feature]] == cats[1],]$firm_bac_ratio
cat2 <- cT[cT[[feature]] == cats[2],]$firm_bac_ratio

test_result <- wilcox.test(cat1, cat2)

bplot_ratio <- ggplot(cT, aes(x = !!sym(feature), y = firm_bac_ratio, fill= !!sym(feature))) +
  geom_boxplot(outlier.shape = NA, coef = 1.98) +
  geom_jitter(width = 0.1, size = 0.4) +
  ylim(0, 10) +
  geom_text(x = 1.5, y = 10, 
            label = paste("wilc-test p-value:", format.pval(test_result$p.value, digits = 2))) +
  labs(y = "Firmicutes Bacteroidota ratio", x = NULL)

ggsave(paste0(folder, "Firm_bac_ratio.png"), 
       plot = bplot_ratio, width = 10, height = 7)

message("COMPOSITIONAL: Firmicutes Bacterio Ratio compared with 
        wilcoxon and boxplot plotted")

# ALPHA DIVERSITY ANALYSIS .....................................................

# use the prepared function to calculate diversity metrics with vegan 
# on absolute count values
alph_div_table <- calc_alpha_div(phyl@otu_table@.Data, phen)

# plot the results
create_plot_alpha_div(alph_div_table, folder)
message("ALPHA-DIV: Results succesfully plotted to ", folder)


# BETA DIVERSITY ANALYSIS ......................................................

# initialize distance metrics and results table
anosim_table <- tribble(
  ~distance_metric, ~anosim_R, ~anosim_q, ~permanova_R2, ~permanova_F, ~permanova_pval
)

# Define the list of distance metrics
distances <- c("euclidean", "bray", "jaccard", "uunifrac", "wunifrac")
# need to convert phen to dataframe for permanova
phen_tmp <- as.data.frame(phen@.Data)
names(phen_tmp) <- "feature"

for (d in distances) {
  # Calculate distance matrix
  dDist <- phyloseq::distance(phyl.rel, method = d)
  # Visualize with NMDS
  dMDS <- phyloseq::ordinate(phyl.rel, "NMDS", distance = dDist)
  p <- phyloseq::plot_ordination(phyl.rel, dMDS, color = feature) +
    stat_ellipse(aes(color = feature), level = 0.9) +
    ggtitle(paste0("NMDS of ", d, " distance matrix"))
  ggsave(paste0(folder, "NMDS_", d, "distance.png"), 
         plot = p, width = 10, height = 10)
  
  # Run ANOSIM and save results in the table
  features <- as.data.frame(t(otu_table(phyl.rel)@.Data))
  anosim_result <- vegan::anosim(dDist, phen$feature, permutations = 250)
  permanova_result <- adonis2(dDist ~ feature, data = phen_tmp, permutations = 250)
  
  # save results in table
  anosim_table <- add_row(anosim_table,
                          distance_metric = d,
                          anosim_R = anosim_result$statistic,
                          anosim_q = anosim_result$signif,
                          permanova_R2 = permanova_result$R2[1],
                          permanova_F = permanova_result$F[1],
                          permanova_pval = permanova_result$`Pr(>F)`[1]
                          )
}

# Display the ANOSIM results table
print(anosim_table)

# SIAMCAT MODELING AND ASSOCIATIONTESTING  .....................................

# for visualisation purposes the rownames needs to be renamed to the species only
old_row_names <- rownames(phyl.rel@otu_table)
rownames(phyl.rel@otu_table) <- lapply(old_row_names, function(x)
  strsplit(x, "s__")[[1]][2])

old_row_names <- rownames(phyl.rel@tax_table)
taxa_names(phyl.rel) <- lapply(old_row_names, function(x)
  strsplit(x, "s__")[[1]][2])

cats <- unique(phyl.rel@sam_data[[feature]])
# make sure the case is yes --> so yes is minority class
sc <- SIAMCAT::siamcat(phyloseq = phyl.rel, label = feature, case = "current smoker")

# is already filtered from previous analyis, however needed to proceed
sc <- filter.features(sc, filter.method = 'abundance', cutoff=0.001)

sc <- normalize.features(sc, norm.method = "log.clr", 
                         norm.param = list(log.n0 = 1e-05, sd.min.q = 1))

sc <- check.associations(sc, feature.type = 'normalized', log.n0 = 10e-5)



sc <- create.data.split(sc, num.folds = 5, num.resample = 1)

# create visualisation for the statistical associations
result <- try({
  # Code that might produce an error
  association.plot(sc, fn.plot=paste0(folder,"asso.pdf"))
  volcano.plot(sc, fn.plot=paste0(folder,"vulc.pdf"))
}, silent = TRUE)

if (inherits(result, "try-error")) {
  message("An error occurred, very likely no differencially abundant features
          were found, no output figure produced")
  result <- NA
}

assos <- sc@associations$assoc.results
sig_assos <- assos[assos$p.adj <= 0.05,]
message("ASSOCIATIONS: ",dim(sig_assos)[1], " significant features (differentially 
        abundant species) were found.")

sc <- train.model(sc, method = "randomForest")

sc <- make.predictions(sc)

sc <- evaluate.predictions(sc)

# create visualisation to interpret models
result <- try({
  # Code that might produce an error
  model.evaluation.plot(sc, fn.plot=paste0(folder,"roc.pdf"))
  model.interpretation.plot(sc, consens.thres = 0.01, 
                            fn.plot=paste0(folder,"int.pdf"))
}, silent = TRUE)

if (inherits(result, "try-error")) {
  message("An error occurred, check whats wrong, maybe not enough features?")
  result <- NA
}

# print evaluation metrics
message("Classifier Performance, AUROC: ", sc@eval_data$auroc)
message("Classifier performance, auprc: ", sc@eval_data$auprc)

# SAVING AND CLEANUP............................................................

message("FINISHING: Saving the data and exiting...")

# Specify the objects you want to keep
keep_objects <- c("alpha_div_table", "anosim_table", "sc", "phyl", "phyl.raw",
                  "phyl.rel", "folder", "wd")

# Remove all objects except those you want to keep
rm(list = setdiff(ls(), keep_objects))

save.image(paste0(wd, "data/workspace.RData"))

sessionInfo()

