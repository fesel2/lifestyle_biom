suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(phyloseq)
  library(SIAMCAT)
  library(MintTea)
  library(pheatmap)
  
  
  library(igraph)
  library(tidygraph)
  library(ggraph)
})



wd <- "/home/lukas/Documents/Master/biomML/smoking_microbiome/"
load(paste0(wd, "data/workspace.RData"))

# prepare the taxa data
mask <- sc@associations$assoc.results$p.adj < 0.05
da_species <- rownames(sc@associations$assoc.results[mask,])

OTU_subset <- sc@phyloseq@otu_table[rownames(sc@phyloseq@otu_table) %in% da_species,]
Tax <- as.data.frame(t(sc@phyloseq@otu_table))

# string formating the column names
colnames(Tax) <- lapply(colnames(Tax), function(x)
  paste0("T__", x))

Tax$ID <- rownames(Tax)

# prepare the metabolome data m
metabolome <- read.xlsx(paste0(wd, "data/FG500_metabolomics.xlsx"),
                        sheet = 2, startRow = 2, rowNames = TRUE)

# get the annotation info
M_description <- c(metabolome$annotation)
names(M_description) <- rownames(metabolome)

# get rid of annotation and transpose
# data is already z-scored, --> fill Nas with 0 
m <- metabolome %>%
  select(-annotation) %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
colnames(m) <- sapply(colnames(m), function(x) 
  paste0("M__", x))

m$ID <- rownames(m)

# prepare immunological cytokine data
cytokines <- read.xlsx(paste0(wd, "data/FG500_ELISA_2024-05-03_11_25_06.xlsx"),
                       sheet = 1)
rownames(cytokines) <- cytokines$ID_500FG

# Function to calculate mode
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# replace Nas with mode to tackle detection limit, if Na is mode impute with mean
# divide by standard deviation as normalization

c <- cytokines %>%
  select(-auto_id, -id_computed, -ID_500FG) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), get_mode(.), .))) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  mutate(across(everything(), ~ . / sd(., na.rm = TRUE))) 

colnames(c) <- paste0("C__", colnames(c))
# unfortunately the labels used in cmd are different, however in the original 
# study origin, the label pairs may be found

c$ID <- rownames(c)

# correlate cytokines optional try to reduce dimensions
cor_matrix <- c %>% select(-ID) %>% cor()
heatmap(cor_matrix, main = "Heatmap Correlation")

# # Convert the correlation matrix to a distance matrix, 
# dist_matrix <- as.dist(1 - cor_matrix)
# 
# # Perform hierarchical clustering
# hclust_result <- hclust(dist_matrix, method = "ward.D2")
# 
# clusters <- cutree(hclust_result, h = 0.5)

# Assuming `phyl.rel@sam_data` is correctly structured and extract relevant columns
sam_data <- as.data.frame(sample_data(phyl.rel))[, c("subject_id", "smoking_status")]
sam <- sam_data %>%
  t() %>%
  t() %>%
  as.data.frame()

sam$ID <- rownames(sam)

df_T_sam <- inner_join(sam, Tax, by = "ID")

df_T_M_sam <- inner_join(df_T_sam, m, by = c("subject_id" = "ID"))

df_T_M_C_sam <- inner_join(df_T_M_sam, c, by = c("subject_id" = "ID"))

final_df <- df_T_M_C_sam %>%
  select(-ID)


# Correlation Analyis of Cytokines and DA species
# Select columns starting with "C" or "T"
T_c_final_df <- final_df[, grep("^C|^T", names(final_df))]
T_c_corr <- as.data.frame(cor(T_c_final_df))

DA_species_T <- unlist(lapply(da_species, function(x)
  paste0("T__", x)))

T_c_corr_da_species <- T_c_corr[rownames(T_c_corr) %in% DA_species_T, grepl("^C", names(T_c_corr))]

# create an appropriate heatmap

my_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
# Create breaks to ensure the colorbar is zero-centered
max_val <- max(abs(T_c_corr_da_species))
breaks <- seq(-max_val, max_val, length.out = length(my_palette) + 1)

pheatmap(T_c_corr_da_species, 
         main = "Correlation of immunological features with differentially abundant species",
         cluster_rows = FALSE, 
         cluster_cols = TRUE,
         scale = "none", 
         color = my_palette, 
         breaks = breaks,
         show_rownames = TRUE, 
         show_colnames = TRUE,
         fontsize_col = 5)

# Save the current graphics parameters
par(mar = c(5, 4, 4, 6) + 0.1)  # Adjust the margins for the color bar label

# Run MintTea, the integration software to find biological modules
mt_results <- MintTea(proc_data = final_df, 
                      study_group_column = "smoking_status", 
                      control_group_name = "no smoker",
                      case_group_name = "current smoker",
                      sample_id_column = "subject_id",
                      view_prefixes = c("T", "M", "C"))

save.image(paste0(wd, "backup.RData"))

mt_results$`keep_10//des_0.5//nrep_10//nfol_5//ncom_5//edge_0.8`$module1$features

network <- mt_results$`keep_10//des_0.5//nrep_10//nfol_5//ncom_5//edge_0.8`$module3$module_edges


# Load packages
# run this for every module and adapt graph visualisations

edges <- mt_results$`keep_10//des_0.5//nrep_10//nfol_5//ncom_5//edge_0.8`$module2$module_edges
nodes <- as.data.frame(mt_results$`keep_10//des_0.5//nrep_10//nfol_5//ncom_5//edge_0.8`$module2$features)
names(nodes) <- "feature"

# Apply function to split each string and extract the second part

# function to split all strings in a vector

split_strings <- function(vector, n){
  as.vector(sapply(vector, function(x) {
    split_result <- strsplit(x, "__")[[1]]
    return(split_result[n])
  }))
}
nodes$name <- split_strings(nodes$feature, 2)
nodes$omic <- split_strings(nodes$feature, 1)


edges$feature.x <- split_strings(edges$feature.x, 2)
edges$feature.y <- split_strings(edges$feature.y, 2)
filtered_edges <- edges %>%
  filter(edge_weight > 0.8)

nodes <- nodes %>%
  select(c("name", "omic")) %>%
  filter(name %in% filtered_edges$feature.x | name %in% filtered_edges$feature.y)

# Create igraph object from edge list
graph <- graph_from_data_frame(d = filtered_edges, vertices = nodes, directed = FALSE)

ggraph(graph, layout = "graphopt") + 
  geom_edge_link(aes(edge_alpha = 0.1, edge_width = 10), edge_colour = "grey") +
  geom_node_point(aes(color = omic, size = 10)) + 
  geom_node_text(aes(label = name), repel = TRUE, point.padding = unit(0.2, "lines")) +
  theme_void()

