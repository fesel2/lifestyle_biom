# check distributions of sample and bacteria
quality_checks <- function (feature.table, folder) {
  message("QC: Generating Histograms")
  # Quality check Reads per sample
  png(paste0(folder, "reads_per_sample.png"))
  reads_sample <- hist(colSums(feature.table), xlab="Reads for Sample")
  dev.off()
  # Quality check Reads per Bacteria
  png(paste0(folder, "reads_per_feature.png"))
  reads_bac <- hist(log1p(rowSums(feature.table)), xlab="log10+1 of Reads for Feature")
  dev.off()
}

create_compositional_table <- function(phylo, group="phylum") {
  # temporary merge the grouping column to the otu_table
  motu.rel.fil.tax <- merge(phylo@otu_table@.Data,
                            phylo@tax_table@.Data[, group, drop = FALSE],
                            by = "row.names", all.x = TRUE)
  # remove index information to support tibble usage
  rownames(motu.rel.fil.tax) <- motu.rel.fil.tax$row.names
  motu.rel.fil.tax$Row.names <- NULL  # Removing the redundant column
  # group, sum up and transpose
  summed_phylum <- 
    motu.rel.fil.tax %>%
    group_by(!!sym(group)) %>%
    summarise_all(sum) %>%
    filter(!is.na(!!sym(group))) %>%
    column_to_rownames(var = group) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    as_tibble()
  return(summed_phylum)
}

select_comp_features <- function(comp_table, lab) {
  tmp_label_df <- phyl.rel@sam_data[, c(lab)]
  tmp_label_df$ID <- rownames(tmp_label_df)
  summed_phylum_long <-
    comp_table %>%
    gather(key = "columns", value = "values", -sample) %>%
    left_join(tmp_label_df, by = c("sample" = "ID"))
  return(summed_phylum_long)
}

create_figure <- function(comp_features, f) {
  n <- length(unique(comp_features$sample))
  summed_phylum_fig <- comp_features %>%
    ggplot(aes(x = sample, y = values, fill = columns)) +
    geom_bar(stat = "identity", width = 1) +
    labs(x = paste0("Samples, n=", n), y = "Phylum (%)", fill = "Phyla") +
    scale_x_discrete(breaks = NULL) +
    theme(axis.text.x = element_blank()) +
    theme_minimal() +
    guides(
      fill = guide_legend(
        ncol = 1  # Set the number of columns in the legend
      )
    )
  ggsave(paste0(f, "bacterial_composition.png"), 
         plot = summed_phylum_fig,
         width = 10,
         height = 4)
  return(summed_phylum_fig)
}

# create the heatmap
create_heatmap <- function (heatmap_data, feat) {
  # Extract unique categories
  cats <- unique(heatmap_data[[feat]])
  # Create a named vector of colors
  color_mapping <- setNames(c("red", "#F0F0F0", "blue", "green"), cats)
  
  hm <- heatmap_data %>%
    ggplot(aes(x = rownames(heatmap_data), y = 1, fill = factor(!!sym(feat)))) +
    geom_tile() + 
    scale_fill_manual(values = color_mapping) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(fill = feat)+
    scale_x_discrete(breaks = NULL) +
    theme(axis.text.x = element_blank())
  return(hm)
}

# calculate alpha diversity
calc_alpha_div <- function (feature.table, phenotype) {
  alpha.div <- feature.table %>%
    t() %>%
    as_tibble(rownames = "ID") %>%
    pivot_longer(-ID) %>%
    group_by(ID) %>%
    summarize(richness = specnumber(value),
              shannon = diversity(value, index="shannon"),
              simpson = diversity(value, index="simpson"),
              invsimpson = 1/simpson,
              n = sum(value))
  phenotype$ID <- rownames(phenotype)
  alpha.div.meta <- left_join(alpha.div, phenotype, by="ID", keep=FALSE)
  return(alpha.div.meta)
}

# plot alpha div results
create_plot_alpha_div <- function (alpha.div.table, folder) {
  plot_alpha_div <- alpha.div.table %>%
    pivot_longer(cols=c(richness, shannon, invsimpson, simpson), 
                 names_to="metric") %>%
    ggplot(aes(metric, value, fill = feature)) + 
    geom_boxplot() + 
    stat_compare_means() + # wilcox.test p-value
    facet_wrap(. ~ metric, scale = "free") +
    ggtitle(paste("Alpha Diversity")) +
    theme_classic() 
  ggsave(paste0(folder, "alpha_div.png"), 
         plot = plot_alpha_div, width = 10, height = 7)
}


run_anosim <- function (feature.table, phenotype) {
  phenotype <- na.omit(phenotype)
  features <- as.data.frame(t(feature.table))
  features <- features[rownames(features) %in% rownames(phenotype),]
  dist <- vegan::vegdist(features)
  ano <- vegan::anosim(dist, phenotype$feature)
  message("BETA-DIV: Analysis of Similarities Statistic: ", ano$statistic)
  message("BETA-DIV: Analysis of Similarities Significance: ", ano$signif)
  return(ano)
}

# run statistics to calculate differential abundances
l2fc_wilc <- function(feature.table, phenotype, folder, p_cutoff = 0.05) {
  message("DA-WILC: Testing for associations with wilcoxon")
  phenotype$ID <- rownames(phenotype)
  #make sure to drop the NAs
  phenotype <- na.omit(phenotype)
  uniq_cats <- unique(phenotype$feature)
  
  # initialize a empty dataframe to hold the results
  p.cal <- tribble(~ID, ~pval, ~adj.pval, ~log10.adj, ~aucs.mat, ~fc, ~sig)
  for (rowname in row.names(feature.table)) {
    
    # get the values indeces corresponding to label
    x_IDs <- phenotype %>%
      filter(feature == uniq_cats[1]) %>%
      pull(ID)
    y_IDs <- phenotype %>%
      filter(feature == uniq_cats[2]) %>%
      pull(ID)
    # make sure only existing columns are included
    x_IDs <- x_IDs[x_IDs %in% colnames(feature.table)]
    y_IDs <- y_IDs[y_IDs %in% colnames(feature.table)]
    x <- as.numeric(feature.table[rowname, x_IDs])
    y <- as.numeric(feature.table[rowname, y_IDs])
    
    q.p <- quantile(log10(x+1e-05), probs=seq(.1, .9, .05))
    q.n <- quantile(log10(y+1e-05), probs=seq(.1, .9, .05))
    
    # create matrix
    p.cal=add_row(p.cal, 
                  ID = rowname,
                  pval = wilcox.test(x, y, exact=FALSE)$p.value,
                  aucs.mat = roc(controls=y, cases=x, direction='<', ci=TRUE, auc=TRUE)$ci[2], 
                  fc = sum(q.p - q.n)/length(q.p))
  }
  # p.adjust
  message("DA-WILC: Adjusting p-values")
  p.cal <- p.cal %>% 
    mutate(adj.pval = p.adjust(pval, method = "BH")) %>%
    mutate(sig = ifelse(adj.pval < p_cutoff, "p.adj < 0.05", "not sig")) %>%
    mutate(log10.adj = -log10(adj.pval))
  
  # save the dataframe
  name <- paste0(folder, "DA_table_", uniq_cats[1], "_vs_", uniq_cats[2], ".tsv")
  write.table(p.cal, file = name, sep="\t", row.names = TRUE, col.names = TRUE)
  return(p.cal)
}

vulcano_plot <- function(DA_table, folder, p_cutoff=0.05) {
  # Visualize with volcano plots
  DA_table <- na.omit(DA_table)
  p.wilcox <- ggplot(DA_table, aes(x = fc, y = log10.adj)) +
    geom_point(aes(color = sig), alpha = 0.7) +
    scale_color_manual(values = c("black", "red")) +
    geom_hline(yintercept = -log10(p_cutoff)) +  # Add a line for significance
    geom_text_repel(data = subset(DA_table, adj.pval < p_cutoff), aes(label = ID), size = 3) +  # Add ID for significant species
    ggtitle("Differential Abundance") +
    xlab("Generalized fold change") + 
    ylab("Log10 adjusted p-value") +
    theme_classic()
  ggsave(paste0(folder,"DA_vulcano.png"), 
         plot = p.wilcox,
         width = 10,
         height = 7)
}