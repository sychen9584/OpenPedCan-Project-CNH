# Author: Sam Chen
# Date: 2025-08-08
# Compare immune cell proportions between medulloblastoma subtypes

####### Analysis summary #######
# - Used quanTIseq results as its scores can be compared across cancer subtypes 
#   unlike XCell per the documentation.

# - Most of the deconvoluted cells are labeled as uncharacterized (~90%). I dropped them and 
#   renormalized the fractions of remaining cell types.

# - PCA, correlation, and PERMNOVA test suggest that immune cells significantly distribute
#   differently across mb subtypes

# - Pairwise wilcoxon tests and visualizations (violin+boxplot & stacked barplot) indicate
#   higher distribution of:
#   - myeloid dendritic cells in WNT
#   - B cells in SHH
#   - monocytes in group3
#   - NK cells in group4

####### Output files ########
#   medulloblastoma_pca.png 
#   medulloblastoma_cor_heatmap.png
#   medulloblastoma_stacked_barplot.png
#   medulloblastoma_violin_boxplot.png 
#   medulloblastoma_parwise_significance.png - effect sizes for each comparison with significance annotated
#   medulloblastoma_violin_boxplot_significance.pdf - violin boxplots with significance annotated
#   medulloblastoma_pairwise_stats.csv

### MAIN ANALYSIS ###
####### Load libraries #######

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(ggplot2)
  library(pheatmap)
  library(ggfortify)
  library(vegan)
  library(rstatix)
})


####### Load data #######
setwd("~/OpenPedCan-Project-CNH")

# xCell scores can be compared between samples but not across cell types or cancer types per the documentation. 
# Therefore in this analysis, we'll rely on the results from quanTIseq deconvolution. 
quantiseq_output <- readRDS("analyses/immune-deconv/results/quantiseq_output.rds")
# confirm Medulloblastoma is in output
table(quantiseq_output$cancer_group) 

mb_output <- quantiseq_output %>% filter(cancer_group == "Medulloblastoma") %>% 
  mutate(molecular_subtype = stringr::str_replace(molecular_subtype, "^MB, ", "")) %>%
  mutate(molecular_subtype = factor(molecular_subtype, levels = c("WNT", "SHH", "Group3", "Group4")))

# tabulate medulloblastoma subtypes and deconvoluted cell types 
table(mb_output[, c("molecular_subtype", "cell_type")])

###### Per subtype stats summary #############
summary_tbl <- mb_output %>%
  group_by(molecular_subtype, cell_type) %>%
  summarise(
    n = n(),
    mean_fraction = mean(fraction, na.rm = TRUE),
    median_fraction = median(fraction, na.rm = TRUE),
    sd_fraction = sd(fraction, na.rm = TRUE),
    .groups = "drop"
  )

# For uncharacterized cells, their fractions are consistently high among all four subtypes (0.85 - 0.92).
# Therefore, it's likely safe to drop them and recalculate the fractions for other real cell types.
mb_filtered <- mb_output %>% filter(cell_type != "uncharacterized cell") %>% 
  group_by(Kids_First_Biospecimen_ID) %>%
  mutate(fraction_norm = fraction /sum(fraction, na.rm = TRUE)) %>% ungroup()

###### Visualizations of cell type distributions across medulloblastoma subtypes ##########
# Viz 1: PCA based on cell fractions across patients, colored by subtypes
mb_wide <- mb_filtered %>%
  select(Kids_First_Biospecimen_ID, molecular_subtype, cell_type, fraction_norm) %>%
  pivot_wider(names_from = cell_type, values_from = fraction_norm, values_fill = 0)
meta <- mb_wide %>% select(Kids_First_Biospecimen_ID, molecular_subtype)

# Perform PCA on cell-type fractions
pca_res <- prcomp(mb_wide %>% select(-Kids_First_Biospecimen_ID, -molecular_subtype),
                  scale. = TRUE)

viz1 <- autoplot(pca_res, data = meta, colour = 'molecular_subtype') +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))+
  labs(title = "PCA of immune cell composition")

viz1
ggsave("analyses/immune-deconv/results/medulloblastoma_pca.png",
       plot= viz1, width=8, height=6, dpi = 300)

# Viz 2: Correlation heatmap of cell types
mat  <- mb_wide %>%
  select(-Kids_First_Biospecimen_ID, -molecular_subtype) %>%
  as.data.frame()
rownames(mat) <- meta$Kids_First_Biospecimen_ID

# correlation of immune cell fractions between patients
cor_pat <- cor(t(as.matrix(mat)), method = "spearman", use = "pairwise.complete.obs")

# reorder metadata and correlation matrix to group by subtypes
subtype_levels <- c("WNT","SHH","Group3","Group4")
meta$molecular_subtype <- factor(meta$molecular_subtype, levels = subtype_levels)

order_idx <- order(meta$molecular_subtype)
meta_ordered <- meta[order_idx, ]
cor_pat <- cor_pat[meta_ordered$Kids_First_Biospecimen_ID, meta_ordered$Kids_First_Biospecimen_ID]

# set up annotation df and colors
ann <- data.frame(molecular_subtype = meta$molecular_subtype)
rownames(ann) <- meta$Kids_First_Biospecimen_ID
ann$molecular_subtype <- factor(ann$molecular_subtype, levels = subtype_levels)

subtype_cols <- RColorBrewer::brewer.pal(4, "Set2")
names(subtype_cols) <- subtype_levels
ann_colors <- list(molecular_subtype = subtype_cols)

ph <- pheatmap(
  cor_pat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  clustering_method = "average",
  annotation_row = ann,
  annotation_col = ann,
  annotation_colors = ann_colors,
  show_rownames = F,
  show_colnames = F,
  main = "Patient × patient correlation of immune composition (Spearman)",
  border_color = NA,
  annotation_names_row = FALSE
)
ph

png("analyses/immune-deconv/results/medulloblastoma_cor_heatmap.png", width = 2000, height = 1600, res = 300)
print(ph)
dev.off()

# Viz 3: Cell type compositions across mb subtypes and patients (boxplot + violinplot)
viz3 <- mb_filtered %>% 
  ggplot(aes(x=molecular_subtype, y=fraction_norm, fill = molecular_subtype)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.15, outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2", name="Molecular subtype")+
  facet_wrap(~cell_type, scales = "free_y") +
  labs(x = "Medulloblastoma subtype", y = "Proportion of cell types") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
viz3
ggsave("analyses/immune-deconv/results/medulloblastoma_violing_boxplot.png",
       plot = viz3, width=10, height=6, dpi = 300)


# Viz 4: Stacked barplot that summerizes each subtype for cleaner visualization
mb_filtered_medians <- mb_filtered %>% 
  group_by(molecular_subtype, cell_type) %>%
  summarise(median_fraction = median(fraction_norm), .groups = "drop") %>%
  group_by(molecular_subtype) %>%
  mutate(median_fraction = median_fraction / sum(median_fraction)) %>%
  ungroup()

viz4 <- mb_filtered_medians %>%
  ggplot(aes(x = molecular_subtype, y = median_fraction, fill = cell_type)) + 
  geom_col() + scale_fill_brewer(palette="Set3")+
  labs(x = "Medulloblastoma subtype", 
       y = "Proportion of cell types",
       fill = "Cell type") + theme_classic()
viz4
ggsave("analyses/immune-deconv/results/medulloblastoma_stacked_barplot.png",
       plot = viz4, width = 8, height = 6, dpi = 300)

####### Statistical Significance ########

# 1. Apply PERMANOVA to test if overall immune cell composition differs by MB subtype
## needed to install vegan package
X <- mb_wide %>%
  select(-Kids_First_Biospecimen_ID, -molecular_subtype) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

vegan::adonis2(X ~ molecular_subtype, data =meta, method = "bray", permutations = 999)

# 2. Kruskal-Wallis (non-parametric) tests to perform per-cell type global comparisons
kw_tbl <- mb_filtered %>%
  group_by(cell_type) %>%
  rstatix::kruskal_test(fraction_norm ~ molecular_subtype) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

# 3. Pairwise Wilcoxon for every cell type between mb subtypes
pw_tbl <- mb_filtered %>%
  group_by(cell_type) %>%
  pairwise_wilcox_test(
    fraction_norm ~ molecular_subtype,
    p.adjust.method = "none", 
    detailed = TRUE
  ) %>% ungroup() %>%
  mutate(p.adj.global = p.adjust(p, method = "BH")) %>%
  select(cell_type, group1, group2, n1, n2, p, p.adj.global) %>%
  add_significance("p.adj.global")

# wilcoxon effect sizes
# requires coin package
eff_tbl <- mb_filtered %>%
  group_by(cell_type) %>%
  wilcox_effsize(
    fraction_norm ~ molecular_subtype,
    paired = FALSE,
    comparisons = combn(levels(mb_filtered$molecular_subtype), 2, simplify = FALSE)
  ) %>%
  ungroup()

pairwise_stats <- pw_tbl %>%
  left_join(eff_tbl %>% select(cell_type, group1, group2, effsize, magnitude),  
            by = c("cell_type","group1","group2")) %>%
  arrange(cell_type, -effsize)

readr::write_csv(pairwise_stats, file = "analyses/immune-deconv/results/medulloblastoma_pairwise_stats.csv")

##### Additional visualizations after statistical tests ######
# viz 5: Heatmap of effect sizes with significance overlay
viz_5_df <- pairwise_stats %>%
  mutate(comp = paste(group1, "vs", group2))

viz5 <- ggplot(viz_5_df, aes(x = comp, y = cell_type, fill = effsize)) +
  geom_tile(color = "white") +
  geom_text(aes(label = p.adj.global.signif), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Wilcoxon effect sizes"
  ) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Subtype comparison", y = "Cell type")

viz5

ggsave("analyses/immune-deconv/results/medulloblastoma_pairwise_significance.png",
       plot = viz5, width = 8, height = 6, dpi = 300)

# viz 6: update viz 3 with significance
sig_pairs <- pairwise_stats %>%
  mutate(label = ifelse(p.adj.global.signif == "ns", "", p.adj.global.signif)) %>%
  filter(p.adj.global.signif != "ns") %>%
  arrange(cell_type, p.adj.global)

# Function to plot one cell type with significance labels
plot_celltype <- function(ct, yvar = "fraction_norm", outdir = "analyses/immune-deconv/results") {
  df_ct  <- mb_filtered %>% filter(cell_type == ct)
  sig_ct <- sig_pairs %>% filter(cell_type == ct) %>%
    mutate(molecular_subtype = NA)
  
  if (nrow(df_ct) == 0) return(invisible(NULL))
  
  y_base <- max(df_ct[[yvar]], na.rm = TRUE)
  if (!is.finite(y_base)) y_base <- 0.5
  if (nrow(sig_ct) > 0) {
    sig_ct <- sig_ct %>%
      mutate(y.position = y_base + seq_len(n()) * 0.05) %>%  
      mutate(xmin = group1, xmax = group2)                   
  }
  
  p <- ggplot(df_ct, aes(x = molecular_subtype, y = .data[[yvar]], fill = molecular_subtype)) +
    geom_violin(trim = TRUE, alpha = 0.8) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    scale_fill_brewer(palette = "Set2", name="Molecular subtype")+
    labs(
      title = paste0(ct, " across MB molecular subtypes"),
      x = "Molecular subtype",
      y = if (yvar == "fraction_norm") "Proportion of characterized immune cells" else yvar
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  if (nrow(sig_ct) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      data = sig_ct,
      label = "label",
      xmin  = "xmin", xmax  = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.4,
      size = 3, inherit.aes=FALSE
    )
  }
  return(p)
}

unique_cell_types <- unique(mb_filtered$cell_type)
plots <- lapply(unique_cell_types, plot_celltype)
pdf("analyses/immune-deconv/results/medulloblastoma_violin_boxplot_significance.pdf",
    width = 7.5, height = 5.5)
for (p in plots) {
  if (!is.null(p)) print(p)
}
dev.off()
