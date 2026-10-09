# MPN to sAML: paired RNA-seq analysis
# Author: Enzo Marques Panepucci
# Run from the project root: Rscript --vanilla analysis.R
# In RStudio: open the project using the conda R installation, then source("analysis.R").
# The existing inputs/ directory is required. All settings are defined below.
# Outputs: Figures 1A-D and 2A-B, Figures S1A-B/S2/S3/S4A-B/S5, and Tables S1-S11.
# Checkpoint: load("checkpoints/analysis_objects.rda") restores named objects.
# Graphs are always rebuilt from data in the current session; old plot RDS files are not read.

# 1. Settings and packages -------------------------------------------------------
cfg <- list(
  seed = 2026L,
  reuse_checkpoints = TRUE,
  min_count = 10L,
  min_samples = 6L,
  remove_xy = TRUE,
  padj_cutoff = 0.05,
  lfc_cutoff = 0.5,
  results_alpha = 0.10,
  gsea_min_size_kegg = 15L,
  gsea_min_size_hallmark = 10L,
  gsea_max_size = 500L,
  gsea_eps = 1e-10,
  gsea_method = "hybrid",
  msigdb_version = "7.5.1",
  tsne_perplexity = 5L,
  png_dpi = 300L
)

# Restrict an activated conda session to the packages in that environment.
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (nzchar(conda_prefix)) {
  if (!startsWith(normalizePath(R.home()), normalizePath(conda_prefix))) {
    stop("This R session does not use the activated conda environment. Run its Rscript.")
  }
  .libPaths(.Library)
}
required <- c("DESeq2", "fgsea", "dplyr", "tidyr", "tibble", "ggplot2",
              "ggrepel", "pheatmap", "cowplot", "patchwork", "scales",
              "Rtsne", "yaml", "digest", "data.table", "matrixStats")
installed <- installed.packages()
missing <- setdiff(required, rownames(installed))
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       ". Create and activate the environment from environment.yml before running.")
}
for (pkg in intersect(required, loadedNamespaces())) {
  if (as.character(getNamespaceVersion(pkg)) != installed[pkg, "Version"]) {
    stop("Restart R before running: the loaded version of ", pkg,
         " differs from the version installed in this environment.")
  }
}
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})
# The remaining packages are called with package::function.
set.seed(cfg$seed)
for (folder in c("outputs/figures", "outputs/supplementary_figures",
                 "outputs/tables", "outputs/reports", "checkpoints")) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}
input_files <- c("inputs/count_matrix.csv", "inputs/coldata.csv",
                 "inputs/cibersortx_job40.csv", "inputs/gene_modules.csv",
                 "inputs/reference/ensembl_annotation_raw.tsv",
                 "inputs/reference/msigdb_7.5.1_ensembl.rds")
if (!all(file.exists(input_files))) {
  stop("Missing inputs: ", paste(input_files[!file.exists(input_files)], collapse = ", "))
}
input_checksums <- tools::md5sum(input_files)
package_versions <- as.data.frame(installed.packages()[, c("Package", "Version")])
package_versions <- package_versions[order(package_versions$Package), ]
rownames(package_versions) <- NULL
runtime_key <- digest::digest(list(R = R.version.string, packages = package_versions))
checkpoint_file <- "checkpoints/analysis_objects.rda"
checkpoint_env <- new.env(parent = emptyenv())
checkpoint_names <- character()
if (isTRUE(cfg$reuse_checkpoints) && file.exists(checkpoint_file)) {
  checkpoint_names <- load(checkpoint_file, envir = checkpoint_env)
}

theme_set(theme_classic(base_size = 12) +
            theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10, colour = "grey35"),
                  plot.tag = element_text(face = "bold", size = 18),
                  legend.title = element_text(face = "bold"),
                  plot.margin = margin(12, 14, 12, 12)))
group_colors <- c("Pre-leukaemic" = "#9ECAE1", "Post-leukaemic" = "#CB181D",
                  "NP first" = "#BDBDBD", "NP second" = "#525252")
cell_colors <- c("Monocyte" = "#4C78A8", "NK" = "#59A14F", "GMP" = "#76B7B2",
                 "Progenitor" = "#E6B84A", "Malignant" = "#D85759",
                 "Promonocyte" = "#AA78A7")

# 2. Inputs and sample quality control ----------------------------------------
counts_raw <- as.matrix(read.csv("inputs/count_matrix.csv", row.names = 1,
                                check.names = FALSE))
col_data <- read.csv("inputs/coldata.csv", row.names = 1, check.names = FALSE)
cibersort_raw <- read.csv("inputs/cibersortx_job40.csv", row.names = 1,
                         check.names = FALSE)
gene_modules <- read.csv("inputs/gene_modules.csv", stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(rownames(counts_raw)), !anyDuplicated(colnames(counts_raw)),
          !anyDuplicated(rownames(col_data)), !anyNA(counts_raw),
          all(is.finite(counts_raw)), all(counts_raw >= 0),
          all(counts_raw == floor(counts_raw)),
          all(rownames(col_data) %in% colnames(counts_raw)),
          all(col_data$Condition %in% c("pre", "pos")))

col_data$sample <- rownames(col_data)
col_data$Progression_Group <- ifelse(grepl("^NP", col_data$Patient_ID),
                                     "Non-progressor", "Progressor")
col_data$Sampling_point <- factor(dplyr::case_when(
  col_data$Progression_Group == "Progressor" & col_data$Condition == "pre" ~ "Pre-leukaemic",
  col_data$Progression_Group == "Progressor" & col_data$Condition == "pos" ~ "Post-leukaemic",
  col_data$Condition == "pre" ~ "NP first",
  TRUE ~ "NP second"
), levels = names(group_colors))
col_data <- col_data[order(col_data$Sampling_point, col_data$Patient_ID), ]
col_data$sample_label <- col_data$Patient_ID
pair_check <- table(col_data$Patient_ID, col_data$Condition)
stopifnot(all(pair_check == 1))
col_data_prog <- droplevels(col_data[col_data$Progression_Group == "Progressor", ])
col_data_prog$Patient_ID <- factor(col_data_prog$Patient_ID)
col_data_prog$Condition <- factor(col_data_prog$Condition, levels = c("pre", "pos"))
stopifnot(nrow(col_data_prog) == 24, nrow(col_data) == 34)

# All HTSeq technical categories start with __, including __too_low_aQual.
counts_clean <- counts_raw[!grepl("^__", rownames(counts_raw)), , drop = FALSE]
rownames(counts_clean) <- sub("\\..*$", "", rownames(counts_clean))
stopifnot(!anyDuplicated(rownames(counts_clean)))
sample_qc <- data.frame(sample = colnames(counts_clean),
                        gene_assigned_counts = colSums(counts_clean),
                        detected_genes = colSums(counts_clean > 0)) |>
  left_join(col_data, by = "sample") |>
  mutate(included = !is.na(Patient_ID),
         reason = ifelse(included, "Included", "No metadata; not analysed"))
write.csv(sample_qc, "outputs/tables/Table_S1_sample_QC.csv", row.names = FALSE)

# 3. Annotation, expression filter, and paired DESeq2 -------------------------
annotation_raw <- read.delim("inputs/reference/ensembl_annotation_raw.tsv",
                             stringsAsFactors = FALSE, check.names = FALSE)
names(annotation_raw) <- c("Ensembl_ID", "gene_symbol", "chromosome")
annotation_raw$gene_symbol[annotation_raw$gene_symbol == ""] <- NA_character_
xy_ids <- unique(annotation_raw$Ensembl_ID[annotation_raw$chromosome %in% c("X", "Y")])
annotation <- annotation_raw |>
  arrange(Ensembl_ID, is.na(gene_symbol), gene_symbol, chromosome) |>
  distinct(Ensembl_ID, .keep_all = TRUE)
# This key changes when inputs, model settings, R, or installed packages change.
analysis_key <- digest::digest(list(
  revision = "single-script-v1",
  files = input_checksums[c("inputs/count_matrix.csv", "inputs/coldata.csv",
                            "inputs/reference/ensembl_annotation_raw.tsv")],
  design = "~ Patient_ID + Condition",
  parameters = cfg[c("min_count", "min_samples", "remove_xy", "results_alpha",
                     "padj_cutoff", "lfc_cutoff")],
  runtime = runtime_key))
de_objects <- c("dds", "res", "res_df", "degs", "counts_filtered", "counts_prog",
                "dds_all", "vsd_all", "vsd_prog", "vst_matrix",
                "normalized_counts_all", "keep_expression", "keep_chromosome")
checkpoint_objects <- c("analysis_key", de_objects)
if (isTRUE(cfg$reuse_checkpoints) &&
    identical(checkpoint_env$analysis_key, analysis_key) &&
    all(de_objects %in% checkpoint_names)) {
  for (object in de_objects) {
    assign(object, get(object, envir = checkpoint_env), envir = environment())
  }
  message("Paired DESeq2 and VST objects restored from the .rda checkpoint.")
} else {
  message("Fitting paired DESeq2 and calculating VST matrices.")
  prog_counts_raw <- counts_clean[, rownames(col_data_prog), drop = FALSE]
  keep_expression <- rowSums(prog_counts_raw >= cfg$min_count) >= cfg$min_samples
  keep_chromosome <- !cfg$remove_xy | !rownames(counts_clean) %in% xy_ids
  counts_filtered <- counts_clean[keep_expression & keep_chromosome,
                                  rownames(col_data), drop = FALSE]
  counts_prog <- counts_filtered[, rownames(col_data_prog), drop = FALSE]
  dds <- DESeqDataSetFromMatrix(counts_prog, col_data_prog,
                               design = ~ Patient_ID + Condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("Condition", "pos", "pre"), alpha = cfg$results_alpha)
  res_df <- as.data.frame(res) |>
    tibble::rownames_to_column("Ensembl_ID") |>
    left_join(annotation, by = "Ensembl_ID") |>
    mutate(direction = case_when(
      !is.na(padj) & padj < cfg$padj_cutoff & log2FoldChange > cfg$lfc_cutoff ~ "Up",
      !is.na(padj) & padj < cfg$padj_cutoff & log2FoldChange < -cfg$lfc_cutoff ~ "Down",
      TRUE ~ "Not DEG")) |>
    arrange(padj, Ensembl_ID)
  degs <- filter(res_df, direction != "Not DEG")

  # Joint VST for visualisation; this does not modify the DESeq2 contrast above.
  dds_all <- DESeqDataSetFromMatrix(counts_filtered, col_data, design = ~ Sampling_point)
  dds_all <- estimateSizeFactors(dds_all)
  vsd_all <- varianceStabilizingTransformation(dds_all, blind = FALSE)
  vsd_prog <- varianceStabilizingTransformation(dds, blind = FALSE)
  vst_matrix <- assay(vsd_all)
  normalized_counts_all <- counts(dds_all, normalized = TRUE)
  save(list = checkpoint_objects, file = checkpoint_file,
       compress = "gzip", compression_level = 1, version = 3)
}
write.csv(res_df, "outputs/tables/Table_S2_DESeq2_all_genes.csv", row.names = FALSE)
write.csv(degs, "outputs/tables/Table_S3_DEGs.csv", row.names = FALSE)
annotation_audit <- data.frame(Ensembl_ID = rownames(counts_clean),
                              passes_expression = keep_expression,
                              passes_chromosome = keep_chromosome) |>
  left_join(annotation, by = "Ensembl_ID")
write.csv(annotation_audit, "outputs/tables/Table_S11_annotation_filter_audit.csv", row.names = FALSE)
filter_summary <- data.frame(step = c("Input rows", "Gene rows", "Expression filter",
                                      "Expression filter + exclusion of X/Y", "DEGs",
                                      "DEGs with symbols", "Unique DEG symbols"),
                             n = c(nrow(counts_raw), nrow(counts_clean), sum(keep_expression),
                                   nrow(counts_filtered), nrow(degs),
                                   sum(!is.na(degs$gene_symbol)),
                                   n_distinct(na.omit(degs$gene_symbol))))
n_significant <- sum(res_df$padj < cfg$padj_cutoff, na.rm = TRUE)
message(n_significant, " genes with adjusted P < ", cfg$padj_cutoff, "; ",
        nrow(degs), " also meet the fold-change threshold; ",
        sum(degs$direction == "Down"), " downregulated and ",
        sum(degs$direction == "Up"), " upregulated.")

# 4. Figure 1A: volcano -------------------------------------------------------
volcano_df <- res_df |>
  mutate(
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
    plot_category = case_when(
      !is.na(padj) & padj < cfg$padj_cutoff & log2FoldChange > cfg$lfc_cutoff ~ "Up",
      !is.na(padj) & padj < cfg$padj_cutoff & log2FoldChange < -cfg$lfc_cutoff ~ "Down",
      !is.na(padj) & padj < cfg$padj_cutoff & abs(log2FoldChange) <= cfg$lfc_cutoff ~ "Significant, small effect",
      TRUE ~ "Not significant"
    ),
    plot_category = factor(
      plot_category,
      levels = c("Down", "Significant, small effect", "Not significant", "Up")
    )
  )
label_genes <- volcano_df |>
  filter(plot_category %in% c("Down", "Up"), !is.na(gene_symbol), abs(log2FoldChange) >= 1) |>
  arrange(padj) |> distinct(gene_symbol, .keep_all = TRUE) |> slice_head(n = 20)
p_1A <- ggplot(filter(volcano_df, is.finite(neg_log10_padj)),
               aes(log2FoldChange, neg_log10_padj, colour = plot_category)) +
  geom_point(alpha = 0.65, size = 1.3) +
  geom_vline(xintercept = c(-cfg$lfc_cutoff, cfg$lfc_cutoff), linetype = 2, colour = "grey55") +
  geom_hline(yintercept = -log10(cfg$padj_cutoff), linetype = 2, colour = "grey55") +
  ggrepel::geom_text_repel(data = label_genes, aes(label = gene_symbol),
                          size = 3, max.overlaps = Inf, seed = cfg$seed,
                          box.padding = 0.45, show.legend = FALSE) +
  scale_colour_manual(
    values = c(Down = "#2166AC",
               "Significant, small effect" = "#7B3294",
               "Not significant" = "#C3C8CE",
               Up = "#B2182B"),
    drop = FALSE
  ) +
  labs(tag = "(A)", title = "Differential expression during leukaemic progression",
       subtitle = paste0("12 paired patients | ", sum(degs$direction == "Down"),
                          " downregulated; ", sum(degs$direction == "Up"), " upregulated"),
       x = expression(log[2]~fold~change~"(post / pre)"),
       y = expression(-log[10]~adjusted~italic(P)), colour = NULL) +
  theme(legend.position = "top")

# 5. Figure 1B and Figure S2: expression heatmaps -----------------------------
degs_symbols <- degs |> filter(!is.na(gene_symbol)) |>
  arrange(padj, Ensembl_ID) |> distinct(gene_symbol, .keep_all = TRUE)
module_membership <- gene_modules |>
  left_join(res_df |> filter(!is.na(gene_symbol)) |>
              arrange(padj, Ensembl_ID) |> distinct(gene_symbol, .keep_all = TRUE), by = "gene_symbol") |>
  mutate(in_figure_1B = !is.na(direction) & direction != "Not DEG")
write.csv(module_membership, "outputs/tables/Table_S9_heatmap_gene_membership.csv", row.names = FALSE)
module_degs <- filter(module_membership, in_figure_1B)
stopifnot(nrow(module_degs) > 1)
hm_matrix <- vst_matrix[module_degs$Ensembl_ID, rownames(col_data), drop = FALSE]
rownames(hm_matrix) <- module_degs$gene_symbol
hm_z <- t(scale(t(hm_matrix)))
stopifnot(all(is.finite(hm_z)))
module_short <- c("1. Myeloid identity", "2. Myeloid sensors / effectors",
                   "3. Cytokines / chemokines", "4. Adhesion / antigen presentation",
                   "5. Proliferation / DNA repair", "6. Context-dependent oncogenes",
                   "7. Survival / RNA metabolism")
names(module_short) <- unique(gene_modules$module)
module_colors <- setNames(c("#8C510A", "#01665E", "#542788", "#4D4D4D",
                            "#B2182B", "#762A83", "#1B7837"), module_short)
annotation_col <- data.frame(`Sampling point` = col_data$Sampling_point,
                             `Progression status` = col_data$Progression_Group,
                             row.names = rownames(col_data), check.names = FALSE)
annotation_row <- data.frame(`Biological function` = factor(module_short[module_degs$module],
                                                           levels = module_short),
                             row.names = module_degs$gene_symbol, check.names = FALSE)
ann_colors <- list(`Sampling point` = group_colors,
                   `Progression status` = c(Progressor = "#7F0000", `Non-progressor` = "#08519C"),
                   `Biological function` = module_colors)
gaps_col <- head(cumsum(rle(as.character(col_data$Sampling_point))$lengths), -1)
gaps_row <- head(cumsum(rle(module_degs$module)$lengths), -1)
heat_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
hm_bio <- pheatmap::pheatmap(hm_z, color = heat_colors, annotation_col = annotation_col,
                   annotation_row = annotation_row, annotation_colors = ann_colors,
                   gaps_col = gaps_col, gaps_row = gaps_row,
                   cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA,
                   annotation_names_row = FALSE,
                   labels_col = col_data$sample_label, fontsize = 9,
                   fontsize_row = 8, fontsize_col = 8, angle_col = 90,
                   main = "Selected DEGs organised by biological function", silent = TRUE)
# Keep the annotation legends in one lower, aligned column.
heatmap_1B <- hm_bio$gtable
legend_index <- which(heatmap_1B$layout$name == "legend")
annotation_legend_index <- which(heatmap_1B$layout$name == "annotation_legend")
z_legend <- heatmap_1B$grobs[[legend_index]]
annotation_legend <- heatmap_1B$grobs[[annotation_legend_index]]
legend_columns <- unique(heatmap_1B$layout$l[c(legend_index, annotation_legend_index)])
matrix_column <- heatmap_1B$layout$l[heatmap_1B$layout$name == "matrix"]
heatmap_1B$widths[matrix_column] <- heatmap_1B$widths[matrix_column] +
  sum(heatmap_1B$widths[legend_columns])
heatmap_1B$widths[legend_columns] <- grid::unit(0, "pt")
heatmap_1B$grobs[c(legend_index, annotation_legend_index)] <-
  rep(list(grid::nullGrob()), 2)
p_1B <- cowplot::ggdraw() +
  cowplot::draw_grob(heatmap_1B, x = 0, y = 0, width = 0.78, height = 0.95) +
  cowplot::draw_grob(annotation_legend, x = 0.80, y = 0.05, width = 0.20, height = 0.84) +
  cowplot::draw_grob(z_legend, x = 0.80, y = 0, width = 0.20, height = 0.59) +
  cowplot::draw_label("Row Z-score", x = 0.80, y = 0.61, hjust = 0,
                      fontface = "bold", size = 9) +
  cowplot::draw_label("(B)", x = 0.01, y = 0.995, hjust = 0, vjust = 1,
                      fontface = "bold", size = 18)

all_deg_matrix <- vst_matrix[degs_symbols$Ensembl_ID, rownames(col_data), drop = FALSE]
rownames(all_deg_matrix) <- degs_symbols$gene_symbol
all_deg_z <- t(scale(t(all_deg_matrix)))
hm_all <- pheatmap::pheatmap(all_deg_z, color = heat_colors, annotation_col = annotation_col,
                   annotation_colors = ann_colors, gaps_col = gaps_col,
                   cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE,
                   labels_col = col_data$sample_label, fontsize_col = 8,
                   border_color = NA, main = "All annotated DEGs | row Z-score", silent = TRUE)
p_S2 <- cowplot::ggdraw(hm_all$gtable)

# 6. Figures 1C-D: GSEA -------------------------------------------------------
pathways <- readRDS("inputs/reference/msigdb_7.5.1_ensembl.rds")
# The ranking uses all Ensembl IDs with a calculated P value, without padj/effect filtering
# and without conversion or deduplication by gene symbol.
rank_table <- res_df |>
  filter(is.finite(pvalue), is.finite(log2FoldChange)) |>
  arrange(Ensembl_ID) |>
  mutate(p_safe = pmax(pvalue, .Machine$double.xmin),
         confidence = sign(log2FoldChange) * -log10(p_safe),
         impact = log2FoldChange * -log10(p_safe),
         hybrid = sign(log2FoldChange) * case_when(
           p_safe < 0.05 ~ -log10(0.05) + abs(log2FoldChange),
           p_safe < 0.10 ~ -log10(p_safe) + abs(log2FoldChange),
           TRUE ~ -log10(p_safe)))
gsea_key <- digest::digest(list(
  analysis_key, rank_table,
  input_checksums["inputs/reference/msigdb_7.5.1_ensembl.rds"],
  cfg[c("gsea_min_size_kegg", "gsea_min_size_hallmark", "gsea_max_size",
        "gsea_eps", "seed")], runtime_key))
checkpoint_objects <- c(checkpoint_objects, "gsea_key", "gsea_results", "rank_table", "pathways")
if (isTRUE(cfg$reuse_checkpoints) && identical(checkpoint_env$gsea_key, gsea_key) &&
    "gsea_results" %in% checkpoint_names) {
  gsea_results <- checkpoint_env$gsea_results
  message("GSEA results restored from the .rda checkpoint.")
} else {
  message("Calculating KEGG and Hallmark GSEA for the three ranking metrics.")
  gsea_results <- list()
  for (collection in names(pathways)) {
    for (metric in c("confidence", "impact", "hybrid")) {
      ranks <- sort(setNames(rank_table[[metric]], rank_table$Ensembl_ID), decreasing = TRUE)
      min_size <- if (collection == "KEGG") cfg$gsea_min_size_kegg else cfg$gsea_min_size_hallmark
      set.seed(cfg$seed)
      gsea_results[[paste(collection, metric, sep = "_")]] <- as.data.frame(
        fgsea::fgseaMultilevel(pathways[[collection]], ranks, minSize = min_size,
                               maxSize = cfg$gsea_max_size, eps = cfg$gsea_eps,
                               scoreType = "std", nproc = 1))
    }
  }
  save(list = checkpoint_objects, file = checkpoint_file,
       compress = "gzip", compression_level = 1, version = 3)
}
stopifnot(cfg$gsea_method %in% c("hybrid", "confidence", "impact"))
gsea_kegg <- gsea_results[[paste0("KEGG_", cfg$gsea_method)]]
gsea_hallmark <- gsea_results[[paste0("Hallmark_", cfg$gsea_method)]]
gsea_plots <- list()
gsea_heights <- c()
for (collection in names(pathways)) {
  gsea_tab <- gsea_results[[paste(collection, cfg$gsea_method, sep = "_")]]
  gsea_export <- gsea_tab |> arrange(padj)
  gsea_export$leadingEdge <- vapply(gsea_export$leadingEdge, paste, character(1), collapse = ";")
  table_id <- if (collection == "KEGG") "S4" else "S5"
  write.csv(gsea_export, paste0("outputs/tables/Table_", table_id, "_GSEA_", collection, ".csv"), row.names = FALSE)
  shown <- gsea_tab |> filter(!is.na(padj), padj < cfg$padj_cutoff) |>
    arrange(NES) |>
    mutate(label = gsub("_", " ", sub("^(KEGG_|HALLMARK_)", "", pathway)),
           label = factor(label, levels = unique(label)),
           direction = ifelse(NES > 0, "Post-leukaemic", "Pre-leukaemic"))
  gsea_plots[[collection]] <- ggplot(shown, aes(NES, label, fill = direction)) +
    geom_col(width = 0.78) + geom_vline(xintercept = 0, colour = "grey40") +
    scale_fill_manual(values = c("Pre-leukaemic" = "#2166AC", "Post-leukaemic" = "#B2182B")) +
    labs(tag = if (collection == "KEGG") "(C)" else "(D)",
         title = paste(collection, "pathway enrichment"),
         subtitle = paste0(nrow(shown), " pathways | adjusted P < 0.05 | ", cfg$gsea_method, " ranking"),
         x = "Normalised enrichment score (NES)", y = NULL, fill = NULL) +
    theme(axis.text.y = element_text(size = 8.5), legend.position = "top",
          axis.line.y = element_blank(), axis.ticks.y = element_blank())
  gsea_heights[collection] <- max(5, 2.5 + nrow(shown) * 0.20)
}
p_1C <- gsea_plots$KEGG
p_1D <- gsea_plots$Hallmark

# 7. Figures S4A-B: comparison of gene-ranking metrics ------------------------
rank_comparison <- bind_rows(gsea_results, .id = "analysis") |>
  tidyr::separate(analysis, into = c("collection", "metric"), sep = "_", extra = "merge")
rank_comparison$leadingEdge <- vapply(rank_comparison$leadingEdge, paste, character(1), collapse = ";")
write.csv(rank_comparison, "outputs/tables/Table_S10_GSEA_ranking_comparison.csv", row.names = FALSE)
write.csv(rank_table, "outputs/tables/GSEA_gene_rankings.csv", row.names = FALSE)
upset_plots <- list()
for (collection in names(pathways)) {
  membership <- data.frame(pathway = names(pathways[[collection]]))
  for (metric in c("confidence", "impact", "hybrid")) {
    tab <- gsea_results[[paste(collection, metric, sep = "_")]]
    membership[[metric]] <- membership$pathway %in% tab$pathway[!is.na(tab$padj) & tab$padj < cfg$padj_cutoff]
  }
  membership$pattern <- paste0(as.integer(membership$confidence),
                                as.integer(membership$impact), as.integer(membership$hybrid))
  intersections <- membership |> filter(pattern != "000") |> count(pattern, name = "n") |>
    arrange(desc(n), pattern) |> mutate(set = factor(pattern, levels = pattern))
  matrix_dots <- tidyr::crossing(pattern = intersections$pattern,
                                 metric = c("confidence", "impact", "hybrid")) |>
    mutate(set = factor(pattern, levels = intersections$pattern),
           row = match(metric, c("confidence", "impact", "hybrid")),
           included = substring(pattern, row, row) == "1")
  p_top <- ggplot(intersections, aes(set, n)) + geom_col(fill = "#264653", width = 0.68) +
    geom_text(aes(label = n), vjust = -0.4, size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(title = paste(collection, "| intersections of significant pathways"),
         subtitle = "Adjusted P < 0.05 within each collection and ranking",
         x = NULL, y = "Pathways") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
  p_bottom <- ggplot(matrix_dots, aes(set, row)) +
    geom_point(colour = "#DFE3E7", size = 4) +
    geom_line(data = filter(matrix_dots, included), aes(group = pattern), colour = "#264653") +
    geom_point(data = filter(matrix_dots, included), colour = "#264653", size = 4) +
    scale_y_reverse(breaks = 1:3, labels = c("Confidence", "Impact", "Hybrid")) +
    labs(x = NULL, y = NULL) + theme_void() +
    theme(axis.text.y = element_text(size = 11, colour = "black"))
  upset_plots[[collection]] <- p_top / p_bottom + patchwork::plot_layout(heights = c(3, 1))
}
p_S4A <- upset_plots$KEGG
p_S4B <- upset_plots$Hallmark

# 8. Figures 2A-B: CIBERSORTx cell-state estimates ----------------------------
cell_states <- names(cell_colors)
stopifnot(all(rownames(col_data) %in% rownames(cibersort_raw)),
          !anyDuplicated(rownames(cibersort_raw)),
          all(cell_states %in% names(cibersort_raw)))
cibersort <- cibersort_raw[rownames(col_data), , drop = FALSE]
fractions <- as.matrix(cibersort[, cell_states])
stopifnot(!anyNA(fractions), all(fractions >= 0), all(fractions <= 1),
          max(abs(rowSums(fractions) - 1)) < 1e-6)
cibersort_table <- tibble::rownames_to_column(cibersort, "sample") |>
  left_join(col_data, by = "sample")
write.csv(cibersort_table, "outputs/tables/Table_S6_CIBERSORTx_fractions.csv", row.names = FALSE)
ciber_long <- cibersort_table |>
  pivot_longer(all_of(cell_states), names_to = "CellState", values_to = "Proportion") |>
  mutate(CellState = factor(CellState, levels = cell_states),
         sample = factor(sample, levels = rownames(col_data)))
ciber_summary <- ciber_long |> group_by(Sampling_point, CellState) |>
  summarise(n = n(), mean = mean(Proportion), median = median(Proportion),
            sd = sd(Proportion), .groups = "drop")
write.csv(ciber_summary, "outputs/tables/Table_S7_CIBERSORTx_group_summary.csv", row.names = FALSE)
p_2A <- ggplot(ciber_long, aes(sample, Proportion, fill = CellState)) +
  geom_col(width = 0.88) +
  facet_grid(~ Sampling_point, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cell_colors, drop = FALSE) +
  scale_x_discrete(labels = setNames(col_data$sample_label, rownames(col_data))) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  labs(tag = "(A)", title = "Estimated cell-state composition",
       subtitle = "12 progressor pairs and 5 non-progressor pairs",
       x = "Patient", y = "Estimated proportion", fill = "Cell state") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        legend.position = "bottom", strip.background = element_rect(fill = "#EEF1F4", colour = NA),
        strip.text = element_text(face = "bold"), panel.spacing.x = grid::unit(0.8, "lines"))
p_2B <- ggplot(ciber_summary, aes(CellState, mean, fill = Sampling_point)) +
  geom_col(position = position_dodge(0.84), width = 0.78, alpha = 0.80,
           colour = "grey35", linewidth = 0.2) +
  geom_point(data = ciber_long, aes(y = Proportion, group = Sampling_point),
             position = position_jitterdodge(jitter.width = 0.13, dodge.width = 0.84, seed = cfg$seed),
             shape = 21, colour = "grey30", size = 1.7, alpha = 0.85, stroke = 0.25) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, 0.07))) +
  labs(tag = "(B)", title = "Mean cell-state proportions across sampling points",
       subtitle = "Bars: group means | Points: individual samples | Descriptive comparison",
       x = NULL, y = "Estimated proportion", fill = "Sampling point") +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 11))

# 9. Figure S5: paired CIBERSORTx changes -------------------------------------
paired_ciber <- ciber_long |> filter(Progression_Group == "Progressor") |>
  select(Patient_ID, Condition, CellState, Proportion) |>
  pivot_wider(names_from = Condition, values_from = Proportion) |>
  arrange(CellState, Patient_ID)
paired_results <- list()
for (cell in cell_states) {
  dat <- filter(paired_ciber, CellState == cell)
  no_change <- all(dat$pos == dat$pre)
  p <- if (no_change) NA_real_ else wilcox.test(dat$pos, dat$pre,
                                              paired = TRUE, exact = FALSE, correct = TRUE)$p.value
  paired_results[[cell]] <- data.frame(CellState = cell, n_pairs = nrow(dat),
    mean_pre = mean(dat$pre), mean_post = mean(dat$pos),
    mean_change = mean(dat$pos - dat$pre), median_change = median(dat$pos - dat$pre),
    pvalue = p, status = if (no_change) "All paired differences are zero" else "Paired Wilcoxon")
}
paired_tests <- bind_rows(paired_results) |> mutate(padj = p.adjust(pvalue, "BH"))
write.csv(paired_tests, "outputs/tables/Table_S8_CIBERSORTx_paired_tests.csv", row.names = FALSE)
paired_labels <- paired_tests |>
  mutate(label = ifelse(is.na(padj), "No variation", paste0("BH P = ", format.pval(padj, digits = 2))))
p_S5 <- ggplot(filter(ciber_long, Progression_Group == "Progressor"),
                aes(Sampling_point, Proportion, group = Patient_ID)) +
  geom_line(colour = "grey65", linewidth = 0.5) +
  geom_point(aes(fill = Sampling_point), shape = 21, size = 2.4) +
  geom_text(data = paired_labels, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 3.3) +
  facet_wrap(~ CellState, ncol = 3) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1.07)) +
  labs(title = "Paired changes in estimated cell-state proportions",
       subtitle = "12 progressors | Wilcoxon signed-rank tests with BH adjustment",
       x = NULL, y = "Estimated proportion") + theme(legend.position = "none")

# 10. Figures S1A-B and S3: sample structure and library sizes ----------------
vst_prog <- assay(vsd_prog)[, rownames(col_data_prog), drop = FALSE]
top_var <- head(order(matrixStats::rowVars(vst_prog), decreasing = TRUE), 500)
pca <- prcomp(t(vst_prog[top_var, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca_df <- data.frame(sample = rownames(pca$x), pca$x[, 1:2]) |>
  left_join(col_data_prog, by = "sample")
pca_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
p_S1A <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_line(aes(group = Patient_ID), colour = "grey80") +
  geom_point(aes(fill = Sampling_point), shape = 21, size = 3.8) +
  ggrepel::geom_text_repel(aes(label = Patient_ID), size = 3, seed = cfg$seed, max.overlaps = Inf) +
  scale_fill_manual(values = group_colors) +
  labs(tag = "(A)", title = "Principal component analysis",
       subtitle = "500 most variable genes | VST | Progressors only",
       x = sprintf("PC1 (%.1f%%)", pca_var[1]), y = sprintf("PC2 (%.1f%%)", pca_var[2]), fill = NULL) +
  theme(legend.position = "bottom")
set.seed(cfg$seed)
tsne <- Rtsne::Rtsne(t(vst_prog), perplexity = cfg$tsne_perplexity,
                     pca = TRUE, initial_dims = min(50, ncol(vst_prog) - 1),
                     theta = 0.5, max_iter = 1000, num_threads = 1)
tsne_df <- data.frame(sample = colnames(vst_prog), tSNE1 = tsne$Y[, 1], tSNE2 = tsne$Y[, 2]) |>
  left_join(col_data_prog, by = "sample")
p_S1B <- ggplot(tsne_df, aes(tSNE1, tSNE2)) +
  geom_point(aes(fill = Sampling_point), shape = 21, size = 3.8) +
  ggrepel::geom_text_repel(aes(label = Patient_ID), size = 3, seed = cfg$seed, max.overlaps = Inf) +
  scale_fill_manual(values = group_colors) +
  labs(tag = "(B)", title = "t-SNE of paired progressor samples",
       subtitle = paste0("All filtered genes | VST | Perplexity = ", cfg$tsne_perplexity),
       x = "t-SNE 1", y = "t-SNE 2", fill = NULL) + theme(legend.position = "bottom")
p_S3 <- ggplot(sample_qc, aes(reorder(sample, as.numeric(sub("Amostra_", "", sample))),
                               gene_assigned_counts / 1e6, fill = Sampling_point)) +
  geom_col(width = 0.8) + scale_fill_manual(values = group_colors, na.value = "#E8B36A",
                                           na.translate = TRUE) +
  labs(title = "Gene-assigned counts per sample", subtitle = "Grey/orange unannotated sample is excluded from analysis",
       x = "Sample", y = "Counts (millions)", fill = "Sampling point") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        legend.position = "bottom")
write.csv(pca_df, "outputs/tables/PCA_coordinates.csv", row.names = FALSE)
write.csv(tsne_df, "outputs/tables/tSNE_coordinates.csv", row.names = FALSE)

# 11. Export all figures, tables, and named objects ---------------------------
plots <- list(Figure_1A_volcano = p_1A, Figure_1B_biological_heatmap = p_1B,
              Figure_1C_GSEA_KEGG = p_1C, Figure_1D_GSEA_Hallmark = p_1D,
              Figure_2A_CIBERSORTx_stacked = p_2A, Figure_2B_CIBERSORTx_means = p_2B,
              Figure_S1A_PCA = p_S1A, Figure_S1B_tSNE = p_S1B,
              Figure_S2_all_DEGs_heatmap = p_S2, Figure_S3_library_sizes = p_S3,
              Figure_S4A_GSEA_KEGG_intersections = p_S4A,
              Figure_S4B_GSEA_Hallmark_intersections = p_S4B,
              Figure_S5_CIBERSORTx_paired = p_S5)
dimensions <- data.frame(name = names(plots),
                          width = c(10, 13, 12, 12, 12, 12, 9, 9, 12, 12, 10, 10, 12),
                          height = c(7, 12, gsea_heights["KEGG"], gsea_heights["Hallmark"],
                                     5.6, 6.2, 6.5, 6.5, 10, 6, 6, 6, 7))
for (i in seq_along(plots)) {
  folder <- if (grepl("^Figure_S", names(plots)[i])) "supplementary_figures" else "figures"
  stem <- file.path("outputs", folder, names(plots)[i])
  ggsave(paste0(stem, ".pdf"), plots[[i]], width = dimensions$width[i],
         height = dimensions$height[i], device = cairo_pdf, limitsize = FALSE)
  ggsave(paste0(stem, ".png"), plots[[i]], width = dimensions$width[i],
         height = dimensions$height[i], dpi = cfg$png_dpi, bg = "white", limitsize = FALSE)
}
# One .rda file contains analysis data and individually named current-session plots.
checkpoint_objects <- unique(c(checkpoint_objects,
  "cfg", "runtime_key", "input_checksums", "package_versions", "col_data", "col_data_prog",
  "sample_qc", "annotation", "annotation_audit", "filter_summary", "n_significant", "degs_symbols",
  "group_colors", "cell_colors", "cell_states", "module_short", "module_colors",
  "gaps_col", "gaps_row",
  "gene_modules", "module_membership", "module_degs", "hm_z", "all_deg_z",
  "annotation_col", "annotation_row", "ann_colors", "heat_colors", "hm_bio", "hm_all",
  "volcano_df", "label_genes", "gsea_kegg", "gsea_hallmark", "rank_comparison",
  "cibersort_table", "ciber_long", "ciber_summary", "paired_ciber", "paired_tests",
  "vst_prog", "pca", "pca_df", "pca_var", "tsne", "tsne_df", "plots", "dimensions",
  "p_1A", "p_1B", "p_1C", "p_1D", "p_2A", "p_2B", "p_S1A", "p_S1B", "p_S2",
  "p_S3", "p_S4A", "p_S4B", "p_S5"))
save(list = checkpoint_objects, file = checkpoint_file,
     compress = "gzip", compression_level = 1, version = 3)
write.csv(filter_summary, "outputs/tables/Filtering_summary.csv", row.names = FALSE)
write.csv(dimensions, "outputs/reports/figure_manifest.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "outputs/reports/sessionInfo.txt")
write.csv(package_versions, "outputs/reports/package_versions.csv", row.names = FALSE)
yaml::write_yaml(cfg, "outputs/reports/parameters_used.yml")
input_paths <- list.files("inputs", recursive = TRUE, full.names = TRUE)
write.csv(data.frame(file = input_paths, md5 = unname(tools::md5sum(input_paths))),
          "outputs/reports/input_checksums.csv", row.names = FALSE)

message("Finished. Figures and tables are in outputs/; objects are in ", checkpoint_file)
rm(checkpoint_env)

# Explore or resize figures in RStudio using their named objects:
# plot(p_1A)
# plot(p_1B)
# plot(p_2A)
# ggplot2::ggsave("Figure_1B_resized.pdf", p_1B, width = 14, height = 15, units = "in")
# Restore the objects in a later session in the same conda environment:
# library(DESeq2)
# library(ggplot2)
# load("checkpoints/analysis_objects.rda")
# plot(p_1A)

# Optional: freeze the activated conda environment after a successful run.
# conda env export | sed '/^prefix:/d' > environment.yml
