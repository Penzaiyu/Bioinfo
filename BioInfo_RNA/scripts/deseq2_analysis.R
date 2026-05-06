suppressPackageStartupMessages({
  required_pkgs <- c("DESeq2", "ggplot2", "pheatmap")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      paste0(
        "Missing R packages: ", paste(missing_pkgs, collapse = ", "),
        "\nInstall them first."
      )
    )
  }

  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript deseq2_analysis.R sample_info.tsv counts_dir outdir count_column")
}

sample_info_file <- args[[1]]
counts_dir <- args[[2]]
outdir <- args[[3]]
count_column <- as.integer(args[[4]])

if (!count_column %in% c(2L, 3L, 4L)) {
  stop("count_column must be 2, 3, or 4.")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)

sample_info <- read.delim(sample_info_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
rownames(sample_info) <- sample_info$sample

read_star_counts <- function(sample_id) {
  infile <- file.path(counts_dir, paste0(sample_id, ".ReadsPerGene.out.tab"))
  if (!file.exists(infile)) {
    stop("Missing count file: ", infile)
  }

  dat <- read.delim(infile, header = FALSE, stringsAsFactors = FALSE)
  dat <- dat[-seq_len(4), c(1, count_column)]
  colnames(dat) <- c("gene_id", sample_id)
  dat
}

count_tables <- lapply(sample_info$sample, read_star_counts)
count_matrix <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_tables)
count_matrix[is.na(count_matrix)] <- 0
rownames(count_matrix) <- count_matrix$gene_id
count_matrix <- as.matrix(count_matrix[, -1, drop = FALSE])
storage.mode(count_matrix) <- "integer"

keep <- rowSums(count_matrix >= 10) >= 3
count_matrix <- count_matrix[keep, , drop = FALSE]

write.table(
  data.frame(gene_id = rownames(count_matrix), count_matrix, check.names = FALSE),
  file = file.path(outdir, "tables", "raw_count_matrix.filtered.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info[colnames(count_matrix), , drop = FALSE],
  design = ~ group
)

dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized = TRUE)
write.table(
  data.frame(gene_id = rownames(normalized_counts), normalized_counts, check.names = FALSE),
  file = file.path(outdir, "tables", "normalized_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

vsd <- vst(dds, blind = FALSE)
plot_pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(plot_pca_data, "percentVar"))

p <- ggplot(plot_pca_data, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, show.legend = FALSE) +
  labs(
    title = "RNA-seq PCA",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(outdir, "plots", "PCA_plot.pdf"), p, width = 7, height = 5)

gene_vars <- apply(assay(vsd), 1, var)
top_var_genes <- head(order(gene_vars, decreasing = TRUE), 50)
mat <- assay(vsd)[top_var_genes, ]
mat <- t(scale(t(mat)))
pdf(file.path(outdir, "plots", "heatmap_top50_variable_genes.pdf"), width = 8, height = 10)
pheatmap(mat, annotation_col = sample_info[, "group", drop = FALSE], fontsize_row = 6)
dev.off()

write_results <- function(dds_object, contrast_name, numerator, denominator) {
  res <- results(dds_object, contrast = c("group", numerator, denominator))
  res <- lfcShrink(dds_object, contrast = c("group", numerator, denominator), res = res, type = "normal")
  res_df <- as.data.frame(res[order(res$padj), ])
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]

  write.table(
    res_df,
    file.path(outdir, "tables", paste0(contrast_name, ".DESeq2.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  sig_df <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
  write.table(
    sig_df,
    file.path(outdir, "tables", paste0(contrast_name, ".sig.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  volcano <- subset(res_df, !is.na(padj) & !is.na(log2FoldChange))
  volcano$significant <- volcano$padj < 0.05 & abs(volcano$log2FoldChange) >= 1

  g <- ggplot(volcano, aes(log2FoldChange, -log10(padj), color = significant)) +
    geom_point(alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
    theme_bw(base_size = 12) +
    labs(title = contrast_name, x = "log2 fold change", y = "-log10 adjusted p-value")

  ggsave(file.path(outdir, "plots", paste0(contrast_name, ".volcano.pdf")), g, width = 7, height = 5)
}

write_results(dds, "groupB_vs_groupA", levels(sample_info$group)[2], levels(sample_info$group)[1])
