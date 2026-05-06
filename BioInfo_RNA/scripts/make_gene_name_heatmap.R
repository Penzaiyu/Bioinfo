suppressPackageStartupMessages({
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript make_gene_name_heatmap.R <sample_info.tsv> <counts_dir> <gtf_file> <out_pdf>")
}

sample_info_file <- args[[1]]
counts_dir <- args[[2]]
gtf_file <- args[[3]]
out_pdf <- args[[4]]

extract_attr <- function(x, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  m <- regexec(pattern, x)
  regmatches(x, m) |>
    vapply(function(hit) if (length(hit) >= 2) hit[2] else NA_character_, character(1))
}

sample_info <- read.delim(sample_info_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
rownames(sample_info) <- sample_info$sample

read_star_counts <- function(sample_id) {
  infile <- file.path(counts_dir, paste0(sample_id, ".ReadsPerGene.out.tab"))
  dat <- read.delim(infile, header = FALSE, stringsAsFactors = FALSE)
  dat <- dat[-seq_len(4), c(1, 2)]
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

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info[colnames(count_matrix), , drop = FALSE],
  design = ~ group
)
dds <- DESeq2::estimateSizeFactors(dds)
vsd <- DESeq2::vst(dds, blind = FALSE)

gene_vars <- apply(SummarizedExperiment::assay(vsd), 1, var)
top_var_genes <- head(order(gene_vars, decreasing = TRUE), 50)
mat <- SummarizedExperiment::assay(vsd)[top_var_genes, ]
mat <- t(scale(t(mat)))

gtf <- read.delim(
  gtf_file,
  header = FALSE,
  comment.char = "#",
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE
)
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
genes <- gtf[gtf$feature == "gene", "attribute", drop = FALSE]
genes$gene_id <- extract_attr(genes$attribute, "gene_id")
genes$gene_name <- extract_attr(genes$attribute, "gene_name")
genes$gene_id_base <- sub("\\.[0-9]+$", "", genes$gene_id)
genes <- unique(genes[, c("gene_id", "gene_id_base", "gene_name")])

row_ids <- rownames(mat)
row_ids_base <- sub("\\.[0-9]+$", "", row_ids)
gene_name_map <- genes$gene_name[match(row_ids_base, genes$gene_id_base)]
display_names <- ifelse(is.na(gene_name_map) | gene_name_map == "", row_ids, gene_name_map)
display_names <- make.unique(display_names)
rownames(mat) <- display_names

pdf(out_pdf, width = 8, height = 10)
pheatmap(mat, annotation_col = sample_info[, "group", drop = FALSE], fontsize_row = 7)
dev.off()
