args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript annotate_deseq_tables.R <gtf_file> <results_dir> <out_summary_file>")
}

gtf_file <- args[[1]]
results_dir <- args[[2]]
out_summary_file <- args[[3]]

extract_attr <- function(x, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  m <- regexec(pattern, x)
  regmatches(x, m) |>
    vapply(function(hit) if (length(hit) >= 2) hit[2] else NA_character_, character(1))
}

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
genes$gene_type <- extract_attr(genes$attribute, "gene_type")
genes$gene_id_base <- sub("\\.[0-9]+$", "", genes$gene_id)
genes <- unique(genes[, c("gene_id", "gene_id_base", "gene_name", "gene_type")])

tables_dir <- file.path(results_dir, "tables")
summary_lines <- character()

for (f in list.files(tables_dir, pattern = "\\.DESeq2\\.tsv$", full.names = FALSE)) {
  path <- file.path(tables_dir, f)
  dat <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  dat$gene_id_base <- sub("\\.[0-9]+$", "", dat$gene_id)
  dat <- merge(dat, genes, by = "gene_id_base", all.x = TRUE, sort = FALSE)
  if ("gene_id.x" %in% colnames(dat)) {
    dat$gene_id <- dat$gene_id.x
  }
  dat <- dat[, c("gene_id", "gene_name", "gene_type", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

  out_path <- file.path(tables_dir, sub("\\.DESeq2\\.tsv$", ".annotated.tsv", f))
  write.table(dat, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

  sig <- subset(dat, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
  sig_path <- file.path(tables_dir, sub("\\.DESeq2\\.tsv$", ".sig.annotated.tsv", f))
  write.table(sig, sig_path, sep = "\t", quote = FALSE, row.names = FALSE)

  summary_lines <- c(summary_lines, paste0("## ", sub("\\.DESeq2\\.tsv$", "", f)), paste0("Significant genes: ", nrow(sig)), "")
}

writeLines(summary_lines, out_summary_file)
