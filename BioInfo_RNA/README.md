# RNA-seq Differential Expression Pipeline

This repository contains a lightweight RNA-seq analysis workflow for paired-end bulk RNA-seq data.

## Included components

- `scripts/run_rnaseq_pipeline.sh`
  End-to-end shell pipeline for:
  - raw FASTQ quality control
  - adapter trimming
  - STAR genome indexing and alignment
  - gene-level counting with `STAR --quantMode GeneCounts`
  - downstream differential expression analysis

- `scripts/deseq2_analysis.R`
  DESeq2 workflow that produces:
  - filtered count matrix
  - normalized counts
  - PCA plot
  - top variable gene heatmap
  - volcano plots
  - differential expression tables

- `scripts/annotate_deseq_tables.R`
  Adds gene symbols and gene biotypes to DESeq2 result tables using a GTF annotation file.

- `scripts/make_gene_name_heatmap.R`
  Regenerates the top-variable-gene heatmap with gene symbols on the y-axis.

- `config/sample_info.template.tsv`
  Template sample sheet to be filled with your own FASTQ paths.

## Requirements

- Bash
- `fastqc`
- `fastp`
- `STAR`
- `Rscript`
- R packages:
  - `DESeq2`
  - `ggplot2`
  - `pheatmap`

## Usage

1. Edit `config/sample_info.template.tsv` with your sample names, groups, and FASTQ paths.
2. Update reference paths in `scripts/run_rnaseq_pipeline.sh`.
3. Run:

```bash
bash scripts/run_rnaseq_pipeline.sh
```

## Notes

- This archive intentionally excludes:
  - raw sequencing data
  - result files
  - server credentials and connection helpers
  - project-specific output paths

- You should adapt output directories and reference file paths for your own environment.
