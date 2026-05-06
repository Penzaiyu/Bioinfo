#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/path/to/project"
OUTDIR="${PROJECT_ROOT}/RNA_R"
SAMPLE_INFO="${PROJECT_ROOT}/config/sample_info.tsv"

# Update these paths for your environment.
GENOME_FASTA="/path/to/reference/GRCh38.primary_assembly.genome.fa"
GTF_FILE="/path/to/reference/gencode.annotation.gtf"
STAR_INDEX_DIR="${OUTDIR}/reference/star_index"

# STAR GeneCounts columns:
# 2 = unstranded, 3 = forward-stranded, 4 = reverse-stranded
COUNT_COLUMN="${COUNT_COLUMN:-2}"
THREADS="${THREADS:-8}"

for cmd in fastp fastqc STAR Rscript awk sed; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "Missing required command: ${cmd}" >&2
        exit 1
    fi
done

if [[ ! -f "${SAMPLE_INFO}" ]]; then
    echo "Missing sample sheet: ${SAMPLE_INFO}" >&2
    exit 1
fi

if [[ ! -f "${GENOME_FASTA}" ]]; then
    echo "Missing genome FASTA: ${GENOME_FASTA}" >&2
    exit 1
fi

if [[ ! -f "${GTF_FILE}" ]]; then
    echo "Missing annotation GTF: ${GTF_FILE}" >&2
    exit 1
fi

mkdir -p \
    "${OUTDIR}/logs" \
    "${OUTDIR}/qc/fastp" \
    "${OUTDIR}/qc/fastqc_raw" \
    "${OUTDIR}/qc/fastqc_trimmed" \
    "${OUTDIR}/trimmed" \
    "${OUTDIR}/align" \
    "${OUTDIR}/counts" \
    "${OUTDIR}/results" \
    "${OUTDIR}/reference"

if [[ ! -s "${STAR_INDEX_DIR}/Genome" ]]; then
    mkdir -p "${STAR_INDEX_DIR}"
    STAR \
        --runThreadN "${THREADS}" \
        --runMode genomeGenerate \
        --genomeDir "${STAR_INDEX_DIR}" \
        --genomeFastaFiles "${GENOME_FASTA}" \
        --sjdbGTFfile "${GTF_FILE}" \
        --sjdbOverhang 149 \
        2>&1 | tee "${OUTDIR}/logs/star_genome_generate.log"
fi

tail -n +2 "${SAMPLE_INFO}" | while IFS=$'\t' read -r sample group fq1 fq2; do
    sample_dir="${OUTDIR}/align/${sample}"
    trim_r1="${OUTDIR}/trimmed/${sample}_R1.trimmed.fq.gz"
    trim_r2="${OUTDIR}/trimmed/${sample}_R2.trimmed.fq.gz"

    mkdir -p "${sample_dir}"

    fastqc -t "${THREADS}" -o "${OUTDIR}/qc/fastqc_raw" "${fq1}" "${fq2}" \
        2>&1 | tee "${OUTDIR}/logs/${sample}.fastqc_raw.log"

    fastp \
        --thread "${THREADS}" \
        --in1 "${fq1}" \
        --in2 "${fq2}" \
        --out1 "${trim_r1}" \
        --out2 "${trim_r2}" \
        --html "${OUTDIR}/qc/fastp/${sample}.fastp.html" \
        --json "${OUTDIR}/qc/fastp/${sample}.fastp.json" \
        2>&1 | tee "${OUTDIR}/logs/${sample}.fastp.log"

    fastqc -t "${THREADS}" -o "${OUTDIR}/qc/fastqc_trimmed" "${trim_r1}" "${trim_r2}" \
        2>&1 | tee "${OUTDIR}/logs/${sample}.fastqc_trimmed.log"

    STAR \
        --runThreadN "${THREADS}" \
        --genomeDir "${STAR_INDEX_DIR}" \
        --readFilesIn "${trim_r1}" "${trim_r2}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${sample_dir}/${sample}." \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        2>&1 | tee "${OUTDIR}/logs/${sample}.star.log"

    cp "${sample_dir}/${sample}.ReadsPerGene.out.tab" "${OUTDIR}/counts/${sample}.ReadsPerGene.out.tab"
done

Rscript "$(dirname "$0")/deseq2_analysis.R" \
    "${SAMPLE_INFO}" \
    "${OUTDIR}/counts" \
    "${OUTDIR}/results" \
    "${COUNT_COLUMN}" \
    2>&1 | tee "${OUTDIR}/logs/deseq2_analysis.log"
