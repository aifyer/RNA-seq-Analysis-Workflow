# RNA-seq Analysis Workflow

I have developed analysis workflows for a variety of high-throughput sequencing and omics data types, including RNA-seq, spatial transcriptomics, spatial proteomics, single-cell RNA-seq (scRNA-seq), ATAC-seq, ChIP-seq, as well as microarray-based platforms such as ChIP-chip and DNA methylation arrays.

This repository contains a generic RNA-seq analysis workflow for read alignment, quality assessment, gene-level counting, and differential expression analysis.

## Files

- `script0.sh`: aligns paired-end FASTQ files with STAR, indexes BAM files with samtools, and runs selected RSeQC quality-control steps.
- `script1.R`: runs `featureCounts` to generate gene-level read counts from aligned BAM files.
- `script2.R`: performs normalization, PCA, differential expression testing, and visualization using edgeR and related R packages.

## Workflow

1. Prepare sample identifiers in `sample_ids.txt`.
2. Run STAR alignment and QC with `script0.sh`.
3. Prepare a list of BAM files in `bam_files.txt`.
4. Generate gene count data with `script1.R`.
5. Prepare sample metadata in `coldata.csv`.
6. Run downstream differential expression analysis with `script2.R`.

## Required Inputs

- Paired-end FASTQ files
- STAR genome index
- Gene annotation file in GTF format
- RSeQC reference annotation
- Sample metadata file
- BAM file list
- Gene annotation/name table

## Main Outputs

- Sorted and indexed BAM files
- RSeQC quality-control results
- `feature_counts.rda`
- `featureCounts.summary`
- Differential expression summary tables
- PCA plots
- Volcano plots
- Heatmaps

## Software Requirements

- STAR
- samtools
- RSeQC
- R
- Rsubread
- edgeR
- tidyverse-compatible R packages
- plotly
- htmlwidgets
- RColorBrewer

## Notes

The scripts are intended as reusable main workflow code. Local paths, sample names, group names, and environment-specific variables should be configured separately before use.
