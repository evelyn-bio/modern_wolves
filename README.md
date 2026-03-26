# Modern Wolf Population Genomics Pipeline 
Snakemake workflow for large-scale population genomic analysis of modern wolf genomes. Designed to run on HPC clusters with SLURM job scheduling. 
## Overview 
This pipeline takes a multi-sample VCF file and runs a suite of population genomic analyses to characterise genetic structure, differentiation, and evolutionary relationships across samples. Developed for a large-scale study of European wolf population history and conservation genomics.
DOI: https://doi.org/10.64898/2026.03.20.712003
## Analyses 
**Population structure** PCA, Haplonet clustering
**Allele sharing** Dstatistics, f4-ratio statistics
**Admixture graph** Orientagraph
See config.yaml for input format and filtering parameters. 
## Requirements
snakemake/8.11.3 miniconda/24.5.0
## Usage
snakemake  --use-conda --use-env --workflow-profile slurm --rerun-incomplete
