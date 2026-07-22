## Imputation pipeline for wolf genomes

A Snakemake workflow for imputation of wolf genomes using a phased canid reference panel using the GLIMPSE1 pipeline.
This workflow is part of a study on wolf population genetics. A preprint of this study can be found here:
https://doi.org/10.64898/2026.03.20.712003

The canid reference panel is described in:

Bougiouri *et al.* (2025). *Proceedings of the National Academy of Sciences*.
https://doi.org/10.1073/pnas.2416980122

## Overview

This repository contains three Snakemake workflows:

* **Reference panel preparation (`refpanel.smk`)**

  * Extracts biallelic SNP sites from a phased reference panel.
  * Generates site lists for genotype calling.
  * Creates GLIMPSE chunk definitions.

* **Downsampling workflow (`downsample.smk`)**

  * Downsamples BAM files to multiple target coverages.
  * Calls variants at reference panel sites.
  * Performs genotype imputation with GLIMPSE.
  * Evaluates imputation accuracy across genotype probability thresholds.

* **Target panel imputation (`imputepanel.smk`)**

  * Calls variants from low-coverage wolf genomes.
  * Performs chromosome-wise genotype imputation using GLIMPSE.
  * Filters and phases imputed genotypes.
  * Produces a merged imputed panel suitable for downstream population genomic analyses.

## Workflow structure
  
- Reference panel preparation
- Downsampling (optional)
- Variant calling
- Imputation
- Filtering & phasing
- Final imputed panel

---

## Repository structure

```text
imputation/
├── Snakefile
├── config.yaml
├── rules/
│   ├── refpanel.smk
│   ├── downsample.smk
│   └── imputepanel.smk
├── envs/
├── scripts/
└── README.md
```

---

## Input data

The workflow requires:

* phased reference panel VCF(s)
* recombination map(s)
* reference genome (FASTA)
* BAM files for target samples
* text files listing BAM paths and sample names

The workflow can be run using snakemake --executor slurm --use-conda --use-env --workflow-profile slurm --dry-run


---

## License

This project is released under the MIT License. See the `LICENSE` file for details.

---

## Author

Evelyn Todd
