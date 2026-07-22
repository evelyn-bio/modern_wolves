# Population genomic analysis of wolf genomes

A Snakemake workflow for population genomic analyses of imputed wolf genomes.
This workflow is part of a study on wolf population genetics. A preprint of this study can be found here:
https://doi.org/10.64898/2026.03.20.712003

The workflow performs quality control, relatedness filtering, population structure analyses, admixture inference, D-statistics, qpAdm modelling, phylogenetic inference, genetic diversity analyses, and haplotype-based analyses from a merged imputed VCF.

The imputation workflow used to generate the input VCF is available in the `imputation/` directory.

---

## Overview

This workflow includes:

* **Quality control**
  * Detects duplicate and closely related individuals using KING.
  * Produces filtered datasets for downstream analyses.
  * Filters variants by minor allele frequency.

* **Population structure**
  * Principal component analysis (smartPCA).
  * ADMIXTURE clustering 
  * HaploNet clustering

* **Gene flow**
  * D-statistics.
  * f4-ratio statistics.

* **Admixture graph**
  * OrientAGraph 

* **Genetic diversity**
  * Per-country nucleotide diversity (p).

---

## Workflow structure

```
Imputed VCF
      ¦
      ?
Quality control
      ¦
      ?
Population structure
      ¦
      +--------? PCA
      +--------? ADMIXTURE
      +--------? HaploNet
      +--------? D-statistics
      +--------? f4-ratio
      +--------? OrientAGraph
      +--------? p Diversity
```

## Repository structure

```text
popgen/
+-- Snakefile
+-- config.yaml
+-- gnomix shell scripts
+-- scripts/
+-- data/
+-- README.md
```

---

## Input data

The workflow requires:

* merged imputed VCF
* reference genome (FASTA)
* sample lists
* recombination maps
* metadata files

---

## Additional analyses

Local ancestry inference using Gnomix was performed using standalone shell scripts provided in this directory. These analyses are not integrated into the Snakemake workflow.

---

## License

This project is released under the MIT License. See the `LICENSE` file for details.

---

## Author

Evelyn Todd