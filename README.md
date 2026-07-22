# Modern wolves

A reproducible workflow for genotype imputation and population genomic analyses of modern wolf genomes.

This repository contains the analysis pipelines accompanying:

> **Todd, E. et al.** (2026). *The population structure and genetic health of European wolves*

Preprint: https://doi.org/10.64898/2026.03.20.712003

---

## Repository overview

The project is divided into two independent workflows.

### `imputation/`

A Snakemake workflow for genotype imputation using a phased canid reference panel.

The workflow includes:

- reference panel preparation
- low-coverage genotype calling
- GLIMPSE imputation
- genotype probability filtering
- imputation accuracy assessment

See `imputation/README.md` for details.

---

### `popgen/`

A Snakemake workflow for downstream population genomic analyses of the imputed dataset.

The workflow includes:

- quality control and relatedness filtering
- principal component analysis (PCA)
- ADMIXTURE
- HaploNet
- D-statistics
- qpAdm and f4-ratio analyses
- OrientAGraph
- nucleotide diversity analyses

Gnomix local ancestry inference analyses are provided as standalone shell scripts.

See `popgen/README.md` for details.

---

## Citation

If you use this workflow, please cite:

E. T. Todd et al., The population structure and genetic health of European wolves. bioRxiv 10.64898/2026.03.20.712003 (2026).


---

## License

This repository is released under the MIT License. See the `LICENSE` file for details.

---

## Author

Evelyn Todd