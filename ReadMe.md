# systemsGWAS
### Formerly **multiomicGWAS**

> **systemsGWAS** is an integrative genome-wide association framework designed to identify the genetic architecture underlying observed phenotypes, latent biological traits, and systems-level interactions across genetics, environments, management, and multi-omic data.

⚠️ **Repository transition**

This repository was originally developed as **multiomicGWAS**. The software is currently transitioning to **systemsGWAS** to reflect its expanded scope beyond conventional multi-omic GWAS.

Repository:
https://github.com/bodeolukolu/multiomicGWAS

The repository name will be updated to **systemsGWAS** following completion of the accompanying methodological manuscripts.

---

# Introduction

Genome-wide association studies have traditionally focused on identifying loci controlling individual phenotypic traits. While highly successful, modern breeding increasingly requires understanding coordinated biological processes involving multiple interacting traits, environments, management practices, and host-associated omics.

**systemsGWAS** extends conventional GWAS by providing a unified framework capable of analysing:

- Conventional phenotypic traits
- Multi-trait analyses
- Secondary traits
- Host-associated microbiome and metagenome profiles
- Latent biological phenotypes
- Environment-specific trait architectures
- Future systems-level phenotypes such as the Developmental Blueprint Index (DBI)

The framework leverages **GWASpoly** for association analysis across ploidy levels (2x–8x) while providing automated data integration, relationship matrix construction, and flexible phenotypic modelling.

Relationship matrices are generated automatically, including:

- Genomic relationship matrices (GRM) using **AGHmatrix**
- Microbiome/metagenome kernels using Aitchison compositional distances
- Additional user-defined kernels (planned)

This architecture enables integrative systems-genetics analyses rather than treating genomic and omic datasets independently.

---

# Current Features

✔ GWAS across diploid through octoploid populations

✔ Automated genotype quality control

✔ Multiple GWASpoly genetic models

✔ Automatic genomic relationship matrix construction

✔ Integration of microbiome/metagenome abundance data

✔ Microbial taxa analysed as either:

- covariates
- independent phenotypes

✔ Automatic microbiome kernel construction

✔ Flexible phenotype processing

✔ Publication-quality figures

---

# Planned systemsGWAS Features

The long-term vision is to evolve systemsGWAS into a comprehensive systems genetics platform.

Planned capabilities include:

### Latent Phenotype GWAS

- Developmental Blueprint Index (DBI)
- user-defined latent phenotypes
- PCA-derived traits
- factor-analysis traits

---

### Multi-environment GWAS

- environment-specific GWAS
- stability GWAS
- plasticity GWAS
- G×E association analyses

---

### Systems Genetics

Integration of

- genomic data
- transcriptomics
- metabolomics
- microbiome
- environmental variables
- management variables

within a unified association framework.

---

### Crop Systems Integration

Future releases will support integration with crop systems modelling, allowing latent developmental phenotypes and environment-specific ideotypes to be analysed alongside genomic data.

---

# Software Roadmap

Version | Focus
--------|-------------------------------------------------
v1.x | multiomicGWAS
v2.x | systemsGWAS architecture
v2.x | latent phenotype framework
v2.x | multi-environment GWAS
v3.x | systems genetics
v4.x | crop systems integration (APSIM support)

---

# Installation

Clone the repository

```bash
git clone https://github.com/bodeolukolu/multiomicGWAS.git
```

After the repository rename:

```bash
git clone https://github.com/bodeolukolu/systemsGWAS.git
```

---

# Running systemsGWAS

Download

```
run_parameters_systemsGWAS.R
```

Edit the analysis parameters and execute the pipeline.

The parameter file controls:

- genotype input
- phenotype input
- covariates
- microbiome/metagenome analyses
- GWAS models
- filtering
- output directories
- plotting options

---

# Dependencies

Core packages include

- GWASpoly
- AGHmatrix
- sommer
- compositions
- mice
- qvalue

Data manipulation

- data.table
- dplyr
- stringr
- reshape2
- zoo

Visualisation

- ggplot2
- qqplotr
- heatmaply
- GGally

Statistics

- ppcor

Additional packages may be required for optional modules.

---

# Publications

Applications of systemsGWAS include:

*(To be updated following publication.)*

---

# Citation

If you use systemsGWAS in your research, please cite the accompanying publication once available.

Until then, please reference the GitHub repository:

https://github.com/bodeolukolu/multiomicGWAS

---

# Contributing

Bug reports, feature requests, and pull requests are welcome.

Suggestions for extending systemsGWAS to additional omic data types, latent phenotype analyses, or systems genetics applications are encouraged.

---

# Contact

Dr. Bode A. Olukolu

Science Leader – Crop Genetics

Department of Primary Industries

Queensland, Australia

Email:
bolukolu@utk.edu

---

# Versioning

This project follows Semantic Versioning:

https://semver.org/

---

# License

Apache License 2.0

https://github.com/bodeolukolu/multiomicGWAS/blob/master/LICENSE
