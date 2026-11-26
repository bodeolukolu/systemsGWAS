<p align="center">
  <img src="misc/systemsGWAS_logo.png" alt="systemsGWAS logo" width="620">
</p>

#### systemsGWAS

**Genome-wide association analysis of measured phenotypes, environmental/geographic variables, microbiome and phenomic features, and composite biological traits.** 

**systemsGWAS** is an R pipeline that combines marker dosage data with measured phenotypes, environmental variables, or traits constructed from microbiome/metagenome and phenomic features (including numeric measurements extracted from hyperspectral or drone imagery). It supports single-trait and structure-calibrated environmental GWAS (sgwas-env), joint multi-trait GWAS of two or more phenotype columns, GWAS of individual features, unsupervised composites, phenotype-informed sPLS proxies, joint measured-trait–feature-score association, and SNP × feature-score interactions. Phenomic feature tables can replace microbiome tables, or both can be combined after source-specific preprocessing, with feature centring and scaling performed internally for PCA and sPLS. 

Legacy association analyses use **GWASpoly** and **AGHmatrix** (gwaspoly-marginal) for single-trait GWAS. Single-trait, environmental, and integrated association analyses (sgwas-marginal, sgwas-env, sgwas-multitrait, sgwas-joint, sgwas-interaction, sgwas-all) are implemented using a dedicated, high-speed Base-R spectral REML/GLS engine (**sgwas**) supporting both mixed linear models (MLM) and fixed-effect ordinary least squares (GLM), top-PC kinship deflation, genomic control calibration, and additive and dosage-threshold dominance models across diploid, tetraploid, hexaploid, and octoploid dosage files (2x, 4x, 6x, 8x). 

The existing metag_data_species and microbial_* argument names are used for both microbiome and phenomic inputs. The pipeline reads extracted numeric features, not raw image files. See [phenomic input requirements](#phenomic-features-from-hyperspectral-or-drone-imaging) for transformations and the current nonnegative-input requirement. 

---

## What can I analyse?

| Analysis | GWAS outcome | Does constructing the outcome use a target phenotype such as FER or DTA? |
|---|---|---|
| Conventional GWAS | A measured trait, such as Fusarium ear rot (FER) or days to anthesis (DTA) | Uses that measured trait directly |
| Multi-trait phenotype GWAS | Two or more measured traits from `phenofile`, fitted together | Uses the selected phenotypes directly; no microbiome required |
| Single-feature GWAS | One microbial abundance or phenomic measurement, such as a spectral band or canopy-cover estimate | No |
| Summed-abundance GWAS | Total abundance of a selected group of species | No |
| PCA-composite GWAS | PC1, PC2, or PC3 of selected microbiome or phenomic features, mapped separately | No |
| sPLS-proxy GWAS | Weighted microbiome or phenomic component trained against a target phenotype | Yes, when learning the weights |
| Saved-weight GWAS | A numeric score calculated using previously exported weights | No new target observations are required; the original weights retain their training history |
| Joint trait–score GWAS | A measured trait and one microbial/phenomic score as two outcomes | Depends on the score constructor; unsupervised scores are preferred |
| SNP × score GWAS | Measured trait; the SNP effect can depend on the microbial/phenomic score | Depends on the score constructor; both main effects are retained |

Current capabilities also include marker filtering, additive and dominance effect models (`Add`, `Dom`), numeric phenotype covariates, optional leave-one-chromosome-out (LOCO) kinship, association tables, and Manhattan, QTL-profile, and QQ plots. A legacy microbiome-adjustment pathway is retained for existing analyses.

**Scope:** Set `association_mode = "sgwas-multitrait"` to fit all selected `trait_names` together in one multivariate mixed model. See [multi-trait GWAS](#multi-trait-gwas-from-phenofile). In `association_mode = "sgwas-marginal"` (default for phenotype-only runs) or `"gwaspoly-marginal"`, multiple `trait_names` request separate single-trait scans. When `microbial_features` or `microbial_weights` is provided and `association_mode` is omitted, `systemsGWAS()` automatically defaults to `"sgwas-all"`, pairing each measured trait with the selected microbial/phenomic score in bivariate and interaction mixed models.

---

## Interpretation of association mode outputs

| Output Directory / Test | What it tests | Test df |
|---|---|---|
| `Trait_marginal` | Is the SNP associated with the measured trait, analysing that trait alone? The microbial score is not included. | 1 |
| `Microbial_marginal` | Is the SNP associated with the microbial score, analysing that score alone? The measured trait is not included. | 1 |
| `sgwas-env` | Is the SNP associated with an environmental/geographic variable, using optional top-PC GRM deflation, `GLM`, and genomic control calibration to balance population-structure confounding? | 1 |
| `sgwas-multitrait` | Are the SNP effects across all selected phenotypes jointly different from zero , alongside per-trait conditional and single-trait tests? | (omnibus) / 1 (per-trait) |
| `Trait_joint` | Is the SNP’s measured-trait effect significant when the trait and microbial score are fitted together as two correlated outcomes? | 1 |
| `Microbial_joint` | Is the SNP’s microbial-score effect significant in that same bivariate model? | 1 |
| `Bivariate_joint` | Are the SNP effects on the measured trait and microbial score jointly different from zero? Either outcome or both can drive significance. | 2 |
| `SNP_at_mean` | Is the SNP associated with the measured trait at the average microbial score, within the interaction model? | 1 |
| `Interaction` | Does the SNP’s effect on the measured trait change with the microbial score, after retaining both main effects? | 1 |
| `SNP_joint` | Does the SNP contribute to the measured trait through its main effect, its interaction with the microbial score, or both? | 2 |

---

## Installation

Clone the repository:

```sh
git clone https://github.com/bodeolukolu/systemsGWAS.git
cd systemsGWAS
```

Install the dependencies in R before running a full GWAS:

```r
cran_packages <- c(
  "devtools", "BiocManager", "AGHmatrix", "qqplotr", "corrplot", "ggcorrplot",
  "ggplot2", "plyr", "dplyr", "tidyr", "car", "MASS", "Hmisc", "data.table",
  "stringr", "heatmaply", "ppcor", "zoo", "GGally", "reshape2", "compositions",
  "sommer", "mice"
)
missing_packages <- cran_packages[
  !vapply(cran_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("GWASpoly", quietly = TRUE)) {
  devtools::install_github("jendelman/GWASpoly", build_vignettes = FALSE)
}
bioc_packages <- c("mixOmics", "qvalue")
missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
```

Use compatible R and Bioconductor releases. The Base-R `sgwas` engine (`sgwas-marginal`, `sgwas-multitrait`, `sgwas-joint`, `sgwas-interaction`, `sgwas-all`) does not require loading the legacy **GWASpoly** dependency set; sPLS score construction requires **mixOmics**, and verbose plotting uses **ggplot2**.

For **preparation only** (`run_gwas = FALSE`), independent composites and saved-weight application use Base R; sPLS additionally needs **mixOmics**. Genotypes and GWAS packages are not required for preparation.

---

## Input files

Use tab-separated files with headers. Sample identifiers must match exactly across the inputs used in an analysis. The legacy `gwaspoly-marginal` engine requires at least **1,000 markers passing its relationship-matrix filters**. The `sgwas` engine has no artificial 1,000-marker cutoff, though reliable genomic kinship still requires adequate genome-wide marker coverage.

### Marker dosage data

Provide a separate file for each requested ploidy through `genofile_2x`, `genofile_4x`, `genofile_6x`, or `genofile_8x`. The first five columns are `SNP`, `CHROM`, `POS`, `REF`, and `ALT`; subsequent columns are sample IDs.

```text
SNP	CHROM	POS	REF	ALT	Line01	Line02	Line03
SNP_001	1	10500	A	G	0	1	2
SNP_002	1	21900	C	T	2	0	1
SNP_003	2	31700	G	A	1	NA	0
```

Dosages range from 0 to the relevant ploidy; use `NA` for missing genotypes. Selecting a ploidy does not convert another dosage file into that ploidy.

### Phenotypes and numeric covariates

The first column contains sample IDs; subsequent columns contain numeric traits or covariates. Use `NA` for missing measurements. Duplicate phenotype IDs are averaged per column.

```text
Plant_ID	FER	DTA	Inoc	NI
Line01	12.4	65	16.7	8.1
Line02	20.8	68	27.3	14.3
Line03	NA	63	NA	10.2
```

`phenofile` is required for conventional trait GWAS, integrated joint/interaction GWAS, sPLS training, or requested phenotype covariates. It is unnecessary for independent microbiome/phenomic outcomes or saved-weight application without covariates.

### Microbial abundance data

For the microbial workflow, `metag_data_species` accepts either orientation below.

**Samples in rows:** the ID column can be `Plant_ID`, `SampleID`, `Sample_ID`, or `Taxa`; all remaining columns are numeric features.

```text
Plant_ID	Fusarium_verticillioides	Fusarium_proliferatum	Other_species
Line01	120	30	850
Line02	80	60	860
Line03	150	10	840
```

**Species in rows:** identify features with a `species` or `Feature` column and use sample IDs as the abundance columns. Recognised taxonomy columns (`genus`, `family`, `order`, `class`, `phylum`, `kingdom`, `domain`) are ignored.

```text
species	genus	Line01	Line02	Line03
Fusarium_verticillioides	Fusarium	120	80	150
Fusarium_proliferatum	Fusarium	30	60	10
Other_species	Other	850	860	840
```

Abundances must be finite and nonnegative, with unique sample IDs and feature names. Pass the complete intended reference microbiome even when selecting a genus: relative abundance and CLR use all input features as their reference.

### Phenomic features from hyperspectral or drone imaging

Supply a feature table through **the same `metag_data_species` argument**, for example `metag_data_species = "phenomic_features.tsv"`. Use one sample per row with `Plant_ID` followed by numeric measurements, or the transposed layout with a `Feature` column.

```text
Plant_ID	HS_450nm	HS_550nm	HS_650nm	HS_750nm	HS_850nm	UAV_canopy_cover
Line01	0.08	0.12	0.09	0.35	0.48	0.82
Line02	0.10	0.15	0.11	0.31	0.42	0.74
Line03	0.07	0.11	0.08	0.39	0.51	0.89
```

For phenomic tables, explicitly set:

```r
microbial_transform = "none",  # preserve the supplied measurement scale
microbial_prevalence = 0       # disable the microbial presence/absence filter
```

**Current input constraint:** The feature reader requires all values to be finite and nonnegative, even with `microbial_transform = "none"`. Signed indices or pre-standardized matrices can be constructed externally and mapped via `phenofile`.

---

## Quick start: conventional trait GWAS (`sgwas-marginal` default)

By default (`association_mode = "sgwas-marginal"`), phenotype-only runs use the fast Base-R REML/GLS engine (`sgwas`):

```r
source("R/systemsGWAS.R")
data_dir <- "/path/to/data"

systemsGWAS(
  wdir = data_dir,
  projname = "DTA_GWAS",
  ploidy_levels = 2,
  genofile_2x = "genotypes.tsv",
  phenofile = "traits.tsv",
  trait_names = "DTA",
  covariate_pheno = NULL,
  method = "MLM",
  model_effect = c("Add", "Dom"),
  maf = 0.02,
  cores = 1
)
```

To run the legacy **GWASpoly** engine instead, set `association_mode = "gwaspoly-marginal"`.

---

## Environmental GWAS (`sgwas-env`) and Integrated GLM Support

Geographic and environmental traits (such as `Env_PC1`, `Env_PC2`, latitude, longitude, or elevation) often correlate strongly with macro-population structure. Standard mixed linear models (`MLM`) with a full genomic relationship matrix (`GRM`) can over-correct for this structure and mask true local-adaptation loci (false negatives), whereas omitting the GRM entirely (`GLM`) can inflate false positives.

`systemsGWAS` addresses this trade-off with **`association_mode = "sgwas-env"`** and native **`method = "GLM"`** support in the base-R REML/GLS engine:

### Key Features
* **Top-PC GRM Deflation (`env_deflate_pcs`):** Removes the leading \(k\) genetic eigenvectors (which capture macro-geographic clines) from the GRM prior to REML variance estimation while retaining mid- and lower-order eigenvectors to control cryptic relatedness and family structure.
* **Genomic Control Calibration (`env_gc = TRUE`):** Automatically computes the empirical genomic inflation factor (\(\lambda_{GC}\)) per trait and marker model and exports calibrated p-values (`P_Trait_GC`) and Benjamini–Hochberg/Storey FDR q-values (`P_Trait_GC_BH`) alongside uncalibrated statistics in `ASSOCIATIONS.tsv`.
* **Integrated `GLM` Mode (`method = "GLM"`):** Runs fixed-effect Ordinary Least Squares (OLS) scans without a GRM inside the fast integrated engine (`sgwas-marginal`, `sgwas-env`, `sgwas-joint`, `sgwas-interaction`), optionally adjusting for fixed covariates (`covariate_pheno`) and leading genetic PCs (`association_n_pcs`).

### Example Usage

```r
# 1. Recommended EnvGWAS: LOCO MLM with top 1 geographic eigenvector deflated + Genomic Control
systemsGWAS(
  wdir = "./",
  projname = "Chickpea_EnvGWAS",
  ploidy_levels = 2,
  genofile_2x = "genotypes.tsv",
  phenofile = "env_pcs.tsv",
  trait_names = c("Env_PC1", "Env_PC2"),
  association_mode = "sgwas-env",
  method = "MLM",
  LOCO = TRUE,
  env_deflate_pcs = 1L,
  env_gc = TRUE
)

# 2. Integrated GLM with fixed genetic PCs (no GRM)
systemsGWAS(
  wdir = "./",
  projname = "Chickpea_EnvGLM",
  ploidy_levels = 2,
  genofile_2x = "genotypes.tsv",
  phenofile = "env_pcs.tsv",
  trait_names = c("Env_PC1", "Env_PC2"),
  association_mode = "sgwas-env",
  method = "GLM",
  association_n_pcs = 2L,
  env_gc = TRUE
)
```
---

## Multi-trait GWAS from phenofile

Use `association_mode = "sgwas-multitrait"` with two or more numeric phenotype columns:

```r
systemsGWAS(
  wdir = data_dir,
  projname = "NI_Inoc_DTA_MT",
  ploidy_levels = 2,
  genofile_2x = "genotypes.tsv",
  phenofile = "traits.tsv",
  trait_names = c("NI", "Inoc", "DTA"),
  association_mode = "sgwas-multitrait",
  multitrait_missing = "impute",          # Or "complete_case".
  multitrait_missing_threshold = 0.30,    # >30% overall OR in any trait: drop incomplete rows.
  method = "MLM",
  model_effect = "Add",
  covariate_pheno = NULL,
  covariate_metag = FALSE,
  LOCO = TRUE,
  association_maxit = 300L,
  association_plot_formats = c("pdf", "tiff")
)
```

The selected traits are fitted together using unstructured genetic and residual
trait covariance matrices. The mode exports `Trait_marginal`, `Trait_joint`, and
`Multitrait_joint` tests. The omnibus has one numerator degree of freedom per
selected trait and can be significant because of only one trait.

By default, `multitrait_missing = "impute"` retains partially observed samples and
exports conditional phenotype imputations with standard errors. Joint association
tests integrate missing outcomes out of the likelihood, so imputed values do not
count as extra observed measurements. Marginal baselines use observed values of
each trait only. Set `multitrait_missing = "complete_case"` to drop incomplete
phenotype rows explicitly.

Missingness is reported overall (missing phenotype cells / all selected phenotype
cells) and per trait, among samples with matching genotypes and complete required
covariates. If either rate is **greater than** `multitrait_missing_threshold`
(default **0.30**), the pipeline warns and reverts that ploidy to complete-case
analysis. Exactly 30% is permitted. The fraction of samples with any missing trait
is reported separately and does not trigger this threshold. Samples with every
selected trait missing, or missing required covariates, are excluded.

Trait effects and SEs are returned in their original units. The output folder is
`<projname>_multitrait_gwas/`, including full trait covariance matrices, sample
exclusions, missingness/fallback reports, conditional imputations, convergence
diagnostics and the existing report/plot styles. Imputations use the selected
traits, fixed covariates and background kinship; LOCO exports are chromosome-specific.

For two traits, use e.g. `trait_names = c("NI", "Inoc")`. Do not also include an
exact or near-exact composite of those outcomes. No microbial inputs or selectors
are needed. Existing multi-name `sgwas-marginal` calls still run separate scans.

See [MULTITRAIT_GUIDE.md](MULTITRAIT_GUIDE.md) for the model, approximation,
missing-data policy, multiplicity families and practical limits, and
[run_multitrait_example.R](run_multitrait_example.R) for a commented call.

---

## Select microbiome or phenomic features and construct an outcome

Two arguments control feature construction: `microbial_features` selects input features, and `microbial_method` defines the numeric outcome.

| Goal | `microbial_features` | `microbial_method` | `microbial_response` |
|---|---|---|---|
| Fv abundance | `"Fusarium_verticillioides"` | `"single"` | `NULL` |
| Total Fusarium abundance | `"Fusarium*"` | `"sum"` | `NULL` |
| Fusarium PC1, PC2, or PC3 | `"Fusarium*"` | `"pc1"`, `"pc2"`, or `"pc3"` | `NULL` |
| Phenotype-informed Fusarium proxy | `"Fusarium*"` | `"spls"` | `"FER"` or `"DTA"` |
| Phenotype-informed whole-microbiome proxy | `"all"` | `"spls"` | `"FER"` or `"DTA"` |
| Hyperspectral PC1 | `"HS_*"` | `"pc1"` | `NULL` |
| Apply saved microbiome/phenomic weights | `NULL` | `"weighted"` | `NULL` |

### Independent species abundance or composite (`gwaspoly-marginal` or `sgwas-marginal`)

To map a constructed microbial score as a single trait without a measured host phenotype (`phenofile = NULL`), specify `association_mode = "gwaspoly-marginal"` (or prepare the trait with `run_gwas = FALSE` and map the generated `*_gwas_phenotypes.tsv` with `association_mode = "sgwas-marginal"`):

```r
systemsGWAS(
  wdir = data_dir,
  projname = "Fusarium_PC1_GWAS",
  ploidy_levels = 2,
  genofile_2x = "genotypes.tsv",
  metag_data_species = "metag_species.tsv",
  trait_names = NULL,
  microbial_features = "Fusarium*",
  microbial_method = "pc1",
  microbial_trait_name = "Fusarium_PC1",
  association_mode = "gwaspoly-marginal",
  method = "MLM"
)
```

### Learn a phenotype-informed proxy, then apply it to an independent cohort

First, estimate stable feature weights without running GWAS (`run_gwas = FALSE`):

```r
proxy <- systemsGWAS(
  wdir = data_dir,
  projname = "DTA_proxy_training",
  phenofile = "traits_2020.tsv",
  metag_data_species = "metag_2020.tsv",
  microbial_features = "all",
  microbial_method = "spls",
  microbial_response = "DTA",
  microbial_trait_name = "DTA_microbial_proxy",
  microbial_bootstraps = 100L,
  microbial_stability = 0.60,
  run_gwas = FALSE
)
```

Then apply the exported weights to a validation/test cohort (`microbial_method = "weighted"`). When running integrated modes (`"sgwas-all"`, `"sgwas-joint"`, `"sgwas-interaction"`), weighted proxy outputs are automatically organized inside `<projname>_integrated_gwas/<projname>_microbial_traits/`:

```r
systemsGWAS(
  wdir = data_dir,
  projname = "DTA_2021_integrated",
  ploidy_levels = 2,
  genofile_2x = "genotypes.tsv",
  phenofile = "traits_2021.tsv",
  trait_names = "DTA",
  metag_data_species = "metag_2021.tsv",
  microbial_method = "weighted",
  microbial_weights = "DTA_proxy_training_microbial_traits/DTA_microbial_proxy_weights.tsv",
  microbial_trait_name = "DTA_microbial_proxy",
  association_mode = "sgwas-all",
  method = "MLM",
  model_effect = c("Add", "Dom")
)
```

---

## Joint trait–microbial GWAS and SNP × microbial-score GWAS (`sgwas`)

Set `association_mode` to one of the following (or omit `association_mode` when supplying `microbial_features` / `microbial_weights` to auto-select `"sgwas-all"`):

| Value | Analyses performed |
|---|---|
| `"sgwas-marginal"` (default for phenotype-only) | Base-R REML/GLS single-trait scan (`*_base_gwas/`) |
| `"sgwas-multitrait"` | All selected phenotype columns in one multivariate mixed model, plus matched marginal scans (`*_multitrait_gwas/`) |
| `"gwaspoly-marginal"` | Legacy **GWASpoly** single-trait scan |
| `"sgwas-joint"` | Bivariate measured-trait + microbial-score mixed model, plus matched marginal scans (`*_integrated_gwas/`) |
| `"sgwas-interaction"` | Measured-trait model with SNP, microbial score, and SNP × score, plus matched marginal scans (`*_integrated_gwas/`) |
| `"sgwas-all"` (default when microbial features/weights are passed) | Both bivariate and interaction analyses plus matched marginal baselines (`*_integrated_gwas/`) |

Example: Maize DTA paired with an unsupervised microbial PC1 across additive and dominance models with LOCO kinship:

```r
source("R/systemsGWAS.R")
systemsGWAS(
  wdir = "./",
  projname = "DTA_microbiome_integrated",
  ploidy_levels = 2,
  genofile_2x = "Zm_2x_rd6_maf0.02_dose.txt",
  phenofile = "Traits_DTA_lsmeans.txt",
  trait_names = "DTA",
  metag_data_species = "metag_species.tsv",
  microbial_features = "all",
  microbial_method = "pc1",
  microbial_trait_name = "Microbiome_PC1",
  microbial_transform = "clr",
  covariate_metag = FALSE,
  association_mode = "sgwas-all",
  association_n_pcs = 0L,
  association_max_missing = 0.20,
  association_chunk_size = 256L,
  association_maxit = 300L,
  association_geno_freq = 0.99,
  association_output_style = "verbose",
  association_plot_formats = c("pdf", "tiff"),
  method = "MLM",
  model_effect = c("Add", "Dom"),
  maf = 0.02,
  LOCO = TRUE
)
```

---

## Key arguments and defaults

| Argument | Default | Purpose |
|---|---|---|
| `microbial_features` | `NULL` | Select microbiome/phenomic feature names, a prefix (`"Fusarium*"`), or `"all"` |
| `microbial_method` | `"auto"` | Choose `single`, `sum`, `pc1`, `pc2`, `pc3`, `spls`, or `weighted` |
| `microbial_response` | `NULL` | Training phenotype column for `spls` |
| `microbial_weights` | `NULL` | Complete exported `*_weights.tsv` file for reuse |
| `microbial_trait_name` | `"Microbial_composite"` | Numeric outcome column and output filename prefix |
| `microbial_transform` | `"auto"` | `"auto"` (relative for `sum`, CLR otherwise), `"none"`, `"relative"`, `"clr"`, or `"log1p"` |
| `microbial_prevalence` | `0.05` | Minimum fraction of samples with positive values (`0` for phenomic tables) |
| `microbial_bootstraps` | `100L` | Number of sPLS bootstrap resamples |
| `microbial_stability` | `0.30` | Minimum bootstrap selection frequency (`0.60` recommended for stringent selection) |
| `microbial_seed` | `123L` | Reproducible sPLS resampling seed |
| `microbial_pseudocount` | `1e-6` | Positive offset for CLR after conversion to proportions |
| `run_gwas` | `TRUE` | Set `FALSE` to construct/export a microbiome or phenomic outcome only |
| `association_mode` | `"sgwas-marginal"` (auto-switches to `"sgwas-all"` with microbial features) | `"sgwas-marginal"`, `"gwaspoly-marginal"`, `"sgwas-joint"`, `"sgwas-interaction"`, `"sgwas-all"`, or `"sgwas-multitrait"` |
| `model_effect` | `c("Add", "Dom")` (`"Add"` default in integrated modes if omitted) | Additive (`"Add"`) and/or dosage-threshold dominance (`"Dom"`) genetic models |
| `method` | `c("GLM", "MLM")` (`"MLM"` default in `sgwas` modes) | Statistical model: mixed linear model with GRM (`"MLM"`) or fixed-effect OLS without GRM (`"GLM"`) |
| `LOCO` | `FALSE` | Use leave-one-chromosome-out kinship matrices in `MLM` scans (`TRUE` recommended for `sgwas-env`) |
| `env_deflate_pcs` | `0L` | Number of leading genetic eigenvectors to deflate from the GRM before REML fitting (prevents macro-structure over-correction in `sgwas-env`) |
| `env_gc` | `TRUE` | Compute genomic inflation factor ($\lambda_{GC}$) and export calibrated `P_Trait_GC` and `P_Trait_GC_BH` columns |
| `association_n_pcs` | `0L` | Additional genotype PCs in integrated fixed effects |
| `association_max_missing` | `0.50` | Maximum missing dosage fraction per marker in `sgwas` scans |
| `multitrait_missing` | `"impute"` | Multi-trait only: conditional imputation with observed-data inference, or `"complete_case"` |
| `multitrait_missing_threshold` | `0.30` | Multi-trait phenotype missingness limit; exceeding it overall or in any trait triggers complete-case fallback |
| `association_geno_freq` | `0.99` | Maximum genotype/dominance class frequency filter |
| `association_chunk_size` | `256L` | Markers per vectorized `sgwas` scan block |
| `association_maxit` | `300L` | Optimizer iteration cap per bivariate or multi-trait REML start |
| `association_output_style` | `"verbose"` | `"verbose"` (hypothesis-specific folders + faceted plots) or `"compact"` |
| `association_plot_formats` | `c("tiff")` | Output plot formats (`"pdf"`, `"tiff"`, both, or `NULL` to skip plots) |

---

## Output files

* **Environmental GWAS (`<projname>_env_gwas/`):** Outputs from `association_mode = "sgwas-env"`, including `ASSOCIATIONS.tsv` (with `Lambda_GC`, `P_Trait_GC`, and `P_Trait_GC_BH`), `NULL_MODELS.tsv`, `SAMPLE_AUDIT.tsv`, and ploidy-specific `Trait_marginal/` report folders.
* **Feature preparation (`<projname>_microbial_traits/`, or nested inside `<projname>_integrated_gwas/<projname>_microbial_traits/` for weighted integrated runs):** Exports `*_composite.tsv`, `*_weights.tsv` (complete reference set), `*_selected_weights.tsv`, `*_gwas_phenotypes.tsv`, `*_summary.tsv`, `*_bootstrap.tsv` (for sPLS), and `*_model.rds`.
* **Multi-trait phenotype GWAS (`<projname>_multitrait_gwas/`):** Long-form `ASSOCIATIONS.tsv`, `TRAIT_SET.tsv`, `TRAIT_COVARIANCES.tsv`, `NULL_MODELS.tsv`, `REML_ATTEMPTS.tsv`, `MULTITRAIT_MISSINGNESS.tsv`, sample audits, metadata and numbered trait/omnibus report folders. Each ploidy also has `IMPUTATION_MANIFEST.tsv` and `imputation/` exports (`COMPLETED_TRAITS.tsv`, `IMPUTED_VALUES.tsv`, `IMPUTATION_MASK.tsv`). BH FDR and Bonferroni are calculated across markers, encodings, traits and ploidies separately for each of the three test families; separate calls are not pooled.
* **Base-R univariate GWAS (`<projname>_base_gwas/`):** Outputs from `association_mode = "sgwas-marginal"`.
* **Integrated GWAS (`<projname>_integrated_gwas/`):** Outputs `ASSOCIATIONS.tsv`, `NULL_MODELS.tsv`, `SAMPLE_AUDIT.tsv`, `MODEL_DETAILS.rds`, `REPORT_MANIFEST.tsv`, `RUN_STATUS.tsv`, and ploidy/hypothesis-specific subdirectories (`2x/Trait_marginal/`, `2x/Bivariate_joint/`, `2x/Interaction/`, etc.) with wide/long tables and Manhattan, QTL-profile, and QQ plots.

---

## Development and AI Disclosure

**systemsGWAS** was designed, architected, and validated by **Dr. Bode A. Olukolu** and collaborators. Large language model (LLM) coding assistants were used during software development as productivity aids for code refactoring, drafting unit test assertions, and editing technical documentation. All statistical formulations, matrix algebra implementations, software verification checks, and biological analyses were specified, reviewed, tested, and verified by the authors.

---

## Citation and Contributing

Until the methodological manuscript is published, please cite this repository and record the commit hash, R version, and package versions used. Please also cite underlying methods and packages where applicable, including [GWASpoly](https://jendelman.github.io/GWASpoly/GWASpoly.html), [AGHmatrix](https://cran.r-project.org/package=AGHmatrix), [mixOmics](https://mixomics.org/), and [qvalue](https://bioconductor.org/packages/qvalue/).

Bug reports and feature requests are welcome via [GitHub Issues](https://github.com/bodeolukolu/systemsGWAS/issues).

## Contact

**Dr. Bode A. Olukolu**  
Science Leader – Crop Genetics  
Department of Primary Industries, Queensland, Australia  
Email: [bolukolu@utk.edu](mailto:bolukolu@utk.edu) | [bode.olukolu@dpi.qld.gov.au](mailto:bode.olukolu@dpi.qld.gov.au)

## License

Apache License 2.0. See the [Apache License 2.0 text](https://www.apache.org/licenses/LICENSE-2.0).
