<p align="center">
  <img src="misc/systemsGWAS_logo.png" alt="systemsGWAS logo" width="620">
</p>

# systemsGWAS

**Genome-wide association analysis of measured phenotypes, microbiome and phenomic features, and composite biological traits.**

**systemsGWAS** is an R pipeline that combines marker dosage data with measured phenotypes or traits constructed from microbiome/metagenome and phenomic features, including numeric measurements extracted from hyperspectral or drone imagery. It supports GWAS of individual features, unsupervised composites, phenotype-informed sPLS proxies, joint measured-trait–feature-score association, and SNP × feature-score interactions. Phenomic feature tables can replace microbiome tables; microbiome and phenomic features can also be combined after source-specific preprocessing, with feature centring and scaling performed internally for PCA and sPLS. Standard association analyses use **GWASpoly** and **AGHmatrix** for marker-trait association analysis and computation of genetic relatedness matrix (GRM). Integrated analyses use a separate base-R REML/GLS engine with additive dosage effects and genomic GRM/kinship. Both accept diploid, tetraploid, hexaploid, and octoploid dosage files.

The existing `metag_data_species` and `microbial_*` argument names are used for both microbiome and phenomic inputs. The pipeline reads extracted numeric features, not raw image files. See [phenomic input requirements](#phenomic-features-from-hyperspectral-or-drone-imaging) for transformations and the current nonnegative-input requirement.

Formerly **multiomicGWAS**. The repository is now [bodeolukolu/systemsGWAS](https://github.com/bodeolukolu/systemsGWAS).

## What can I analyse?

| Analysis | GWAS outcome | Does constructing the outcome use a target phenotype such as FER? |
|---|---|---|
| Conventional GWAS | A measured trait, such as Fusarium ear rot (FER) severity | Uses that measured trait directly |
| Single-feature GWAS | One microbial abundance or phenomic measurement, such as a spectral band or canopy-cover estimate | No |
| Summed-abundance GWAS | Total abundance of a selected group of species | No |
| PCA-composite GWAS | PC1, PC2, or PC3 of selected microbiome or phenomic features, mapped separately | No |
| sPLS-proxy GWAS | Weighted microbiome or phenomic component trained against a target phenotype | Yes, when learning the weights |
| Saved-weight GWAS | A numeric score calculated using previously exported weights | No new target observations are required; the original weights retain their training history |
| Joint trait–score GWAS | A measured trait and one microbial/phenomic score as two outcomes | Depends on the score constructor; unsupervised scores are preferred |
| SNP × score GWAS | Measured trait; the SNP effect can depend on the microbial/phenomic score | Depends on the score constructor; both main effects are retained |

Current capabilities also include marker filtering, additive and dominance effect models, numeric phenotype covariates, optional leave-one-chromosome-out (LOCO) kinship, association tables, and Manhattan and QQ plots. A legacy microbiome-adjustment pathway is retained for existing analyses.

**Scope:** in `association_mode = "standard"`, multiple `trait_names` request separate single-trait scans. In `"joint"` or `"both"`, each measured trait is paired with the selected microbial/phenomic score in a bivariate mixed model. General multi-trait sets and plot-level G×E models remain roadmap items.

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

Use compatible R and Bioconductor releases. The microbial-composite module was tested with R 4.5.0 and mixOmics 6.37.0 on Linux. The standard GWASpoly engine loads its broader dependency set even when some packages are not needed. The integrated engine uses base R and does not install or load the standard-engine packages; sPLS score construction still requires mixOmics.

For **preparation only** (`run_gwas = FALSE`), independent composites and saved-weight application use base R; sPLS additionally needs mixOmics. Genotypes and the GWAS package set are not required for preparation.

## Input files

Use tab-separated files with headers. Sample identifiers must match exactly across the inputs used in an analysis. The tables below illustrate layout only. The standard engine requires at least **1,000 markers passing its relationship-matrix filters**. Integrated scans have no artificial 1,000-marker cutoff, but reliable genomic kinship still requires adequate genome-wide marker coverage.

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
Plant_ID	FER	NI	Inoc	DTA
Line01	12.4	8.1	16.7	65
Line02	20.8	14.3	27.3	68
Line03	NA	10.2	NA	63
```

`phenofile` is required for conventional trait GWAS, sPLS training, or requested phenotype covariates. It is unnecessary for independent microbiome/phenomic outcomes or saved-weight application without those covariates.

### Microbial abundance data

For the new microbial workflow, `metag_data_species` accepts either orientation below.

**Samples in rows:** the ID column can be `Plant_ID`, `SampleID`, `Sample_ID`, or `Taxa`; all remaining columns are numeric features.

```text
Plant_ID	Fusarium_verticillioides	Fusarium_proliferatum	Other_species
Line01	120	30	850
Line02	80	60	860
Line03	150	10	840
```

**Species in rows:** identify features with a `species` or `Feature` column and use sample IDs as the abundance columns. Recognised taxonomy columns, such as `genus`, `family`, `order`, `class`, `phylum`, `kingdom`, and `domain`, are ignored.

```text
species	genus	Line01	Line02	Line03
Fusarium_verticillioides	Fusarium	120	80	150
Fusarium_proliferatum	Fusarium	30	60	10
Other_species	Other	850	860	840
```

Abundances must be finite and nonnegative, with unique sample IDs and feature names. Zero means measured absence; missing measurements are not automatically replaced with zero. Pass the complete intended reference microbiome even when selecting a genus: relative abundance and CLR use all input features as their reference.

### Phenomic features from hyperspectral or drone imaging

Supply a feature table through **the same `metag_data_species` argument**, for example `metag_data_species = "phenomic_features.tsv"`. Use one sample per row with `Plant_ID` followed by numeric measurements, or the transposed layout with a `Feature` column. Features may represent spectral reflectances, nonnegative vegetation indices, canopy cover, height, or other measurements extracted from imagery.

```text
Plant_ID	HS_450nm	HS_550nm	HS_650nm	HS_750nm	HS_850nm	UAV_canopy_cover
Line01	0.08	0.12	0.09	0.35	0.48	0.82
Line02	0.10	0.15	0.11	0.31	0.42	0.74
Line03	0.07	0.11	0.08	0.39	0.51	0.89
```

Image calibration, segmentation, feature extraction, and preparation of a sample-level table happen **outside systemsGWAS**. Do not supply JPEGs, TIFFs, orthomosaics, or hyperspectral image cubes directly. Match the feature-table IDs to the genotypes; resolve repeated plots, dates, or environments appropriately before producing the single row per ID required by this feature reader. It does not fit a plot-design or G×E model during feature construction.

For phenomic tables, explicitly set:

```r
# Include these arguments in the systemsGWAS() call:
# microbial_transform = "none",  # preserve the supplied measurement scale
# microbial_prevalence = 0       # disable the microbial presence/absence filter
```

`"auto"` would otherwise apply CLR to PCA/sPLS inputs. Closing heterogeneous image-derived measurements to a common total is not part of this phenomic workflow. With `"none"`, PCA and sPLS still centre and scale each candidate feature internally; single-feature mode preserves the supplied values. Use `sum` only when summing features has a meaningful interpretation and compatible units.

**Current input constraint:** this feature reader requires all values to be finite and nonnegative, even with `microbial_transform = "none"`. Resolve missing measurements before import. Some vegetation indices, centred/z-scored features, and signed residuals can be negative and therefore cannot be passed directly through this reader. Do not replace negative values with zero merely to satisfy validation. Such individual traits can be analysed through `phenofile`; for a composite of signed features, construct the score externally and import that numeric outcome through `phenofile`. The reader's nonnegative constraint does not apply to ordinary phenotype columns or to the PCA/sPLS scores it generates.

Use consistent feature names, units, extraction settings, and preprocessing when reapplying weights. Learn any data-dependent external preprocessing on the training set when evaluating prediction, and apply those same parameters to validation/new samples. The weights file stores the pipeline's transformation, centring, and scaling; it does not store external image-processing or imputation steps.

## Quick start: conventional trait GWAS

Run the following from the repository root after editing `data_dir` and the filenames:

```r
source("R/systemsGWAS.R")
data_dir <- "/path/to/data"

systemsGWAS(
  wdir = data_dir,
  projname = "FER_GWAS",
  ploidy_levels = 2,
  genofile_2x = "genotypes.tsv",
  phenofile = "traits.tsv",
  trait_names = "FER",
  covariate_pheno = "DTA",
  method = "MLM",
  model_effect = "Add",
  maf = 0.02,
  cores = 1
)
```

Set `covariate_pheno = NULL` if no phenotype covariates are required. For separate scans of several traits, use a vector such as `trait_names = c("NI", "Inoc", "FER")`.

Alternatively, edit [R/run_parameters_systensGWAS.R](R/run_parameters_systensGWAS.R), which describes all public arguments and their defaults. Set its `wdir` to your data directory, then run it from the `R` directory:

```sh
cd R
Rscript run_parameters_systensGWAS.R
```

The filename above retains the existing `systensGWAS` spelling for compatibility. Its `source("systemsGWAS.R")` line assumes the two scripts are in the same working directory.

## Select microbiome or phenomic features and construct an outcome

Two arguments control the main choice: `microbial_features` selects the input features, and `microbial_method` defines the numeric outcome. Their names are retained for compatibility when the table contains phenomic measurements.

| Goal | `microbial_features` | `microbial_method` | `microbial_response` |
|---|---|---|---|
| Fv abundance | `"Fusarium_verticillioides"` | `"single"` | `NULL` |
| Total Fusarium abundance | `"Fusarium*"` | `"sum"` | `NULL` |
| Fusarium PC1, PC2, or PC3 | `"Fusarium*"` | `"pc1"`, `"pc2"`, or `"pc3"` | `NULL` |
| FER-informed Fusarium proxy | `"Fusarium*"` | `"spls"` | `"FER"` |
| FER-informed whole-microbiome proxy | `"all"` | `"spls"` | `"FER"` |
| Hyperspectral PC1 | `"HS_*"` | `"pc1"` | `NULL` |
| FER-informed phenomic proxy | `"all"` or `"HS_*"` | `"spls"` | `"FER"` |
| Apply saved microbiome/phenomic weights | `NULL` | `"weighted"` | `NULL` |

Exact vectors, for example `c("Fusarium_verticillioides", "Fusarium_proliferatum")`, and comma-separated lists are accepted. `"Fusarium*"` matches names starting with `Fusarium`. The compatible selectors `"Fusarium_spp"` and `"Fusarium spp"` match names starting with `Fusarium_` or `Fusarium `. Matching is case-sensitive. Phenomic prefixes work identically: `"HS_*"` selects names beginning with `HS_`, while `c("HS_*", "UAV_*")` selects both groups. Use `"all"` for every feature in the supplied table. For phenomic analyses, also specify `microbial_transform = "none"` and `microbial_prevalence = 0`.

The default `microbial_method = "auto"` chooses `single` for one matched feature and `pc1` for multiple matched features. It does not train on FER. One new feature-derived outcome is constructed per call.

### Choosing PC1, PC2, or PC3

Set `microbial_method` to `"pc1"`, `"pc2"`, or `"pc3"` (case-insensitive) to map one PCA component. Components are ordered by variance explained in the transformed, standardised feature matrix. All three are unsupervised and do not use FER to choose the axis.

Each call constructs one numeric phenotype. To map another component, use another call with a distinct `projname` and `microbial_trait_name`, for example `"Fusarium_PC2_GWAS"` and `"Fusarium_PC2"`. Keep the samples, candidate features, and transformation the same when comparing components from the same PCA. The default `"auto"` continues to choose PC1 when multiple features match.

The requested component must have nonzero numerical variance. Its availability depends on sample size and the effective rank after transformation and feature filtering; having three columns alone does not guarantee PC3. The run stops with the available rank if the component cannot be estimated.

The summary and full weights file report `PC_Component` and `PC_Variance_Explained` (a fraction from 0 to 1). Saved weights reproduce the selected component. The legacy summary column `PC1_Variance_Explained` is retained for compatibility and is `NA` for PC2/PC3.

### Independent species abundance or composite

This example maps Fusarium PC1 without using FER values:

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
  method = "MLM"
)
```

For Fv alone, change `microbial_features` to `"Fusarium_verticillioides"`, `microbial_method` to `"single"`, and the project/trait names accordingly. For total Fusarium abundance, use `microbial_method = "sum"`. Use a distinct project prefix for each analysis.

### Learn an FER-informed proxy, then apply it to GWAS

First, select stable species, estimate weights, and save the numeric proxy without running GWAS:

```r
proxy <- systemsGWAS(
  wdir = data_dir,
  projname = "FER_proxy_training",
  phenofile = "traits.tsv",
  metag_data_species = "metag_species.tsv",
  microbial_features = "Fusarium*",
  microbial_method = "spls",
  microbial_response = "FER",
  microbial_trait_name = "Fusarium_FER_proxy",
  microbial_bootstraps = 100L,
  microbial_stability = 0.60,
  run_gwas = FALSE
)
```

Use `microbial_features = "all"` to search the full microbiome, or provide a list/prefix to restrict the candidate species. Training uses samples with finite observed FER values, but the fitted component scores all microbiome samples, including those without FER observations.

Then apply the exported weights and map the composite:

```r
systemsGWAS(
  wdir = data_dir,
  projname = "FER_proxy_GWAS",
  ploidy_levels = 2,
  genofile_2x = "genotypes.tsv",
  metag_data_species = "metag_species.tsv",
  microbial_method = "weighted",
  microbial_weights = "FER_proxy_training_microbial_traits/Fusarium_FER_proxy_weights.tsv",
  microbial_trait_name = "Fusarium_FER_proxy",
  method = "MLM"
)
```

The second call does not require FER. The weights can also be applied to another cohort with the same reference features and compatible measurement/preprocessing. To learn a proxy and run GWAS in one call, instead set `run_gwas = TRUE` in the training call and supply the required genotype file and ploidy.

### Phenomic PCA components and phenotype-informed proxies

Map a hyperspectral composite without using FER values:

```r
systemsGWAS(
  wdir = data_dir,
  projname = "Hyperspectral_PC1_GWAS",
  ploidy_levels = 2,
  genofile_2x = "genotypes.tsv",
  metag_data_species = "phenomic_features.tsv",
  trait_names = NULL,
  microbial_features = "HS_*",
  microbial_method = "pc1",
  microbial_transform = "none",
  microbial_prevalence = 0,
  microbial_trait_name = "Hyperspectral_PC1",
  method = "MLM"
)
```

For a single measurement, use `microbial_features = "UAV_canopy_cover"` and `microbial_method = "single"`, keeping `microbial_transform = "none"` and `microbial_prevalence = 0`.

To learn a phenomic proxy for FER instead, use sPLS with the observed disease ratings:

```r
phenomic_proxy <- systemsGWAS(
  wdir = data_dir,
  projname = "Phenomic_FER_training",
  phenofile = "traits.tsv",
  metag_data_species = "phenomic_features.tsv",
  microbial_features = "all",
  microbial_method = "spls",
  microbial_response = "FER",
  microbial_transform = "none",
  microbial_prevalence = 0,
  microbial_trait_name = "Phenomic_FER_proxy",
  microbial_bootstraps = 100L,
  microbial_stability = 0.60,
  run_gwas = FALSE
)
```

This exports selected feature weights and the numeric proxy just as for microbiome inputs. Multi-feature sPLS has the same minimum of five variable candidates after filtering; it also requires at least five observed training samples and a variable response. The three-row input illustration above is not a training dataset.

Reuse `Phenomic_FER_training_microbial_traits/Phenomic_FER_proxy_weights.tsv` with `microbial_method = "weighted"` and a compatible phenomic table, following the saved-weight GWAS example above. The stored transformation is `"none"`; it is restored automatically along with the learned centres/scales. A newly supplied transformation does not override the saved model's preprocessing. FER observations are not needed when applying the saved weights unless separately requested as a covariate.

### Combining microbiome and phenomic features

Compatible nonnegative feature blocks can be joined by sample ID into one table and supplied through `metag_data_species`. Use distinct column prefixes, for example `MB_`, `HS_`, and `UAV_`, to retain each feature's source and support prefix selection. Exclude IDs, dates, treatment labels, and other metadata from the numeric feature columns.

Apply appropriate **source-specific preprocessing** before joining the blocks, then use `microbial_transform = "none"` and `microbial_prevalence = 0`. Apply any desired microbiome prevalence filtering separately before constructing the combined table. PCA/sPLS perform **feature standardisation** internally (centring and unit-variance scaling). This does not impose equal total weight on each data source: the pipeline fits the concatenated features and has no separate block-weight control. Retain the same feature definitions and external preprocessing when applying saved weights.

Do not apply CLR to the entire combined table. If the intended preprocessing produces signed values—for example a CLR microbiome block or pre-standardised phenomic block—the current reader cannot accept that combined matrix. Construct the composite externally and map its numeric score through `phenofile` instead.

## Transformations and weights

| `microbial_transform` | Behaviour |
|---|---|
| `"auto"` | Relative abundance for `sum`; CLR for `single`, `pc1`, `pc2`, `pc3`, and `spls` |
| `"none"` | Retain supplied values; use explicitly for compatible phenomic or combined feature tables |
| `"relative"` | Divide abundances by each sample's total |
| `"clr"` | Close to proportions, add the pseudocount, then centre log abundances within each sample |
| `"log1p"` | Apply log(1 + abundance) |

The CLR pseudocount defaults to `1e-6`. CLR measures abundance relative to the reference community; it is not an estimate of absolute microbial load. Do not supply a pre-transformed CLR matrix to the abundance reader, which requires nonnegative values. With only one species in the entire input matrix, use `none` or `log1p`: CLR cannot produce a variable single-column score.

The `"auto"` defaults remain microbiome-oriented. For phenomic features, use `"none"` and disable prevalence filtering as shown above.

PCA and sPLS centre and scale candidate features to unit variance. Each selected PCA component is oriented so that its largest absolute loading is positive; sPLS is oriented towards a positive Spearman correlation with its training response. A negative weight therefore describes an opposing contribution to the composite, not evidence of a protective or causal biological effect.

The exported score is:

```text
score_i = sum_j Weight_j * (transformed_feature_ij - Center_j) / Scale_j
```

sPLS weights are component loadings, not regression coefficients in FER units. A single-feature sPLS component only standardises and orients that feature. Use `sum` when the intended outcome is total abundance rather than a weighted pattern of community variation.

### sPLS feature selection

The default procedure uses 100 row bootstraps, a minimum selection frequency of 0.60, and one final component fitted to the stable subset. Each bootstrap draws `keepX` from `c(5, 10, 20, 30, 50, 100, 200)`, restricted to the number of available candidates; it is not a cross-validated tuning search.

Multi-feature selection, for microbiome or phenomic inputs, requires **at least five variable candidates after filtering**. Five to nine candidates use `keepX = 5`; with exactly five, the screen cannot select a smaller subset. A single-feature score is supported separately. The 5% prevalence default filters candidate features; for sPLS, prevalence and scaling are estimated using the observed training samples.

Selection frequency uses all attempted bootstraps as the denominator. Failed fits do not inflate stability. If no stable features remain, the run stops with a diagnostic message; the bootstrap log is saved.

An sPLS proxy trained on FER is an **FER-informed microbiome or phenomic phenotype**. Its training correlation is descriptive, not cross-validated predictive ability or independent confirmation of FER-associated loci. Evaluate proxy performance using held-out samples or independent data. The GWAS treats the constructed score as fixed and does not propagate uncertainty from estimating the weights.

## Joint trait–microbial GWAS and SNP × microbial-score GWAS

Set `association_mode` to one of the following. The same controls also accept compatible phenomic scores.

| Value | Analyses performed |
|---|---|
| `"standard"` (default) | Existing GWASpoly workflows |
| `"joint"` | Bivariate measured-trait + microbial-score mixed model, plus separate marginal scans on the same samples |
| `"interaction"` | Measured-trait model with SNP, microbial score, and SNP × score, plus marginal scans |
| `"both"` | Both integrated analyses and their matched marginal baselines |

Example: FER paired with an unsupervised Fusarium PC1, with both joint and interaction tests:

```r
source("systemsGWAS.R")
systemsGWAS(
  wdir = "./",
  projname = "FER_Fusarium_integrated",
  ploidy_levels = 2,
  genofile_2x = "Zm_2x_rd6_maf0.02_dose.txt",
  phenofile = "Traits_FER_lsmeans.txt",
  trait_names = "FER",
  metag_data_species = "FER_metagenome_lsmeans.txt",
  microbial_features = "Fusarium_spp",
  microbial_method = "pc1",
  microbial_trait_name = "Fusarium_PC1",
  microbial_transform = "clr",
  microbial_response = NULL,
  trait_microbial_proxy = NULL,
  covariate_metag = FALSE,
  covariate_pheno = "DTA",
  association_mode = "both",
  association_n_pcs = 0L,
  association_max_missing = 0.20,
  association_chunk_size = 256L,
  association_maxit = 300L,
  method = "MLM",
  model_effect = "Add",
  maf = 0.02,
  LOCO = TRUE,
  perm = 1,
  cores = 1
)
```

Use phenotype covariates such as DTA only when appropriate to the scientific question. To use Fv abundance, change `microbial_features` to `"Fusarium_verticillioides"`, `microbial_method` to `"single"`, and `microbial_trait_name` to `"Fv"`. The transform remains explicit: `"none"` uses input abundance, `"relative"` uses proportions, and `"clr"` uses log ratios against the full input reference. PC2 and PC3 are also supported, one component per call. For phenomic inputs, use `microbial_transform = "none"` and `microbial_prevalence = 0`; the existing finite, nonnegative input requirement still applies.

`trait_names = c("NI", "Inoc", "FER")` fits three separate measured-trait–score pairs. It does not fit a four-trait model. `run_gwas = FALSE` still prepares the microbial score only. The measured phenotype is retained as an outcome in integrated modes; it is not replaced by the score or residualised against the microbiome beforehand.

### Models and test interpretation

The bivariate model estimates an unstructured 2 × 2 genetic covariance and an unstructured 2 × 2 residual covariance. Its observation covariance is `Sigma_G %x% K + Sigma_E %x% I`, where `K` is genomic kinship. Each trait has its own fixed-effect coefficients and marker effect. A significant bivariate test may be driven entirely by the microbial trait; inspect the measured-trait effects and p-values separately.

The interaction model is `trait ~ covariates + SNP + microbial_score + SNP:microbial_score + genomic_random_effect`. The microbial score is centred and divided by its analysed-cohort SD, and dosage is centred. `Beta_SNP_at_mean` is the per-dosage-unit SNP effect at mean microbial score. `Beta_SNPxM` is the change in that SNP effect per one SD of the score. Effect signs assume dosage counts the ALT allele; microbiome-mediated causation is not established by either model.

| P-value column | Null hypothesis |
|---|---|
| `P_Trait` | No SNP effect on the measured trait in the marginal model |
| `P_Microbial` | No SNP effect on the microbial score in the marginal model |
| `P_Bivariate` | Both bivariate SNP effects are zero (2-df test) |
| `P_Trait_joint` | Measured-trait SNP effect is zero in the bivariate model |
| `P_Microbial_joint` | Microbial-score SNP effect is zero in the bivariate model |
| `P_SNP_at_mean` | SNP main effect at mean microbial score is zero in the interaction model |
| `P_SNPxM` | Interaction coefficient is zero, retaining both main effects |
| `P_SNP_joint` | SNP main effect and interaction coefficient are both zero (2-df test) |

### Fitting, speed, and safeguards

- **Covariance fitting:** univariate REML is profiled over the genetic variance fraction. Bivariate REML uses three deterministic optimizer starts. Failed bivariate fits stop the scan; the pipeline never substitutes an uncorrected model. Boundary fits are labelled, and genetic correlations are suppressed for flagged boundary fits.
- **Marker tests:** covariance shape is estimated under the null and reused for GLS scans (P3D). A scalar residual variance is re-estimated in each full marker model. Reported F tests are approximate because covariance was estimated; these are not per-marker covariance-refitted likelihood-ratio tests or finite-sample exact tests.
- **Efficiency:** genotype files are read once per ploidy; spectral decompositions and null fits are reused within each trait pair/chromosome. `association_chunk_size` limits temporary scan matrices. Integrated scans are serial and vectorised; `cores` remains a GWASpoly control. The full genotype matrix is still held in memory.
- **Matched cohorts:** IDs are matched by name. All analyses within a pair use the same samples and marker filters. Missing measured outcomes or numeric covariates cause complete-case exclusion, not phenotype imputation. Duplicate phenotype IDs are averaged, as in the existing preparation workflow; these are genotype-level analyses, not plot-design models.
- **Marker QC:** MAF and missingness are calculated on the analysed cohort. Monomorphic/all-missing markers are excluded. Passing missing dosages are imputed with the cohort marker mean. Kinship is the centred-dosage cross-product normalised to mean diagonal one. No minimum of 1,000 markers is imposed in this branch.
- **LOCO:** `LOCO = TRUE` excludes the tested chromosome from kinship, with fresh null fits per chromosome. At least two retained chromosomes are needed. Optional genotype PCs come from whole-genome kinship and remain the same across chromosomes.
- **Rank checks:** constant outcomes/covariates, dependent outcome pairs, unidentifiable bivariate covariance, and collinear covariates fail clearly. Marker-specific rank or residual failures receive `NA` tests and a status column.
- **Current scope:** additive dosage effects only, `method = "MLM"`. If `method` and `model_effect` are omitted, integrated modes choose MLM/Add; explicitly supplied unsupported settings fail. `perm >= 100` is rejected because permutation calibration is not implemented for these modes.
- **Supervised proxies:** sPLS-based scores, including saved sPLS weights, trigger an exploratory-inference label and warning. Standard p-values do not propagate weight-learning uncertainty. Independent training or an appropriate full-procedure calibration is needed for confirmatory use. Unknown saved-weight provenance is explicitly labelled.

Raw p-values have `_BH` and `_Bonferroni` counterparts, adjusted across markers, measured-trait pairs, and ploidies in the current call **separately for each test column**. This does not correct across different test families or separate calls/components. Plan those additional comparisons explicitly. Improved p-values or additional hits alone are not proof of improved discovery power.

### Integrated outputs

Outputs are written under `<projname>_integrated_gwas/`:

| File | Contents |
|---|---|
| `ASSOCIATIONS.tsv` | Combined marker effects, standard errors, raw/adjusted p-values, sample counts and status |
| `NULL_MODELS.tsv` | Convergence, boundary flags, variance/covariance estimates and correlations |
| `SAMPLE_AUDIT.tsv` | Candidate microbial sample IDs, availability checks and cohort inclusion |
| `MODEL_DETAILS.rds` | Model fits, covariance information, construction provenance, settings and R session |
| `RUN_STATUS.tsv` | `running`, `failed`, or `complete`; check before using outputs from a reused prefix |
| `pair*/associations.tsv` | Results for one trait/score/ploidy pair |
| `pair*/marker_qc.tsv` | Marker MAF, missingness, variation and retention |
| `pair*/analysed_samples.tsv` | Matched responses, standardised score and fixed-effect design |
| `pair*/Manhattan_QQ.pdf` | Manhattan and QQ plots for each available test |
| `pair*/FAILED.tsv` | Fit error if a chromosome/pair fails; removed when that pair later succeeds |

The existing `<projname>_microbial_traits/` files still record the microbial/phenomic score construction. Use distinct project names for distinct analyses; the new modes do not implement checkpoint/resume.

The model formulation follows the multi-trait mixed-model approach described by [Korte et al. (2012)](https://doi.org/10.1038/ng.2376). [Zhou and Stephens (2014)](https://doi.org/10.1038/nmeth.2848) discuss efficient multivariate mixed models and the limitations of covariance approximations. This implementation uses base-R REML/GLS, not GEMMA or ASReml.

## Key arguments and defaults

| Argument | Default | Purpose |
|---|---|---|
| `microbial_features` | `NULL` | Select microbiome/phenomic feature names, a prefix, or `"all"`; activate construction |
| `microbial_method` | `"auto"` | Choose `single`, `sum`, `pc1`, `pc2`, `pc3`, `spls`, or `weighted` |
| `microbial_response` | `NULL` | Training phenotype column for sPLS |
| `microbial_weights` | `NULL` | Complete exported weights file for reuse |
| `microbial_trait_name` | `"Microbial_composite"` | Numeric outcome column and output filename prefix |
| `microbial_transform` | `"auto"` | Transformation; explicitly use `"none"` for phenomic features; saved weights restore their original setting |
| `microbial_prevalence` | `0.05` | Minimum fraction of samples with positive values; set `0` for phenomic tables |
| `microbial_bootstraps` | `100L` | Number of sPLS resamples |
| `microbial_stability` | `0.60` | Minimum bootstrap selection frequency |
| `microbial_seed` | `123L` | Reproducible sPLS resampling seed |
| `microbial_pseudocount` | `1e-6` | Positive offset for CLR after conversion to proportions |
| `run_gwas` | `TRUE` | Set `FALSE` to construct/export a microbiome or phenomic outcome only |
| `association_mode` | `"standard"` | `standard`, `joint`, `interaction`, or `both` |
| `association_n_pcs` | `0L` | Additional genotype PCs in integrated fixed effects |
| `association_max_missing` | `0.20` | Maximum missing dosage fraction per marker in integrated scans |
| `association_chunk_size` | `256L` | Markers per vectorised integrated scan block |
| `association_maxit` | `300L` | Optimizer iteration cap per bivariate REML start |

Use letters, digits, and underscores for `microbial_trait_name`, beginning with a letter. `projname` is a simple directory prefix. File paths can be absolute or relative to `wdir`.

For both engines, `maf = 0.02` retains markers with MAF at least 2%; `LOCO = FALSE` uses whole-genome kinship. `cores = 1` controls GWASpoly workers only. Full argument descriptions are in the parameter script.

## Output files

Feature preparation writes to `<projname>_microbial_traits/`. This directory name is also used for phenomic inputs for compatibility. Each filename begins with `microbial_trait_name`.

| Suffix | Contents |
|---|---|
| `_composite.tsv` | Sample ID and numeric constructed outcome |
| `_weights.tsv` | Complete reference feature set, weights, centres, scales, and transformation metadata |
| `_selected_weights.tsv` | Nonzero-weight microbial or phenomic features for interpretation |
| `_gwas_phenotypes.tsv` | Outcome plus explicitly requested numeric phenotype covariates |
| `_summary.tsv` | Method, sample/feature counts, bootstrap success, and selected-component variance or sPLS training correlation |
| `_bootstrap.tsv` | sPLS resampling settings, success, and error details |
| `_model.rds` | Construction metadata, training IDs, and weights |

**Reuse the full `_weights.tsv` file for both microbiome and phenomic models.** The loader requires the complete saved feature set, including zero-weight rows, even with `"none"`. Zero-weight reference taxa are also needed to reproduce CLR and relative abundance. The selected-only table is for interpretation. Saved centres/scales are reused; reordered columns are matched by name, extra columns are ignored, and missing reference taxa cause an error.

GWAS outputs are organised into project/model/ploidy directories containing association scores and effects, threshold-specific result tables, phenotype summaries, and Manhattan/QQ plots. Empty output directories are removed. Reusing project and trait names can overwrite outputs; use distinct prefixes for distinct analyses.

## Compatibility and implementation notes

- Use the existing `metag_data_species` and `microbial_*` controls for both microbiome and compatible phenomic feature tables. No `phenomic_*` argument aliases are currently implemented. Leave `alternate_trait`, `trait_microbial_proxy`, `pheno_taxa_species`, and `pheno_taxa_strain` as `NULL` when using them. Legacy names remain available for existing workflows.
- Keep `covariate_metag = FALSE` for the new feature-derived outcome workflow, including phenomic inputs. The legacy adjustment mode changes the response by fitting a microbiome model; it is not merely an input-file switch.
- Examples use `method = "MLM"`. The legacy `"GLM"` branch modifies the supplied kinship and should not be interpreted as a conventional no-kinship GLM. LOCO rebuilds kinship; use it with MLM.
- In standard mode, `model_effect` selects diploid additive/dominance families; higher-ploidy branches use predefined model sets. Integrated modes currently support additive dosage effects only.
- Numeric phenotype covariates can be given as a vector or comma-separated string. In standard mode, missing values use mice imputation and the genotype-PC count follows covariate count. Integrated modes use complete cases and an explicit `association_n_pcs` setting.
- `fdr` and `bonferroni` control threshold lines in plots; the corresponding tables are still calculated when the switches are `FALSE`. `suggestive` is a −log10(p) cutoff applied after an initial FDR 0.5 candidate filter. `perm >= 100` activates permutation thresholds.
- `cores` controls GWASpoly parallelism, not the sPLS bootstrap loop.
- If R reports unused `microbial_*` arguments, source the updated `R/systemsGWAS.R` rather than an earlier copy.

### Validation status

The composite module, including PC2/PC3 support, passed 42 targeted synthetic checks using real mixOmics fits. The integrated module passed 38 checks, including agreement with independent dense-matrix REML/GLS calculations, additive and interaction signal recovery, LOCO, phenotype-scale invariance, complete-case alignment, multiple trait/ploidy handling, exported tables/PDFs, and failure guards. End-to-end integrated scans used 96 synthetic samples and 84 input markers; these exercised the real new fitting engine, not a test double.

An additional null simulation used 80 independent phenotype replicates on a fixed 80-sample genotype matrix, with 40 tested markers per replicate. At nominal alpha 0.05, observed rejection rates were approximately 0.0513 (trait marginal), 0.0466 (bivariate), 0.0469 (interaction), and 0.0456 (joint SNP). This is a limited calibration smoke test, not proof of calibration for every dataset or at genome-wide tail thresholds. The packaged scripts reproduce the numerical checks and calibration.

Validation was on Linux with R 4.5.0 and mixOmics 6.37.0. A complete legacy GWASpoly scan, Windows execution, and validation on the biological FER dataset were not performed. Inspect convergence, cohort/QC tables, and QQ plots before interpreting study results.

## Roadmap

The following are development directions, not current native model capabilities or committed release dates:

- General joint multi-trait sets beyond the implemented measured-trait–feature-score pairs.
- Multi-environment, stability, plasticity, and G×E association models.
- Dedicated constructors for additional latent phenotypes, including the proposed Developmental Blueprint Index (DBI).
- Native signed-feature input, explicit multi-block weighting, and broader integration of transcriptomic, metabolomic, phenomic, environmental, and management data beyond the compatible numeric-table workflow described above.
- User-defined kernels and links to crop systems models such as APSIM.

Externally prepared numeric traits can already be supplied through `phenofile`; this is distinct from implementing a dedicated model or constructor for each data type.

## Citation and contributing

A methodological publication will be added when available. Until then, cite this repository and record the commit used, together with the R and package versions for your analysis. Also cite the relevant underlying methods and packages, including [GWASpoly](https://jendelman.github.io/GWASpoly/GWASpoly.html), [AGHmatrix](https://cran.r-project.org/package=AGHmatrix), and [mixOmics](https://mixomics.org/) where used.

Bug reports, feature requests, and pull requests are welcome through [GitHub Issues](https://github.com/bodeolukolu/systemsGWAS/issues). For a reproducible report, include the call, error/traceback, `sessionInfo()`, and a small example that preserves the relevant input structure.

## Contact

**Dr. Bode A. Olukolu**  
Science Leader – Crop Genetics  
Department of Primary Industries, Queensland, Australia  
Email: [bolukolu@utk.edu](mailto:bolukolu@utk.edu)

## License

Apache License 2.0. See the [Apache License 2.0 text](https://www.apache.org/licenses/LICENSE-2.0).
