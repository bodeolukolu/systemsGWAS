#!/usr/bin/env Rscript
# Source the updated pipeline, then edit the input filenames below.
source("systemsGWAS.R")

systemsGWAS(
  # Data directory and output prefix. Defaults: "./", "GWAS".
  wdir = "./",
  projname = "Fusarium_PC1",
  
  # Dosage scales to scan; each needs its matching genotype file. Default: c(2,4,6,8).
  ploidy_levels = 2,
  # Marker-row TSV: SNP, CHROM, POS, REF, ALT, then sample dosages 0..ploidy.
  # Each file defaults to NULL; only selected ploidies need a file.
  genofile_2x = "Zm_2x_rd6_maf0.02_dose.txt",
  genofile_4x = NULL,
  genofile_6x = NULL,
  genofile_8x = NULL,
  
  # Numeric phenotype TSV: sample IDs first, then traits/covariates. Default: NULL.
  # Needed for integrated GWAS, sPLS or covariates; independent microbial GWAS can omit it.
  # For spls, supply the measured FER file here, not a generated composite filename.
  # Standard proxy GWAS automatically reads <projname>_microbial_traits/*_gwas_phenotypes.tsv.
  phenofile = NULL,
  # Measured outcomes to map; required for joint/interaction/both.
  # In standard microbial mode the constructed score replaces these outcomes.
  # Default: c("trait1","trait2"). Use NULL in the microbial examples below.
  trait_names = NULL,
  
  # Species TSV: species/Feature rows with sample columns, OR Plant_ID/SampleID
  # rows with numeric species columns. Taxonomy metadata is allowed in row format.
  # Abundances must be nonnegative; IDs/names must be unique. Default: NULL.
  metag_data_species = "metag_species.tsv",
  
  # NEW MICROBIAL PHENOTYPE CONTROLS
  # Exact species, vector of names, "Fusarium*", "Fusarium_spp", or "all".
  # Prefix matching is case-sensitive and starts at the beginning. Default: NULL.
  # NULL keeps ordinary/legacy mode unless microbial_weights is supplied.
  microbial_features = "Fusarium*",
  # "single": one feature; "sum": total abundance; "pc1"/"pc2"/"pc3": one PCA axis;
  # Each PCA run maps one component; use distinct project/trait names for separate axes.
  # "spls": phenotype-trained, bootstrap-selected component; "weighted": reuse weights.
  # Default "auto" chooses single for one match, otherwise pc1; never trains on FER.
  microbial_method = "pc1",
  # Training phenotype column for spls only, e.g. "FER". Default: NULL.
  microbial_response = NULL,
  # Full *_weights.tsv from an earlier run; method must be "weighted". Default: NULL.
  # Reuses transformation, reference taxa, centring and scaling without retraining.
  microbial_weights = NULL,
  # Output column/file prefix; letters, digits, underscores. Default: "Microbial_composite".
  microbial_trait_name = "Fusarium_PC1",
  # "none", "relative", "clr" or "log1p". Default "auto": relative for sum, clr otherwise.
  # Relative/CLR use the whole input microbiome as reference, not just Fusarium.
  # For an already prepared abundance/transform use "none" (nonnegative inputs only).
  microbial_transform = "auto",
  # Minimum fraction of samples with abundance >0. Default: 0.05.
  # sPLS uses observed training samples; weighted mode reuses the saved reference set.
  microbial_prevalence = 0.05,
  # Number of sPLS bootstraps. Default: 100L. Ignored by other methods.
  microbial_bootstraps = 100L,
  # Minimum selection frequency across attempted bootstraps. Default: 0.30.
  # 0.60 is more conservative; neither threshold controls false discoveries.
  # Changing this needs a fresh spls run; weighted mode retains its saved feature set.
  # Multi-species sPLS needs >=5 variable candidates; 5-9 candidates use keepX=5.
  microbial_stability = 0.30,
  # sPLS bootstrap seed. Default: 123L; the caller's random-number state is restored.
  microbial_seed = 123L,
  # Positive offset added to relative abundance before CLR. Default: 1e-6.
  microbial_pseudocount = 1e-6,
  # TRUE builds the trait then runs GWAS; FALSE saves preparation only. Default: TRUE.
  # Preparation needs no genotype file and saves weights, scores and diagnostics.
  run_gwas = TRUE,
  
  # INTEGRATED TRAIT + MICROBIAL/PHENOMIC ASSOCIATION
  # "gwaspoly-marginal": existing scans. Runs single trait without microbiome using GWASpoly; 
  # "sgwas-marginal (default)": Runs single trait without microbiome using Base-R REML/GLS;
  # "sgwas-multitrait": two or more traits
  # "sgwas-joint": two-trait mixed model;
  # "sgwas-interaction": SNP x microbial score; 
  # "sgwas-all": run both, with matched baselines.
  # New modes require phenofile + trait_names and one constructed microbial score.
  # Use method="MLM", model_effect="Add", "Dom", or both; covariate_metag=FALSE.
  association_mode = "sgwas-marginal",
  # Extra genomic PCs as fixed covariates; kinship is always fitted. Default: 0L.
  association_n_pcs = 0L,
  # Maximum missing dosage fraction per marker in the analysed cohort. Default: 0.20.
  # Missing calls use each encoding's mean; missing outcomes/covariates are excluded.
  association_max_missing = 0.50,
  # Markers processed together; lower values reduce temporary memory. Default: 256L.
  association_chunk_size = 256L,
  # Optimizer iteration cap per bivariate REML start; failed fits stop. Default: 300L.
  association_maxit = 300L,
  association_geno_freq=0.99,
  # Integrated reporting: "gwaspoly" (default) uses familiar folders and faceted plots.
  # "compact" keeps the earlier combined plots; this option does not change model fitting.
  association_output_style = "verbose",
  # Integrated plot formats; PDF and 300-dpi TIFF by default. NULL skips plots, not tables.
  association_plot_formats = c("pdf", "tiff"),
  
  ## "impute" or "complete_case".
  multitrait_missing = "impute",
  # Missingness ≤30%: retain partially observed samples and estimate missing phenotypes using the observed traits, covariates and kinship.
  # Missingness >30% overall or in any selected trait: issue a warning and revert to dropping samples missing any selected trait. Exactly 30% is permitted.
  # >30% overall OR in any trait: drop incomplete rows.
  multitrait_missing_threshold = 0.30,
  
  # GWAS CONTROLS
  # "MLM", "GLM", or both. Default: c("GLM","MLM"). Prefer MLM for kinship correction.
  # Legacy GLM modifies the supplied kinship; it is not a conventional no-kinship fit.
  method = "MLM",
  # Standard diploid effect families: "Add" or c("Add","Dom"). Default: c("Add","Dom").
  # Integrated: "Add", "Dom", or both; omitted arguments choose MLM/Add.
  # Dom scans d-dom-alt/ref separately, d=1..ploidy/2 (not simultaneous Add + Dom).
  # Fractional dosages are rounded for Dom, as in GWASpoly; Add retains dosage values.
  # Legacy higher-ploidy branches run their predefined additive/dominance sets.
  model_effect = c("Add","Dom"),
  # Numeric phenotype covariates, vector or comma-separated; default: NULL.
  # Standard mode: mice imputes missing values; PC count follows covariate count.
  # Integrated modes: complete cases; explicit association_n_pcs controls genomic PCs.
  covariate_pheno = NULL,
  # Microbiome adjustment of the response; default: FALSE. Keep FALSE in new mode.
  covariate_metag = FALSE,
  # Retain markers with MAF >= this value. Default: 0.02.
  maf = 0.02,
  # Rebuild kinship leaving out the tested chromosome. Default: FALSE; use with MLM.
  LOCO = FALSE,
  # GWASpoly workers. Default: 1. Integrated scans and sPLS bootstraps are serial.
  # Larger values do not parallelise these branches or set BLAS thread counts.
  cores = 1,
  # Plot FDR (0.05) / Bonferroni (0.05) lines. Defaults: TRUE.
  # Tables are still calculated when these plot switches are FALSE.
  fdr = TRUE,
  bonferroni = TRUE,
  # Suggestive -log10(p) cutoff; legacy table first filters at FDR 0.5. Default: 5.
  suggestive = 5,
  # Standard-mode permutations; >=100 activates, 0/1 skips. Default: 1.
  # Integrated modes reject perm>=100; their F-test p-values are approximate P3D tests.
  perm = 1,
  
  # LEGACY COMPATIBILITY: leave these NULL when using the new controls above.
  # Legacy proxy training response override. Default: NULL. Prefer microbial_response.
  alternate_trait = NULL,
  # Legacy taxa/"auto"/"auto-null" proxy selector. Default: c("NULL").
  # Prefer microbial_features plus microbial_method for new analyses.
  trait_microbial_proxy = NULL,
  # Legacy separate strain/species scans (not composite builders). Defaults: NULL.
  pheno_taxa_strain = NULL,
  pheno_taxa_species = NULL,
  # Legacy strain-table input. Default: NULL.
  metag_data_strains = NULL
)

# EDIT THE CALL ABOVE FOR THESE COMMON ANALYSES (one outcome per call):
# 1. Fv abundance, no FER:
#    microbial_features = "Fusarium_verticillioides", microbial_method = "single",
#    microbial_trait_name = "Fv", phenofile = NULL, microbial_response = NULL.
#    Use microbial_transform = "none" for input abundance, "relative" for proportions,
#    or the default "auto" for CLR abundance.
# 2. Fusarium total abundance, no FER:
#    microbial_features = "Fusarium*", microbial_method = "sum",
#    microbial_trait_name = "Fusarium_total". Default transform gives summed proportions.
# 3. Fusarium PC1, no FER: the active configuration above.
#    For PC2/PC3, set microbial_method = "pc2"/"pc3" and change projname and
#    microbial_trait_name accordingly. The requested component must have nonzero variance.
# 4. FER-informed Fusarium proxy:
#    phenofile = "traits.tsv", microbial_features = "Fusarium*",
#    microbial_method = "spls", microbial_response = "FER",
#    microbial_trait_name = "Fusarium_FER_proxy".
#    Use microbial_features = "all" for the full microbiome, or c("species1", ...).
# 5. Learn/export a proxy separately: use configuration 4 with run_gwas = FALSE.
# 6. Apply exported weights without FER and run GWAS:
#    microbial_features = NULL, microbial_method = "weighted", microbial_response = NULL,
#    phenofile = NULL, microbial_weights = "training_microbial_traits/Fusarium_FER_proxy_weights.tsv".
#    Supply the complete weights file, including zero-weight reference taxa.
# 7. Ordinary FER GWAS: microbial_features = NULL, microbial_weights = NULL,
#    microbial_method = "auto", microbial_response = NULL, metag_data_species = NULL,
#    phenofile = "traits.tsv", trait_names = "FER".
# 8. Joint FER + Fusarium PC1 GWAS:
#    phenofile = "traits.tsv", trait_names = "FER", microbial_features = "Fusarium*",
#    microbial_method = "pc1", microbial_trait_name = "Fusarium_PC1",
#    association_mode = "joint", method = "MLM", model_effect = c("Add", "Dom").
# 9. SNP x Fusarium score: configuration 8 with association_mode = "interaction".
#    Use "both" for both analyses. Supply several trait_names for separate pairs.
#    LOCO=TRUE excludes the tested chromosome from kinship; use a distinct projname.
#    Prefer unsupervised scores; same-cohort sPLS weights yield exploratory p-values.
