#!/usr/bin/env Rscript
# Source the updated pipeline, then edit the input filenames below.

# load systemsGWAS directly from  github
systemsGWAS_url <- paste0("https://raw.githubusercontent.com/","bodeolukolu/systemsGWAS/main/R/systemsGWAS.R")
tryCatch(source(systemsGWAS_url),error = function(e) {stop("Unable to source systemsGWAS.R: ", conditionMessage(e))})

# load systemsGWAS directory, if downloaded
# source("systemsGWAS.R")


systemsGWAS(
  wdir = "./",
  projname = "Fusarium_PC1",
  
  ploidy_levels = 2,
  genofile_2x = "Zm_2x_rd6_maf0.02_dose.txt",

  phenofile = NULL,
  trait_names = NULL,

  metag_data_species = "metag_species.tsv",
  microbial_features = "Fusarium*",           # Exact species, vector of names, "Fusarium*", "Fusarium_spp", or "all".
  microbial_method = "pc1",                   # "single", "sum", "PC1", "PC2", "PC3", or "spls" (specify only one at a time)
  microbial_response = NULL,                  # specify phenotype spls uses e.g. "FER"
  microbial_weights = NULL,                   # Full *_weights.tsv from an earlier run; method must be "weighted". Default: NULL.
  microbial_trait_name = "composite_Fusarium_PC1",
  microbial_transform = "auto",               # "none", "relative", "clr" or "log1p". Default "auto": relative for sum, clr otherwise.
  microbial_prevalence = 0.05,
  microbial_bootstraps = 100L,
  microbial_stability = 0.30,
  microbial_seed = 123L,
  microbial_pseudocount = 1e-6,
  run_gwas = TRUE,
  
  # INTEGRATED TRAIT + MICROBIAL/PHENOMIC ASSOCIATION
  # "gwaspoly-marginal": single trait scan
  # "sgwas-marginal" (default):single trait scan
  # `sgwas-multitrait`: two or more traits
  # "sgwas-joint": two-trait mixed model
  # "sgwas-interaction": SNP x microbial score
  # "sgwas-all": run both, with matched baselines.
  # New modes require phenofile + trait_names and one constructed microbial score.
  # Use method="MLM", model_effect="Add", covariate_metag=FALSE below.
  association_mode = "sgwas-marginal",             # "gwaspoly-marginal", "sgwas-marginal", "sgwas-joint", "sgwas-interaction", "sgwas-all"
  association_n_pcs = 0L,
  association_max_missing = 0.50,
  association_chunk_size = 256L,
  association_maxit = 300L,
  association_geno_freq=0.99,
  association_output_style = "verbose",
  association_plot_formats = c("tiff"),
  
  multitrait_missing = "impute",                    # "impute" or "complete_case".
  multitrait_missing_threshold = 0.30,              # >30% overall OR in any trait: drop incomplete rows.

  method = "MLM",
  model_effect = c("Add","Dom"),
  covariate_pheno = NULL,
  covariate_metag = FALSE,
  maf = 0.02,
  LOCO = FALSE,
  cores = 3,
  fdr = TRUE,
  bonferroni = TRUE,
  suggestive = 5,
  perm = 1,
  
  alternate_trait = NULL,
  trait_microbial_proxy = NULL,
  pheno_taxa_strain = NULL,
  pheno_taxa_species = NULL,
  metag_data_strains = NULL
)
