#!/usr/bin/env Rscript
# systemsGWAS: microbial trait construction and genome-wide association.

# Microbial phenotype construction. All selection and scoring uses sample IDs.
.sgwas_names <- function(x) {
  if(is.null(x)) return(NULL)
  x <- trimws(unlist(strsplit(as.character(x),",",fixed=TRUE),use.names=FALSE))
  x <- unique(x[!is.na(x) & nzchar(x) & x != "NULL"])
  if(length(x)) x else NULL
}

.sgwas_numeric <- function(x,label) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  if(any(is.na(y) & !is.na(x) & nzchar(trimws(as.character(x)))))
    stop("Non-numeric values in ",label,call.=FALSE)
  y
}

.sgwas_read_microbiome <- function(file) {
  d <- read.delim(file,check.names=FALSE,stringsAsFactors=FALSE,colClasses="character")
  if(ncol(d)<2L || !nrow(d) || anyDuplicated(names(d)))
    stop("Microbiome input needs unique column names and at least one feature/sample")
  sample_key <- intersect(c("Plant_ID","SampleID","Sample_ID","Taxa"),names(d))
  feature_key <- intersect(c("species","Species","Feature","feature","strain","Strain"),names(d))
  if(length(sample_key)) {
    key <- sample_key[1L]; ids <- trimws(as.character(d[[key]]))
    cols <- setdiff(names(d),key)
    X <- as.matrix(as.data.frame(lapply(cols,function(nm) .sgwas_numeric(d[[nm]],nm))))
    dimnames(X) <- list(ids,cols)
  } else if(length(feature_key)) {
    key <- feature_key[1L]; features <- trimws(as.character(d[[key]]))
    metadata <- c("species","strain","genus","family","order","class","phylum",
                  "kingdom","domain","tax_id","taxid","taxonomy","feature")
    cols <- names(d)[!tolower(names(d)) %in% metadata]
    if(!length(cols)) stop("No sample abundance columns in microbiome input")
    values <- as.matrix(as.data.frame(lapply(cols,function(nm) .sgwas_numeric(d[[nm]],nm))))
    X <- t(values); dimnames(X) <- list(cols,features)
  } else stop("Microbiome input needs a species/Feature column (features in rows) or Plant_ID/SampleID/Taxa column (samples in rows)")
  if(anyNA(rownames(X)) || anyNA(colnames(X)) || any(!nzchar(rownames(X))) ||
     any(!nzchar(colnames(X))) || anyDuplicated(rownames(X)) || anyDuplicated(colnames(X)))
    stop("Microbiome sample IDs and feature names must be nonempty and unique")
  storage.mode(X) <- "double"
  if(any(!is.finite(X)) || any(X<0))
    stop("Microbiome abundances must be finite, nonnegative numbers; use zero for measured absence, not missing measurements")
  X
}

.sgwas_select_features <- function(features, available) {
  selectors <- .sgwas_names(features)
  if(is.null(selectors)) stop("Supply microbial_features: exact names, a prefix ending in *, Genus_spp, or all")
  if(identical(selectors,"all")) return(available)
  selected <- character()
  for(s in selectors) {
    if(grepl("(_spp| spp)$",s)) {
      prefix <- sub("(_spp| spp)$","",s)
      hit <- available[startsWith(available,paste0(prefix,"_")) | startsWith(available,paste0(prefix," "))]
    } else if(endsWith(s,"*")) {
      prefix <- substr(s,1L,nchar(s)-1L)
      if(!nzchar(prefix)) stop("Use all to select the complete microbiome")
      hit <- available[startsWith(available,prefix)]
    } else hit <- intersect(s,available)
    if(!length(hit)) stop("No microbiome features match: ",s)
    selected <- union(selected,hit)
  }
  available[available %in% selected]
}

.sgwas_transform <- function(X,method,pseudocount) {
  if(method=="none") return(X)
  if(method=="log1p") return(log1p(X))
  total <- rowSums(X)
  if(any(total<=0)) stop("CLR/relative transformation requires a positive total abundance in every sample")
  relative <- X/total
  if(method=="relative") return(relative)
  if(method=="clr") {
    # Close to proportions before adding the pseudocount. Reference columns
    # are fixed in the exported weights, including zero-weight reference taxa.
    logx <- log(relative+pseudocount)
    return(logx-rowMeans(logx))
  }
  stop("Unknown microbial transformation: ",method)
}

.sgwas_phenotypes <- function(file) {
  if(is.null(file)) return(NULL)
  d <- read.delim(file,check.names=FALSE,stringsAsFactors=FALSE,colClasses="character")
  if(ncol(d)<2 || anyDuplicated(names(d))) stop("Phenotype table needs unique trait columns and an ID column")
  names(d)[1L] <- "Plant_ID"
  d$Plant_ID <- trimws(as.character(d$Plant_ID))
  if(anyNA(d$Plant_ID) || any(!nzchar(d$Plant_ID))) stop("Phenotype IDs cannot be missing/blank")
  d
}

.sgwas_id_values <- function(d,column,ids) {
  if(is.null(d) || !column %in% names(d)) stop("Phenotype column not found: ",column)
  values <- .sgwas_numeric(d[[column]],column)
  groups <- split(values,d$Plant_ID)
  means <- vapply(groups,function(x) {x <- x[is.finite(x)];if(length(x)) mean(x) else NA_real_},numeric(1))
  unname(means[match(ids,names(means))])
}

.sgwas_load_weights <- function(file) {
  w <- read.delim(file,check.names=FALSE,stringsAsFactors=FALSE)
  required <- c("Feature","Weight","Center","Scale","Transform","Pseudocount","ReferenceSetSize")
  if(!all(required %in% names(w)) || !nrow(w))
    stop("microbial_weights must be a complete *_weights.tsv exported by this pipeline")
  if(anyNA(w$Feature) || any(!nzchar(w$Feature)) || anyDuplicated(w$Feature)) stop("Weight features must be unique")
  for(nm in c("Weight","Center","Scale","Pseudocount","ReferenceSetSize")) w[[nm]] <- .sgwas_numeric(w[[nm]],nm)
  if(any(!is.finite(as.matrix(w[c("Weight","Center","Scale","Pseudocount","ReferenceSetSize")])))) stop("Weight metadata must be finite")
  if(any(w$Scale<=0) || !any(w$Weight!=0) || length(unique(w$Transform))!=1L ||
     length(unique(w$Pseudocount))!=1L || any(w$Pseudocount<=0)) stop("Invalid transformation or scaling in weights")
  if(any(w$ReferenceSetSize!=nrow(w))) stop("Incomplete weights table: use the full *_weights.tsv, not *_selected_weights.tsv")
  w
}

# Public preparation helper, also called automatically by systemsGWAS().
# The score is sum_j Weight_j * (transformed_abundance_j - Center_j) / Scale_j.
systemsGWAS_prepare_microbiome <- function(
    microbiome_file, features="all", method="auto", response=NULL, phenofile=NULL,
    weights_file=NULL, trait_name="Microbial_composite", transform="auto",
    prevalence=0.05, bootstraps=100L, stability=0.60, seed=123L,
    pseudocount=1e-6, outdir="microbial_traits", covariates=NULL) {
  method <- match.arg(tolower(method),c("auto","single","sum","pc1","pc2","pc3","spls","weighted"))
  transform <- match.arg(tolower(transform),c("auto","none","relative","clr","log1p"))
  if(length(trait_name)!=1L || is.na(trait_name) || !grepl("^[A-Za-z][A-Za-z0-9_]*$",trait_name) || trait_name=="Plant_ID")
    stop("microbial_trait_name must start with a letter and contain only letters, digits and underscores")
  check_fraction <- function(x,name) if(length(x)!=1 || !is.numeric(x) || !is.finite(x) || x<0 || x>1) stop(name," must be between 0 and 1")
  check_fraction(prevalence,"microbial_prevalence");check_fraction(stability,"microbial_stability")
  if(length(bootstraps)!=1L || !is.numeric(bootstraps) || !is.finite(bootstraps) || bootstraps<1 || bootstraps!=floor(bootstraps)) stop("microbial_bootstraps must be a positive integer")
  if(length(seed)!=1L || !is.numeric(seed) || !is.finite(seed) || seed<0 || seed>.Machine$integer.max || seed!=floor(seed)) stop("microbial_seed must be a nonnegative integer")
  if(length(pseudocount)!=1L || !is.numeric(pseudocount) || !is.finite(pseudocount) || pseudocount<=0) stop("microbial_pseudocount must be positive")
  covariates <- .sgwas_names(covariates)
  if(trait_name %in% covariates) stop("The GWAS response cannot be its own covariate")
  response <- .sgwas_names(response)
  if(method=="spls" && length(response)!=1L) stop("sPLS requires one microbial_response phenotype column")
  if(method!="spls" && length(response)) stop("microbial_response is used only with microbial_method = spls")
  if(method=="weighted" && is.null(weights_file)) stop("weighted requires microbial_weights")
  if(method!="weighted" && !is.null(weights_file)) stop("Set microbial_method = weighted to reuse microbial_weights")
  X <- .sgwas_read_microbiome(microbiome_file)
  ph <- if(method=="spls" || length(covariates)) .sgwas_phenotypes(phenofile) else NULL
  ids <- rownames(X); train <- rep(TRUE,nrow(X)); Y <- NULL
  pca_methods <- c("pc1","pc2","pc3"); pca_component <- NA_integer_
  variance_explained <- NA_real_; success <- 0L; bootstrap_log <- data.frame()
  source_method <- method; source_response <- if(length(response)) response else ""
  if(method=="weighted") {
    w <- .sgwas_load_weights(weights_file)
    absent <- setdiff(w$Feature,colnames(X))
    if(length(absent)) stop("Microbiome input is missing saved reference features: ",paste(absent,collapse=", "))
    X <- X[,w$Feature,drop=FALSE]
    transform <- w$Transform[1];pseudocount <- w$Pseudocount[1]
    source_method <- if("Method" %in% names(w)) as.character(w$Method[1]) else "weighted"
    source_response <- if("Response" %in% names(w)) as.character(w$Response[1]) else ""
    if(source_method %in% pca_methods) {
      pca_component <- match(source_method,pca_methods)
      if("PC_Variance_Explained" %in% names(w))
        variance_explained <- as.numeric(w$PC_Variance_Explained[1L])
    }
    w$Selected <- w$Weight!=0
    T <- .sgwas_transform(X,transform,pseudocount)
  } else {
    candidates <- .sgwas_select_features(features,colnames(X))
    if(method=="auto") method <- if(length(candidates)==1L) "single" else "pc1"
    if(method=="single" && length(candidates)!=1L) stop("single requires exactly one selected species")
    if(method=="spls") {
      Y <- .sgwas_id_values(ph,response,ids); train <- is.finite(Y)
      if(sum(train)<5L || !is.finite(sd(Y[train])) || sd(Y[train])==0) stop("sPLS requires at least five observed training samples and a variable response")
    }
    frequency <- colMeans(X[train,candidates,drop=FALSE]>0)
    candidates <- candidates[frequency>=prevalence]
    if(!length(candidates)) stop("No requested features pass microbial_prevalence")
    if(transform=="auto") transform <- if(method=="sum") "relative" else "clr"
    if(method=="sum" && !transform %in% c("none","relative")) stop("sum requires microbial_transform = relative or none (sum abundances, not CLR values)")
    T <- .sgwas_transform(X,transform,pseudocount)
    w <- data.frame(Feature=colnames(X),Weight=0,Center=0,Scale=1,Selected=FALSE,
                    Candidate=colnames(X) %in% candidates,Selection_Frequency=NA_real_,
                    Transform=transform,Pseudocount=pseudocount,ReferenceSetSize=ncol(X),
                    Method=method,Response=source_response,stringsAsFactors=FALSE)
    if(method %in% c("single","sum")) w$Weight[match(candidates,w$Feature)] <- 1
    else {
      SD <- apply(T[train,candidates,drop=FALSE],2,sd)
      candidates <- candidates[is.finite(SD) & SD>sqrt(.Machine$double.eps)]
      if(!length(candidates)) stop("No variable microbial candidates remain after transformation")
      centers <- colMeans(T[train,candidates,drop=FALSE]);scales <- apply(T[train,candidates,drop=FALSE],2,sd)
      if(method %in% pca_methods) {
        pca_component <- match(method,pca_methods)
        fit <- prcomp(T[train,candidates,drop=FALSE],center=TRUE,scale.=TRUE)
        # Check numerical rank, not just column count: CLR and correlated
        # features can leave a requested component with zero variance.
        rank_tol <- max(fit$sdev)*sqrt(.Machine$double.eps)
        effective_rank <- min(sum(fit$sdev>rank_tol),nrow(T[train,,drop=FALSE])-1L)
        if(pca_component>effective_rank || pca_component>ncol(fit$rotation))
          stop("Requested PC",pca_component," cannot be estimated; effective PCA rank is ",
               effective_rank," after filtering and transformation")
        loading <- fit$rotation[,pca_component]
        if(loading[which.max(abs(loading))]<0) loading <- -loading
        variance_explained <- fit$sdev[pca_component]^2/sum(fit$sdev^2)
        chosen <- candidates
      } else {
        # One species needs no sparse selection; 2-4 predictors remain below
        # the five-candidate minimum for the multi-species bootstrap screen.
        if(length(candidates)>1L && length(candidates)<5L)
          stop("Bootstrap sPLS requires at least 5 variable candidate species after filtering; found ",length(candidates))
        if(!requireNamespace("mixOmics",quietly=TRUE)) stop("Install mixOmics to construct an sPLS proxy")
        had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
        if(had_seed) old_seed <- get(".Random.seed",envir=.GlobalEnv)
        on.exit({if(had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv)
                 else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)},add=TRUE)
        set.seed(as.integer(seed))
        grid <- if(length(candidates)==1L) 1L else c(5L,10L,20L,30L,50L,100L,200L)
        grid <- grid[grid<=length(candidates)]
        counts <- setNames(integer(length(candidates)),candidates)
        logs <- vector("list",as.integer(bootstraps)); Xtrain <- T[train,candidates,drop=FALSE];Ytrain <- Y[train]
        for(b in seq_len(as.integer(bootstraps))) {
          idx <- sample.int(nrow(Xtrain),nrow(Xtrain),replace=TRUE)
          keep <- grid[sample.int(length(grid),1L)]
          fit <- tryCatch({
            if(sd(Ytrain[idx])==0 || any(apply(Xtrain[idx,,drop=FALSE],2,sd)==0)) stop("constant bootstrap response or predictor")
            Xboot <- Xtrain[idx,,drop=FALSE]
            # Resampling repeats biological IDs; mixOmics requires unique row names.
            rownames(Xboot) <- NULL
            mixOmics::spls(Xboot,matrix(Ytrain[idx],ncol=1),ncomp=1,keepX=keep,scale=TRUE,mode="regression")
          },error=function(e)e)
          selected <- character(); problem <- ""
          if(inherits(fit,"error")) problem <- conditionMessage(fit)
          else {
            load <- fit$loadings$X[,1L]
            if(length(load)!=length(candidates) || any(!is.finite(load))) problem <- "invalid bootstrap loadings"
            else selected <- candidates[abs(load)>0]
            if(!length(selected) && !nzchar(problem)) problem <- "no selected species"
          }
          if(!nzchar(problem)) {counts[selected] <- counts[selected]+1L;success <- success+1L}
          logs[[b]] <- data.frame(Bootstrap=b,KeepX=keep,Success=!nzchar(problem),N_Selected=length(selected),Message=problem)
        }
        bootstrap_log <- do.call(rbind,logs)
        if(!dir.exists(outdir) && !dir.create(outdir,recursive=TRUE)) stop("Cannot create microbial output directory")
        write.table(bootstrap_log,file.path(outdir,paste0(trait_name,"_bootstrap.tsv")),
                    sep="\t",quote=FALSE,row.names=FALSE,na="NA")
        if(!success) stop("No successful sPLS bootstrap fits; see the saved bootstrap log")
        if(success<bootstraps) message("[systemsGWAS] sPLS bootstraps succeeded: ",success,"/",bootstraps,
                                     "; selection frequency uses all attempted bootstraps.")
        proportions <- counts/as.integer(bootstraps)
        w$Selection_Frequency[match(candidates,w$Feature)] <- proportions
        chosen <- names(proportions)[proportions>=stability & counts>0]
        if(!length(chosen)) stop("No species pass microbial_stability; inspect data or choose a lower threshold")
        fit <- mixOmics::spls(T[train,chosen,drop=FALSE],matrix(Ytrain,ncol=1),ncomp=1,
                              keepX=length(chosen),scale=TRUE,mode="regression")
        loading <- fit$loadings$X[,1L];names(loading) <- chosen
        if(any(!is.finite(loading))) stop("Non-finite final sPLS loadings")
        component <- as.numeric(scale(T[train,chosen,drop=FALSE],center=centers[chosen],scale=scales[chosen]) %*% loading)
        # Verify that exported weights reproduce the actual fitted X component.
        if(!isTRUE(all.equal(component,as.numeric(fit$variates$X[,1L]),tolerance=1e-7)))
          stop("sPLS score reconstruction does not match fitted component; check mixOmics version")
        orientation <- cor(component,Ytrain,method="spearman")
        if(is.finite(orientation) && orientation<0) loading <- -loading
      }
      ix <- match(chosen,w$Feature)
      w$Weight[ix] <- loading;w$Center[ix] <- centers[chosen];w$Scale[ix] <- scales[chosen]
    }
    w$Selected <- w$Weight!=0;source_method <- method
  }
  # Preserve the chosen PCA axis and training variance fraction for reuse.
  # Older exported weights without these columns still produce identical scores.
  w$PC_Component <- pca_component
  w$PC_Variance_Explained <- variance_explained
  Z <- sweep(sweep(T,2,w$Center,"-"),2,w$Scale,"/")
  score <- as.numeric(Z %*% w$Weight)
  if(any(!is.finite(score)) || !is.finite(sd(score)) || sd(score)==0) stop("The constructed microbial trait is non-finite or constant")
  composite <- data.frame(Plant_ID=ids,Value=score,check.names=FALSE);names(composite)[2L] <- trait_name
  gwas_pheno <- composite
  for(nm in covariates) gwas_pheno[[nm]] <- .sgwas_id_values(ph,nm,ids)
  # Export preparation separately from the phenotype table used by GWAS.
  if(!dir.exists(outdir) && !dir.create(outdir,recursive=TRUE)) stop("Cannot create microbial output directory")
  prefix <- file.path(outdir,trait_name)
  paths <- c(composite=paste0(prefix,"_composite.tsv"),weights=paste0(prefix,"_weights.tsv"),
             selected_weights=paste0(prefix,"_selected_weights.tsv"),
             phenotype=paste0(prefix,"_gwas_phenotypes.tsv"),summary=paste0(prefix,"_summary.tsv"),
             model=paste0(prefix,"_model.rds"))
  write_tsv <- function(d,p) write.table(d,p,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
  write_tsv(composite,paths['composite']);write_tsv(w,paths['weights'])
  write_tsv(w[w$Selected,,drop=FALSE],paths['selected_weights']);write_tsv(gwas_pheno,paths['phenotype'])
  if(nrow(bootstrap_log)) {
    paths <- c(paths,bootstrap=paste0(prefix,"_bootstrap.tsv"));write_tsv(bootstrap_log,paths['bootstrap'])
  }
  summary <- data.frame(Trait=trait_name,Method=method,Source_Method=source_method,
                         Response=source_response,Transform=transform,N_Samples=nrow(X),
                         N_Training=if(method=="weighted") NA_integer_ else sum(train),
                         N_Reference=ncol(X),N_Selected=sum(w$Selected),
                         Bootstrap_Attempts=if(method=="spls") as.integer(bootstraps) else 0L,
                         Bootstrap_Successes=success,Stability=if(method=="spls") stability else NA_real_,
                         PC_Component=pca_component,PC_Variance_Explained=variance_explained,
                         PC1_Variance_Explained=if(isTRUE(pca_component==1L)) variance_explained else NA_real_,
                         Training_Spearman=if(method=="spls") cor(score[train],Y[train],method="spearman") else NA_real_)
  write_tsv(summary,paths['summary'])
  model <- list(weights=w,summary=summary,training_ids=if(method=="weighted") NULL else ids[train],
                requested_features=features,prevalence=prevalence,seed=seed,
                score_formula="sum(Weight * (transformed abundance - Center) / Scale)",
                mixOmics_version=if(method=="spls") as.character(utils::packageVersion("mixOmics")) else NULL)
  saveRDS(model,paths['model'])
  message("[systemsGWAS] Constructed ",trait_name," using ",method,": ",sum(w$Selected)," features; ",nrow(X)," samples. Outputs: ",outdir)
  invisible(list(phenotype=gwas_pheno,weights=w,summary=summary,paths=paths,model=model))
}

.systemsGWAS_engine <- function(
    wdir = "./",
    projname = "GWAS",
    ploidy_levels = c(2,4,6,8),
    trait_names = c("trait1","trait2"),
    alternate_trait=NULL,
    trait_microbial_proxy=c("NULL"),
    model_effect = c("Add","Dom"),
    fdr = TRUE,
    bonferroni = TRUE,
    suggestive = 5,
    perm = 1,
    cores = 1,
    genofile_2x = NULL,
    genofile_4x = NULL,
    genofile_6x = NULL,
    genofile_8x = NULL,
    phenofile = NULL,
    method = c("GLM","MLM"),
    covariate_pheno = NULL,
    covariate_metag = FALSE,
    maf = 0.02,
    LOCO = FALSE,
    pheno_taxa_strain = NULL,
    pheno_taxa_species = NULL,
    metag_data_strains = NULL,
    metag_data_species = NULL
) {
  load_packages <- function(pkgs) {
    for (p in pkgs) {
      if (!requireNamespace(p, quietly = TRUE)) {
        install.packages(p, repos = "https://cloud.r-project.org")
      }
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }

  # List of required packages
  pkgs <- c("GWASpoly", "qqplotr", "AGHmatrix", "corrplot", "ggcorrplot", "ggplot2",
            "plyr", "dplyr", "tidyr", "car", "MASS", "Hmisc", "data.table", "stringr",
            "heatmaply", "ppcor", "zoo", "GGally", "reshape2", "compositions",
            "sommer", "mice", "qvalue", "mixOmics")

  # Call it once at the start of your function
  load_packages(pkgs)

  #############################################################################################################################################################################
  # Specify parameters
  ####################
  previous_wdir <- getwd()
  on.exit(setwd(previous_wdir), add = TRUE)
  wdir <- normalizePath(wdir, mustWork = TRUE)
  input_path <- function(path) .sgwas_path(path, wdir)
  setwd(wdir)
  dir_prefix_name <- projname
  pop <- projname
  ploidy_levels <- as.numeric(ploidy_levels)
  gwas_method <- .sgwas_names(method)
  maf <- as.numeric(maf)
  file_ploidy_2 <- genofile_2x
  file_ploidy_4 <- genofile_4x
  file_ploidy_6 <- genofile_6x
  file_ploidy_8 <- genofile_8x
  phenotype_data <- phenofile
  alternate_trait=alternate_trait
  covariatename <- .sgwas_names(covariate_pheno)
  trait_microbial_proxy <- .sgwas_names(trait_microbial_proxy)
  covariate <- covariate_metag
  corr_coeff <- if(covariate_metag) "full" else NULL
  # Microbiome as phenotypes
  taxa_strain <- pheno_taxa_strain
  taxa_species <- pheno_taxa_species
  metagenome_data_strains <- metag_data_strains
  metagenome_data_species <- metag_data_species
  if (!is.null(model_effect)) {
    model_effect <- .sgwas_names(model_effect)
    if ("Add" %in% model_effect && "Dom" %in% model_effect) {model <- "Add_Dom"}
    if ("Add" %in% model_effect && !("Dom" %in% model_effect)) {model <- "Add"}
    if (!("Add" %in% model_effect) && "Dom" %in% model_effect) {model <- "Add_Dom"}
  } else {
    model <- "Add"
  }
  # Thresholds
  threshold_FDR <- fdr
  threshold_Bonferroni <- bonferroni
  threshold_suggestive <- as.numeric(suggestive)
  permutations <- as.numeric(perm)
  cores <- as.numeric(cores)


  #############################################################################################################################################################################
  # Specify other variables for metagenome-based analysis
  ##################################################
  mincorr <- 0                                 # minimum correlation coefficient (trait vs subsets of metagenome)
  maxcorr <- 0.8                                 # maximum correlation coefficient (trait vs subsets of metagenome)
  pvalue <- 1
  perc <- 5
  metag_method <- "Aitchison"            # Aitchison or pca
  #############################################################################################################################################################################
  # Optional parameters
  ######################
  diversity=TRUE                                         # TRUE (default) or FALSE
  biparental=TRUE                                        # Generates a QTL profile manhattan plots
  # To permutation test, number of test must be >= 100. If time is not a constraint, at least 1,000 test is suggested
  if ("2" %in% ploidy_levels) { Gmethod="VanRaden" }
  if ("4" %in% ploidy_levels) { Gmethod="VanRaden" }
  if ("6" %in% ploidy_levels) { Gmethod="VanRaden" }
  if ("8" %in% ploidy_levels) { Gmethod="VanRaden" }

  #############################################################################################################################################################################
  if (!is.null(trait_names) == (!is.null(taxa_strain) || !is.null(taxa_species))) {
    print("Please provide <trait names> or <taxa names>. YOu can't provide both")
  } else {
    for (gwas_method in gwas_method) {
      for (ploidy in ploidy_levels) {
        metag <- NULL
        if(covariate == TRUE) {no_covariate <- FALSE}
        if(covariate == FALSE) {no_covariate <- TRUE}
        if(!is.null(metagenome_data_strains)) {metag <- read.table(input_path(metagenome_data_strains), header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)}
        if(!is.null(metagenome_data_species)) {metag <- read.table(input_path(metagenome_data_species), header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)}

        if (no_covariate == TRUE && is.null(covariatename) && !is.null(trait_names)) { dir1 <- paste(dir_prefix_name,"_pheno_Mnocov_",ploidy,"x_", gwas_method,  sep="") } else {dir1 <- NULL}
        if (no_covariate == TRUE && !is.null(covariatename) && !is.null(trait_names)) { dir2 <- paste(dir_prefix_name,"_pheno_Pcov_Mnocov_",ploidy,"x_", gwas_method,  sep="") } else {dir2 <- NULL}
        if (covariate == TRUE && is.null(covariatename) && !is.null(metagenome_data_strains) && !is.null(trait_names) && is.null(taxa_strain)) { dir3 <- paste(dir_prefix_name,"_pheno_strain_Mcov_",ploidy,"x_", gwas_method,  sep="") } else {dir3 <- NULL}
        if (covariate == TRUE && is.null(covariatename) && !is.null(metagenome_data_species) && !is.null(trait_names) && is.null(taxa_species)) { dir4 <- paste(dir_prefix_name,"_pheno_species_Mcov_",ploidy,"x_", gwas_method,  sep="") } else {dir4 <- NULL}
        if (covariate == TRUE && !is.null(covariatename) && !is.null(metagenome_data_strains) && !is.null(trait_names) && is.null(taxa_strain)) { dir5 <- paste(dir_prefix_name,"_pheno_Pcov_strain_Mcov_",ploidy,"x_", gwas_method,  sep="") } else {dir5 <- NULL}
        if (covariate == TRUE && !is.null(covariatename) && !is.null(metagenome_data_species) && !is.null(trait_names) && is.null(taxa_species)) { dir6 <- paste(dir_prefix_name,"_pheno_Pcov_species_Mcov_",ploidy,"x_", gwas_method,  sep="") } else {dir6 <- NULL}
        if (no_covariate == TRUE && is.null(covariatename) && !is.null(metagenome_data_strains) && !is.null(taxa_strain) && is.null(trait_names)) { dir7 <- paste(dir_prefix_name,"_strain","_Mnocov_",ploidy,"x_", gwas_method,  sep="") } else {dir7 <- NULL}
        if (no_covariate == TRUE && is.null(covariatename) && !is.null(metagenome_data_species) && !is.null(taxa_species) && is.null(trait_names)) { dir8 <- paste(dir_prefix_name,"_species","_Mnocov_",ploidy,"x_", gwas_method,  sep="") } else {dir8 <- NULL}
        if (no_covariate == TRUE && is.null(covariatename) && !is.null(covariatename) && !is.null(metagenome_data_strains) && !is.null(taxa_strain) && is.null(trait_names)) { dir9 <- paste(dir_prefix_name,"Pcov_strain","_Mnocov_",ploidy,"x_", gwas_method,  sep="") } else {dir9 <- NULL}
        if (no_covariate == TRUE && is.null(covariatename) && !is.null(covariatename) && !is.null(metagenome_data_species) && !is.null(taxa_species) && is.null(trait_names)) { dir10 <- paste(dir_prefix_name,"Pcov_species","_Mnocov_",ploidy,"x_", gwas_method,  sep="") } else {dir10 <- NULL}
        if (covariate == TRUE && !is.null(metagenome_data_strains) && !is.null(taxa_strain) && is.null(trait_names)) { dir11 <- paste(dir_prefix_name,"_strain_Mcov_",ploidy,"x_", gwas_method,  sep="") } else {dir11 <- NULL}
        if (covariate == TRUE && !is.null(metagenome_data_species) && !is.null(taxa_species) && is.null(trait_names)) { dir12 <- paste(dir_prefix_name,"_species_Mcov_",ploidy,"x_", gwas_method,  sep="") } else {dir12 <- NULL}
        if (covariate == TRUE && !is.null(covariatename) && !is.null(metagenome_data_strains) && !is.null(taxa_strain) && is.null(trait_names)) { dir13 <- paste(dir_prefix_name,"Pcov_strain_Mcov_",ploidy,"x_", gwas_method,  sep="") } else {dir13 <- NULL}
        if (covariate == TRUE && !is.null(covariatename) && !is.null(metagenome_data_species) && !is.null(taxa_species) && is.null(trait_names)) { dir14 <- paste(dir_prefix_name,"Pcov_species_Mcov_",ploidy,"x_", gwas_method,  sep="") } else {dir14 <- NULL}
        dir_list <- c(dir1,dir2,dir3,dir4,dir5,dir6,dir7,dir8,dir9,dir10,dir11,dir12,dir13,dir14)

        for(job in c(dir_list)) {
          dir.create(job)
          setwd(job)
          dir.create("correlations")
          dir.create("scores_effects")
          dir.create("sigFDR")
          dir.create("sigBonferroni")
          dir.create("sigSuggestive")
          dir.create("sigpermute")
          dir.create("qqplots")
          dir.create("manplots")
          dir.create("normal_distribution")

          if(grepl("_strain_",job) && !grepl("_pheno_",job)) {taxa_level <- "strain"}
          if(grepl("_species_",job) && !grepl("_pheno_",job)) {taxa_level <- "species"}
          if(grepl("_pheno_",job)) {taxa_level <- "notaxa"}

          if("2" %in% ploidy_levels && ploidy == 2){genotype_data <- file_ploidy_2}
          if("4" %in% ploidy_levels && ploidy == 4){genotype_data <- file_ploidy_4}
          if("6" %in% ploidy_levels && ploidy == 6){genotype_data <- file_ploidy_6}
          if("8" %in% ploidy_levels && ploidy == 8){genotype_data <- file_ploidy_8}

          if(!is.null(metag) && covariate_metag == TRUE) {trait_associated_covariates <- TRUE} else {trait_associated_covariates <- FALSE}
          metagenome_data <- NULL
          if(grepl("_strain_",job)) {metagenome_data <- metagenome_data_strains}
          if(grepl("_species_",job)) {metagenome_data <- metagenome_data_species}
          if(grepl("pheno_nocov",job)) {metagenome_data <- NULL}

          #############################################################################################################################################################################
          # organize and generate required files
          ######################################
          normalize_kinmat <- function(kinmat){
            #normalize kinship so that Kij \in [0,1]
            tmp=kinmat - min(kinmat)
            tmp=tmp/max(tmp)
            tmp[1:9,1:9]
            #fix eigenvalues to positive
            diag(tmp)=diag(tmp)-min(eigen(tmp)$values)
            tmp[1:9,1:9]
            return(tmp)
          }
          geno <- read.table(input_path(genotype_data), header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
          geno <- geno %>%
            distinct(SNP, .keep_all = TRUE)
          taxa_all_list <- NULL
          if (!is.null(metagenome_data)){
            traits <- metag
            traits <- subset(traits, select=-c(genus,family,order,class,phylum,kingdom,domain))
            if(grepl("_strain_",job)) {
              traits <- subset(traits, select=-c(tax_id,species))
            }
            rownames(traits) <- traits[,ncol(traits)]; traits <- traits[,-(ncol(traits))]
            colnames(traits) <-  sub("_mean.*", "", colnames(traits))
            traits$percent <- (rowSums(traits[,1:ncol(traits)] > "0")/ncol(traits))*100
            traits <- subset(traits, percent >= perc)
            traits <- subset(traits, select=-c(percent))
            traits <- data.frame(t(traits))
            taxa_all_list <- colnames(traits)
            for (i in 1:ncol(traits)){
              tiff(paste("./normal_distribution/",colnames(traits)[i],"_quantification.tiff",sep=""), width=5, height=5, units = 'in', res = 300, compression = 'lzw')
              plot <- hist(sqrt(traits[,i]+1), main=paste(colnames(traits)[i]," (n=",length(na.omit(traits[,i])),")"), xlab="log transformation",
                           border="grey", col="cornflowerblue")
              dev.off()
            }
            # traits <- sqrt(traits[]+1)
            traits$Plant_ID <- rownames(traits)
            traits <- subset(traits, select=c(ncol(traits),1:(ncol(traits)-1)))
            traits[, 2:ncol(traits)] <- lapply(traits[, 2:ncol(traits), drop=FALSE], function(x) .sgwas_numeric(x, "microbial abundance"))
            samples <- as.data.frame(colnames(geno[,-c(1:5)]))
            colnames(samples) <- "Plant_ID"
            traits <- merge(traits, samples, by="Plant_ID")
          }
          if (!is.null(metagenome_data)){traits_metag <- traits}

          if (!is.null(phenotype_data)){
            traits_pheno <- .sgwas_phenotypes(input_path(phenotype_data))
            for(nm in names(traits_pheno)[-1L])
              traits_pheno[[nm]] <- .sgwas_numeric(traits_pheno[[nm]],nm)
            all_traits_pheno <- colnames(traits_pheno)[2:ncol(traits_pheno)]
            for (i in 2:ncol(traits_pheno)){
              tiff(paste("./normal_distribution/",colnames(traits_pheno)[i],".tiff",sep=""), width=5, height=5, units = 'in', res = 300, compression = 'lzw')
              plot <- hist(traits_pheno[,i], main=paste(colnames(traits_pheno)[i]," (n=",length(na.omit(traits_pheno[,i])),")"),
                           xlab="Phenotypic Values", border="grey", col="cornflowerblue")
              dev.off()
            }
            if (!is.null(metagenome_data)){ traits <- merge(traits_pheno, traits, by="Plant_ID") }
            if (is.null(metagenome_data)){
              samples <- as.data.frame(colnames(geno[,-c(1:5)]))
              colnames(samples) <- "Plant_ID"
              traits <- merge(traits_pheno, samples, by="Plant_ID")
            }
          }
          names(traits) <- sub("\\.x$", "", names(traits))
          names(traits) <- sub("\\.y$", "", names(traits))
          # Aggregate duplicate IDs per column without dropping samples merely
          # because another trait/covariate is missing.
          traits_hold <- data.frame(Plant_ID=unique(as.character(traits$Plant_ID)))
          for(nm in names(traits)[-1L])
            traits_hold[[nm]] <- .sgwas_id_values(traits,nm,traits_hold$Plant_ID)
          traits <- traits_hold
          names(traits) <- sub("\\.x$", "", names(traits))
          names(traits) <- sub("\\.y$", "", names(traits))
          traits <- traits[, !duplicated(colnames(traits))]
          write.csv(traits, "traits.csv", row.names=F, quote = FALSE)


          #############################################################################################################################################################################
          # Perform Genome-Wide Association Analysis
          ##########################################
          traits <- read.csv("traits.csv", header=T, sep=",", check.names=FALSE, stringsAsFactors=FALSE, colClasses=c(Plant_ID="character"))

          if (!is.null(taxa_level)){
            if(taxa_level == "notaxa") {taxa <- c(trait_names)}
            if(taxa_level == "strain") {taxa <- c(taxa_strain)}
            if(taxa_level == "species") {taxa <- c(taxa_species)}
          }
          taxalist <- match(taxa,names(traits)); taxalist <- taxalist[!is.na(taxalist)]

          for (j in c(taxalist)) {
            i <- names(traits)[j]
            pheno <- subset(traits, select=c(1,j))
            pheno <- na.omit(pheno)
            traitname <- (colnames(pheno))[2]
            if(!is.null(alternate_trait)){ train_traitname <- alternate_trait} else {train_traitname <- traitname}
            if (!is.null(trait_microbial_proxy)){
              pheno_original <- pheno
              taxa_prefix <- NULL
              overlapping_taxa <- NULL
              if (length(trait_microbial_proxy) == 1 && grepl("(_spp| spp)$", trait_microbial_proxy)){
                taxa_prefix <- sub("(_spp| spp)$", "", trait_microbial_proxy)
                trait_microbial_proxy <- metag$species[grepl(paste0("^", taxa_prefix, "(_| )"), metag$species)]
              }
              perc <- as.numeric(perc)
              if(is.null(metag)){stop("Please provide provide metagenomic data")}
              metag_proxy <- subset(metag, select=-c(genus,family,order,class,phylum,kingdom,domain))
              metag_proxy <- subset(metag_proxy, species != "Severe_acute_respiratory_syndrome_related_coronavirus")
              metag_proxy <- subset(metag_proxy, select=c(ncol(metag_proxy),1:(ncol(metag_proxy)-1)))
              metag_proxy <- setNames(data.frame(t(metag_proxy[,-1])), metag_proxy[,1])
              metag_proxy <- data.frame(t(metag_proxy))
              metag_proxy[is.na(metag_proxy)] <- "0"
              for (k in 2:ncol(metag_proxy)){
                metag_proxy[,k] <- as.numeric(as.character(metag_proxy[,k]))
              }
              metag_proxy$percent <- (rowSums(metag_proxy > 0, na.rm = TRUE)/ncol(metag_proxy))*100
              metag_proxy <- subset(metag_proxy, percent >= perc)
              metag_proxy <- subset(metag_proxy, select=-c(percent))
              metag_proxy <- data.frame(t(metag_proxy))
              metag_proxy$percent <- (rowSums(metag_proxy > 0, na.rm = TRUE) / ncol(metag_proxy)) * 100
              metag_proxy <- subset(metag_proxy, percent >= perc)
              metag_proxy <- subset(metag_proxy, select=-c(percent))
              rownames(metag_proxy) <- sub("^X", "", rownames(metag_proxy))
              metag_proxy <- clr(metag_proxy + 1e-6)  # pseudocount to avoid log(0)
              metag_proxy <- as.data.frame(metag_proxy)
              row.names(pheno) <- pheno[,1]
              pheno <- subset(pheno, select=-(Plant_ID))
              metag_proxy <- merge(pheno, metag_proxy, by = 'row.names')
              rownames(metag_proxy) <- metag_proxy[,1]; metag_proxy <- metag_proxy[,-1]
              for (k in 1:ncol(metag_proxy)){
                metag_proxy[,k] <- as.numeric(as.character(metag_proxy[,k]))
              }

              assoc.proxy <- "NA"

              if(trait_microbial_proxy[1] != "auto" && trait_microbial_proxy[1] != "auto-null"){
                missing_taxa <- setdiff(trait_microbial_proxy, colnames(metag_proxy))
                if(length(missing_taxa) > 0){
                  warning(paste("The following taxa were not found:", paste(missing_taxa, collapse = ", ")))
                }
                overlapping_taxa <- intersect(trait_microbial_proxy, colnames(metag_proxy))
                if(length(overlapping_taxa) == 0){
                  stop("None of the pre-listed taxa were found in metag_proxy")
                }
                if(length(overlapping_taxa) == 1){
                  pheno <- metag_proxy[, overlapping_taxa, drop = FALSE]
                  selected_taxa_weight <- as.data.frame(cbind(proxy_trait = overlapping_taxa, value.var = 1))
                }
              }

              if(trait_microbial_proxy[1] == "auto" || length(overlapping_taxa) > 1){
                set.seed(123)
                if(trait_microbial_proxy[1] == "auto"){ X <- metag_proxy[, colnames(metag_proxy) != train_traitname, drop = FALSE]}
                if(length(overlapping_taxa) > 1){X <- metag_proxy[, overlapping_taxa, drop = FALSE]}
                Y <- metag_proxy[, train_traitname, drop = FALSE]
                Y_vec <- as.numeric(as.matrix(Y)[,1])
                n_boot <- 100
                # This bootstrap screen requires at least five candidate predictors.
                # Single-taxon proxies are handled directly in the branch above.
                keepX_candidates <- c(5L,10L,20L,30L,50L,100L,200L)
                keepX_candidates <- keepX_candidates[keepX_candidates <= ncol(X)]
                if(length(keepX_candidates) == 0L){
                  stop(paste0(train_traitname, "\tFAILED: bootstrap sPLS requires at least 5 candidate species/predictors after filtering; found ",
                              ncol(X), "."))
                }
                boot_selected <- vector("list", n_boot)
                for(b in seq_len(n_boot)){
                  idx <- sample(seq_len(nrow(X)), replace = TRUE)
                  Xb <- X[idx, , drop = FALSE]
                  Yb <- Y[idx, , drop = FALSE]
                  rownames(Xb) <- NULL; rownames(Yb) <- NULL
                  # Sample an index: sample(5L, 1) would draw from 1:5,
                  # rather than choose the sole allowed keepX value of five.
                  fit <- try(mixOmics::spls(X = Xb, Y = Yb, ncomp = 1, keepX = keepX_candidates[sample.int(length(keepX_candidates), 1L)]), silent = FALSE)
                  if(inherits(fit, "try-error")) next
                  tmp <- try(selectVar(fit, comp = 1)$X$value, silent = TRUE)
                  if(inherits(tmp, "try-error") || is.null(tmp) || nrow(tmp) == 0) next
                  boot_selected[[b]] <- rownames(tmp)
                }
                boot_selected <- Filter(Negate(is.null), boot_selected)
                if(length(boot_selected) == 0){stop(paste0(train_traitname, "\tFAILED: no taxa selected in bootstrap\n"))}
                taxa_freq <- sort(table(unlist(boot_selected)), decreasing = TRUE)
                taxa_freq_df <- data.frame(taxa = names(taxa_freq), frequency = as.numeric(taxa_freq))
                taxa_freq_df$prop <- taxa_freq_df$frequency / n_boot
                select_proxy_taxa <- taxa_freq_df$taxa[taxa_freq_df$prop >= 0.6]
                if(length(select_proxy_taxa) == 0){stop(paste0(train_traitname, "\tFAILED: no stable taxa selected\n"))}
                X_final <- X[, select_proxy_taxa, drop = FALSE]
                spls_model <- mixOmics::spls(X = X_final, Y = Y, ncomp = 1, keepX = ncol(X_final))
                proxy_trait <- spls_model$variates$X[,1]
                r <- cor(proxy_trait, Y_vec, method = "spearman", use = "complete.obs")
                if(!is.na(r) && r < 0){proxy_trait <- -proxy_trait
                  spls_model$variates$X[,1] <- -spls_model$variates$X[,1]
                  spls_model$loadings$X[,1] <- -spls_model$loadings$X[,1]
                }
                pheno <- data.frame(proxy_trait = proxy_trait)
                tmp <- selectVar(spls_model, comp = 1)$X$value
                selected_taxa_weight <- data.frame(proxy_trait = rownames(tmp), weight = tmp[,1], row.names = NULL)
                select_proxy_taxa <- selected_taxa_weight$proxy_trait
              }


              if(trait_microbial_proxy[1] == "auto-null"){
                set.seed(234)
                X <- as.matrix(metag_proxy[, colnames(metag_proxy) != train_traitname , drop = FALSE])
                pcs <- prcomp(X, center = TRUE, scale. = TRUE)$x[, 1:3, drop = FALSE]
                R <- matrix(rnorm(9), 3, 3)
                R <- qr.Q(qr(R))
                pcs_rot <- pcs %*% R
                pheno <- data.frame(proxy_null_trait = pcs_rot[, 1])
                selected_taxa_weight <- as.data.frame(cbind(proxy_trait = "full_microbiome", value.var = 0))
              }


              shared_samples <- intersect(traits[,1], row.names(pheno))
              pheno <- pheno[row.names(pheno) %in% shared_samples, , drop = FALSE]
              pheno <- as.data.frame(pheno)
              pheno <- pheno[, colnames(pheno) != train_traitname , drop = FALSE]

              comp1_perc <- "NA"
              if(ncol(pheno) == 1){
                pheno <- data.frame(Plant_ID = rownames(pheno), pheno, check.names = FALSE)
                rownames(pheno) <- NULL
                if (trait_microbial_proxy[1] == "auto"){
                  traitname <- paste0(train_traitname ,"_sPLS_proxy_trait")
                }
                if (trait_microbial_proxy[1] == "auto-null"){
                  traitname <- paste0(train_traitname ,"_proxy_null_trait")
                }
                if(trait_microbial_proxy[1] != "auto" && trait_microbial_proxy[1] != "auto-null"){
                  if(is.null(taxa_prefix)){
                    traitname <- paste0(train_traitname ,"_single_proxy_",(colnames(pheno))[2])
                  } else {
                    traitname <- paste0(train_traitname ,"_",taxa_prefix,"_fixed_feature_sPLS_proxy_trait")
                  }
                  colnames(pheno) <- c("Plant_ID", traitname)
                }
              }

              pheno_compPC1 <- merge(pheno_original, pheno, by = "Plant_ID")
              cor_trait_comp1 <- cor(pheno_compPC1[[2]], pheno_compPC1[[3]], method = "spearman", use = "complete.obs")
              if(is.null(taxa_prefix)){
                out_file <- paste0(train_traitname , "_",traitname, "_corr_", round(cor_trait_comp1, 3), ".txt")
              } else {
                out_file <- paste0(train_traitname, "_",taxa_prefix,"_proxy_",traitname, "_corr_", round(cor_trait_comp1, 3), ".txt")
              }
              write.table(selected_taxa_weight, file = out_file, sep = "\t", quote = FALSE, row.names = TRUE)
            }


            # Compute correlation coefficients or extract it from pre-computed correlations.
            if (trait_associated_covariates == TRUE &&  covariate_metag == TRUE) {
              if (grepl(".txt", corr_coeff)) {
                tcorr <- read.table(input_path(corr_coeff), header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
                tcorr <- subset(tcorr, tcorr[,1]==i)
                tcorr <- subset(tcorr, abs(tcorr[,3]) <= maxcorr  & abs(tcorr[,3]) >= mincorr & tcorr[,4] <= pvalue)
                rownames(tcorr) <- tcorr[,2]; tcorr <- subset(tcorr, select=c(3,4)); colnames(tcorr)[1:2] <- c("tcorr","tpmat")
                taxa_intersect <- intersect(rownames(tcorr), colnames(traits))
                tcorr <- tcorr[rownames(tcorr) %in% taxa_intersect, ]
              }
              if (corr_coeff == "full"){
                tdata <- traits
                tdata <- tdata[tdata$Plant_ID %in% pheno$Plant_ID, ]
                tdata <- (tdata[,2:ncol(tdata)])
                tdata_0 <- subset(tdata, select=c(i)); tdata <- tdata[,!(names(tdata) %in% i)]
                tdata <- as.data.frame(t(tdata))
                tdata$count <- rowSums(tdata[,1:ncol(tdata)] > 1, na.rm=TRUE)
                tdata$percent <- (tdata$count / (ncol(tdata)-1))*100
                tdata <- subset(tdata, tdata$count >= 25); tdata <- subset(tdata, tdata$percent >= perc)
                tdata <- subset(tdata, select=-c(count,percent))
                tdata <- as.data.frame(t(tdata)); tdata <- cbind(tdata_0, tdata); tdata_0 <- NULL
                tcorr <- cor(tdata[sapply(tdata, is.numeric)], method="spearman", use="pairwise.complete.obs")
                tpmat <- cor_pmat(tdata[sapply(tdata, is.numeric)], method = "spearman", exact=TRUE)
                tcorr <- subset(tcorr[,i], abs(tcorr[,i]) <= maxcorr  & abs(tcorr[,i]) >= mincorr )
                tcorr <- as.data.frame(tcorr)
                tpmat <- subset(tpmat[,i], tpmat[,i] <= 0.05)
                tpmat <- as.data.frame(tpmat)
                tcorr <- merge (tcorr,tpmat,by="row.names")
                rownames(tcorr) <- tcorr[,1]; tcorr <- subset(tcorr, select=c(2,3));
                taxa_intersect <- intersect(rownames(tcorr), colnames(traits))
                tcorr <- tcorr[rownames(tcorr) %in% taxa_intersect, ]
              }
              if (nrow(tcorr) >= 1) {
                factors <- c(rownames(tcorr))
                factors <- setdiff(factors,names(traits)[c(taxalist)])
                rtcorr <- rownames(tcorr)
                if(length(factors) > 0){
                  formula <- paste(names(traits)[j]," ~ ", paste(factors, collapse=" + ",sep=""))
                  ntraits <- traits
                  fit <- lm(paste(formula), data=ntraits)
                  coeff <- as.data.frame(summary(fit)$coefficients)
                  rsq <- as.data.frame(c(as.data.frame(summary(fit)$r.squared),"na","na","na")); row.names(rsq)[1] <- "r.sq"
                  adjrsq <- as.data.frame(c(as.data.frame(summary(fit)$adj.r.squared),"na","na","na")); row.names(adjrsq)[1] <- "adj.r.sq"
                  names(rsq) <- c("Estimate", "Std. Error", "t value",  "Pr(>|t|)")
                  names(adjrsq) <- c("Estimate", "Std. Error", "t value",  "Pr(>|t|)")
                  summout <- rbind(coeff,rsq); summout <- rbind(summout,adjrsq)
                  write.table(summout,paste("./correlations/",i,"_fit_multiple_regression_fit.txt",sep=""), row.names=T, quote = FALSE, sep="\t")
                }
                write.table(tcorr,paste("./correlations/",i,".txt",sep=""), row.names=T, quote = FALSE, sep="\t")

                rtcorr <- row.names(tcorr)
                plotcor <- tdata[,c(rtcorr,i)]
                plotcor <- plotcor %>% mutate_if(is.numeric, round, digits = 12)
                plotcor <- data.frame(apply(plotcor, 2, function(x) as.numeric(as.character(x))))
                cormat <- round(cor(plotcor, method="spearman", use="pairwise.complete.obs"), 2)
                corp <- cor.mtest(plotcor, conf.level = .95)
                corpvalue <- corp$p; colnames(corpvalue) <- colnames(cormat); rownames(corpvalue) <- colnames(cormat)
                tiff(paste("./correlations/Correlogram_",i,".tiff",sep=""), width=15, height=15, units = 'in', res = 300, compression = 'lzw')
                corrplot(cormat, method="circle", type = "lower", p.mat = corp$p, insig = "label_sig",
                         sig.level = c(.001, .01, .05), pch.cex = 1.0, pch.col = "white", order = "hclust", tl.col = "black", tl.srt = 45)
                dev.off()
              } else {tcorr <- NULL}
            } else { tcorr <- NULL }
            if (trait_associated_covariates == TRUE && length(tcorr) >= 1) {
              covariates <- traits[traits$Plant_ID %in% pheno$Plant_ID, ]
              selcov1 <- setdiff(names(covariates),names(traits_pheno))
              selcov2 <- intersect(names(covariates),names(pheno))
              covariates <- subset(covariates, select=c(selcov2,selcov1))
              row.names(covariates) <- covariates[,1]; covariates <- subset(covariates, select=-c(1,2))
              if (ncol(covariates) > 0){
                covariates <- prcomp(covariates, center = TRUE,scale. = TRUE)
                var_explained <- (covariates$sdev^2/sum(covariates$sdev^2))*100
                var_explained[1:5]
                covariates <- as.data.frame(covariates$x)
                if (min(covariates[!is.na(covariates)])-1 < 0){
                  for (c in c(1:ncol(covariates))) {
                    covariates[,c] <- covariates[,c] + abs(min(covariates[!is.na(covariates)])-1)
                  }
                }
                for (c in c(1:ncol(covariates))) {
                  covariates[,c] <- covariates[,c]/max(covariates[], na.rm=TRUE)
                }
                if (ncol(covariates) > 2) {
                  covariates <- subset(covariates, select=c(1:3))
                  covariates$PC4 <- covariates$PC1 + covariates$PC2 + covariates$PC3
                  covariates[] <- covariates[]/(max(covariates$PC4)+(max(covariates$PC4)*0.1))
                  covariates <- covariates[,-c(4)]
                  covariates$Taxa <- row.names(covariates)
                  covariates[is.na(covariates)] <- 0
                  covariates <- covariates[,c(4,1,2,3)]
                } else {
                  if (ncol(covariates) == 2) {
                    covariates$PC3 <- covariates$PC1 + covariates$PC2
                    covariates[] <- covariates[]/(max(covariates$PC3)+(max(covariates$PC3)*0.1))
                    covariates <- covariates[,-c(3)]
                    covariates$Taxa <- row.names(covariates)
                    covariates <- covariates[,c(3,1,2)]
                  }
                  if (ncol(covariates) == 1) {
                    covariates$PC2 <- covariates$PC1
                    covariates[] <- covariates[]/(max(covariates$PC2)+(max(covariates$PC2)*0.1))
                    covariates <- subset(covariates, select=-c(2))
                    covariates$Taxa <- row.names(covariates)
                    covariates <- covariates[,c(2,1)]
                  }
                }
              }
              covariates$Plant_ID <- row.names(covariates)
              pheno <- merge(pheno, covariates, by=c("Plant_ID"))
              pheno <- subset(pheno, select=-c(Taxa))
            }

            if (!is.null(covariatename)) {
              covariates <- traits[traits$Plant_ID %in% pheno$Plant_ID, ]
              covariates <- subset(covariates, select=c("Plant_ID",i,covariatename))
              imp <- mice(covariates, m = 5, method = "pmm", seed = 123)
              covariates <- complete(imp, 1)   # first imputed dataset
              covariates <- subset(covariates, select=c("Plant_ID",covariatename))
              pheno <- merge(pheno,covariates, by=c("Plant_ID"))
            }
            if (covariate == TRUE){pheno <- pheno[,-c(3:5)]}

            # Run model with metagenome and generate residual
            if(metag_method == "Aitchison" && !is.null(metag) && covariate_metag == TRUE) {
              covariates <- traits[traits$Plant_ID %in% pheno$Plant_ID, ]
              selcov1 <- setdiff(names(covariates),names(traits_pheno))
              selcov2 <- intersect(names(covariates),names(pheno))
              covariates <- subset(covariates, select=c(selcov2,selcov1))
              row.names(covariates) <- covariates[,1]; covariates <- subset(covariates, select=-c(1,2))
              metag_clr <- as.matrix(clr(covariates + 1e-6))
              metagKI <- tcrossprod(scale(metag_clr)) / ncol(metag_clr)
              metagKI <- normalize_kinmat(as.matrix(metagKI))
              phenor <- pheno[,(1:2)]
              phenor$ID_M <- phenor$Plant_ID
              fixed_formula <- as.formula(paste(i, "~ 1"))
              model_metag <- mmer(
                fixed = fixed_formula,
                random = ~ vsr(ID_M, Gu = metagKI),
                rcov = ~ units,
                data = phenor
              )
              summary(model_metag)
              pheno_resid <- resid(model_metag)
              pheno_gwas <- data.frame(Plant_ID = pheno$Plant_ID, resid_trait = pheno_resid[,ncol(pheno_resid)])
              colnames(pheno_gwas)[2] <- i
              pheno <- as.data.frame(pheno[,-2]); colnames(pheno)[1] <- "Plant_ID"
              pheno <- merge(pheno_gwas, pheno, by="Plant_ID")
            }
            write.csv(pheno,'pheno.csv', row.names=F, quote = FALSE)
            if(is.null(alternate_trait)){
              if(!is.null(trait_microbial_proxy)){ write.csv(pheno, paste0(names(traits)[j],"_proxy_pheno.csv"), row.names=F, quote = FALSE)}
            } else {
              if(!is.null(trait_microbial_proxy)){ write.csv(pheno, paste0(names(traits)[j],"_",alternate_trait,"_proxy_pheno.csv"), row.names=F, quote = FALSE)}
            }

            geno <- read.table(input_path(genotype_data), header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
            geno <- geno %>%
              distinct(SNP, .keep_all = TRUE)
            samples <- lapply(pheno[,1], as.character)
            geno_hold <- subset(geno, select=c(1:5))
            for (i in 1:length(samples)){
              sampled <- subset(geno, select=c(samples[[i]]))
              geno_hold <- cbind(geno_hold, sampled)
            }
            geno <- geno_hold
            if (ploidy == 2){
              dG <- geno
              dG$freq0 <- (rowSums(dG == "0", na.rm = TRUE))*2 + (rowSums(dG == "1", na.rm = TRUE))*1
              dG$freq1 <- (rowSums(dG == "1", na.rm = TRUE))*1 + (rowSums(dG == "2", na.rm = TRUE))*2
              maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
              dG$min <- apply(dG[,(ncol(dG)-1):ncol(dG)], 1, function(x)x[maxn(2)(x)])
              dG$sum <- rowSums(dG[,c("freq0","freq1")], na.rm=TRUE)
              dG$maf <- as.numeric(dG$min)/as.numeric(dG$sum)
              dG <- dG[dG$maf >= maf, ]
              geno <- dG[,-c((ncol(dG)-4):ncol(dG))]
            }
            if (ploidy == 4){
              dG <- geno
              dG$freq0 <- (rowSums(dG == "0", na.rm = TRUE))*4 + (rowSums(dG == "1", na.rm = TRUE))*3 +
                (rowSums(dG == "2", na.rm = TRUE))*2 + (rowSums(dG == "3", na.rm = TRUE))*1
              dG$freq1 <- (rowSums(dG == "1", na.rm = TRUE))*1 + (rowSums(dG == "2", na.rm = TRUE))*2 +
                (rowSums(dG == "3", na.rm = TRUE))*3 + (rowSums(dG == "4", na.rm = TRUE))*4
              maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
              dG$min <- apply(dG[,(ncol(dG)-1):ncol(dG)], 1, function(x)x[maxn(2)(x)])
              dG$sum <- rowSums(dG[,c("freq0","freq1")], na.rm=TRUE)
              dG$maf <- as.numeric(dG$min)/as.numeric(dG$sum)
              dG <- dG[dG$maf >= maf, ]
              geno <- dG[,-c((ncol(dG)-4):ncol(dG))]
            }
            if (ploidy == 6){
              dG <- geno
              dG$freq0 <- (rowSums(dG == "0", na.rm = TRUE))*6 + (rowSums(dG == "1", na.rm = TRUE))*5 +
                (rowSums(dG == "2", na.rm = TRUE))*4 + (rowSums(dG == "3", na.rm = TRUE))*3 +
                (rowSums(dG == "4", na.rm = TRUE))*2 + (rowSums(dG == "5", na.rm = TRUE))*1
              dG$freq1 <- (rowSums(dG == "1", na.rm = TRUE))*1 + (rowSums(dG == "2", na.rm = TRUE))*2 +
                (rowSums(dG == "3", na.rm = TRUE))*3 + (rowSums(dG == "4", na.rm = TRUE))*4 +
                (rowSums(dG == "5", na.rm = TRUE))*5 + (rowSums(dG == "6", na.rm = TRUE))*6
              maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
              dG$min <- apply(dG[,(ncol(dG)-1):ncol(dG)], 1, function(x)x[maxn(2)(x)])
              dG$sum <- rowSums(dG[,c("freq0","freq1")], na.rm=TRUE)
              dG$maf <- as.numeric(dG$min)/as.numeric(dG$sum)
              dG <- dG[dG$maf >= maf, ]
              geno <- dG[,-c((ncol(dG)-4):ncol(dG))]
            }
            if (ploidy == 8){
              dG <- geno
              dG$freq0 <- (rowSums(dG == "0", na.rm = TRUE))*8 + (rowSums(dG == "1", na.rm = TRUE))*7 +
                (rowSums(dG == "2", na.rm = TRUE))*6 + (rowSums(dG == "3", na.rm = TRUE))*5 +
                (rowSums(dG == "4", na.rm = TRUE))*4 + (rowSums(dG == "5", na.rm = TRUE))*3 +
                (rowSums(dG == "6", na.rm = TRUE))*2 + (rowSums(dG == "7", na.rm = TRUE))*1
              dG$freq1 <- (rowSums(dG == "1", na.rm = TRUE))*1 + (rowSums(dG == "2", na.rm = TRUE))*2 +
                (rowSums(dG == "3", na.rm = TRUE))*3 + (rowSums(dG == "4", na.rm = TRUE))*4 +
                (rowSums(dG == "5", na.rm = TRUE))*5 + (rowSums(dG == "6", na.rm = TRUE))*6 +
                (rowSums(dG == "7", na.rm = TRUE))*7 + (rowSums(dG == "8", na.rm = TRUE))*8
              maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
              dG$min <- apply(dG[,(ncol(dG)-1):ncol(dG)], 1, function(x)x[maxn(2)(x)])
              dG$sum <- rowSums(dG[,c("freq0","freq1")], na.rm=TRUE)
              dG$maf <- as.numeric(dG$min)/as.numeric(dG$sum)
              dG <- dG[dG$maf >= maf, ]
              geno <- dG[,-c((ncol(dG)-4):ncol(dG))]

            }

            pop_struc <- function() {
              pop_data <- subset(geno, select=-c(1:5))
              pop_data$no_missing <- apply(pop_data, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )
              pop_data <- subset(pop_data, no_missing < ncol(pop_data)*0.3)
              pop_data <- subset(pop_data, select=-c(no_missing))
              if (nrow(pop_data) >= 1000) {
                pop_data <- as.matrix(t(pop_data))
                #Computing the full-autopolyploid matrix based on Slater 2016 (Eq. 8 and 9)
                G_matrix <- Gmatrix(SNPmatrix = pop_data, method = Gmethod, missingValue = NA,
                                    maf = maf, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                    pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                    ratio.check = FALSE)
                G_matrix <- normalize_kinmat(as.matrix(G_matrix))
                write.table(G_matrix, file=paste(pop,"_",Gmethod,"_",ploidy,"x.txt",sep=""), row.names=T, col.names = T, quote = FALSE, sep = "\t")
              } else {
                print ("Not enough markers to compute Gmatrix (i.e. threshold of 1000 markers)")
              }
            }
            pop_struc()
            #### Only strip subgenome prefixes for SNP that are not anchored to a single reference genome.
            if (grepl("_Chr", geno$SNP[1])){
              geno$CHROM <- sub(".*\\_","",geno$CHROM)
              geno$SNP <- sub(".*?_","",geno$SNP)
            }
            names(geno)[names(geno) == "SNP"] <- "Markers"
            names(geno)[names(geno) == "CHROM"] <- "Chrom"
            names(geno)[names(geno) == "POS"] <- "Position"
            geno <- geno[!grepl("Chr00", geno$Chrom),]
            geno <- subset(geno, select = -c(4,5))
            if("pvalue" %in% colnames(geno)) {
              geno <- subset(geno, select = -c(pvalue))
            }
            geno <- geno[order(geno$Chrom, geno$Position), ]
            write.csv(geno,'geno.csv', row.names=F, quote = FALSE)
            # geno <- read.csv("geno.csv", header=T, sep=",", fill=TRUE, check.names=FALSE)
            G_matrix <- read.table(file=paste(pop,"_",Gmethod,"_",ploidy,"x.txt",sep=""), header=T, sep="\t", check.names=FALSE)
            if (tolower(gwas_method) == "glm") {G_matrix[G_matrix >= 0] <- 1}

            data <- read.GWASpoly(ploidy=ploidy, pheno.file="pheno.csv", geno.file="geno.csv", format="numeric", n.traits=1, delim=",")
            Kinship <- set.K(data, K=as.matrix(G_matrix), LOCO=FALSE)
            if(LOCO == TRUE){Kinship <- set.K(data, LOCO = TRUE, n.core=cores)}
            if (length(tcorr) == 0) {
              if (is.null(covariatename) ) {
                params <- set.params(MAF=maf,P3D=TRUE)
                if ( ploidy == "2" ) {
                  if (model == "Add"){
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive"), n.core = cores, quiet = FALSE, params=params)
                    models <- c("additive")
                  }
                  if (model == "Add_Dom"){
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive","1-dom"), n.core = cores, quiet = FALSE, params=params)
                    models <- c("additive","1-dom-ref","1-dom-alt")
                  }
                }
                if ( ploidy == "4" ) {
                  GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom","additive"), n.core = cores, quiet = FALSE, params=params)
                  models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","additive")
                }
                if ( ploidy == "6" ) {
                  GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom","3-dom","additive"), n.core = cores, quiet = FALSE, params=params)
                  models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt","additive")
                }
                if ( ploidy == "8" ) {
                  GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom", "3-dom", "4-dom","additive"),
                                          traits=c(traitname), n.core = cores, quiet = FALSE, params=params)
                  models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt","4-dom-ref","4-dom-alt","additive")
                }
              }
              if (!is.null(covariatename) ){
                params <- set.params(fixed=c(colnames(pheno)[3:ncol(pheno)]), fixed.type=rep("numeric",(ncol(pheno)-2)),n.PC=(ncol(pheno)-2),MAF=maf,geno.freq=0.99,P3D=TRUE)
                if ( ploidy == "2" ) {
                  if (model == "Add"){
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive"), n.core = cores, quiet = FALSE, params=params)
                    models <- c("additive")
                  }
                  if (model == "Add_Dom"){
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive","1-dom"), n.core = cores, quiet = FALSE, params=params)
                    models <- c("additive","1-dom-ref","1-dom-alt")
                  }
                }
                if ( ploidy == "4" ) {
                  GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom","additive"), n.core = cores, quiet = FALSE, params=params)
                  models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","additive")
                }
                if ( ploidy == "6" ) {
                  GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom","3-dom","additive"), n.core = cores, quiet = FALSE, params=params)
                  models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt","additive")
                }
                if ( ploidy == "8" ) {
                  GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom", "3-dom", "4-dom","additive"),
                                          traits=c(traitname), n.core = cores, quiet = FALSE, params=params)
                  models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt","4-dom-ref","4-dom-alt","additive")
                }
              }
            } else {
              if(ncol(pheno) == 2){
                params <- set.params(MAF=maf,geno.freq=0.99,P3D=TRUE)
              }
              if(ncol(pheno) == 3){
                params <- set.params(fixed=c(colnames(pheno)[3]), fixed.type=rep("numeric",1),n.PC=1,MAF=maf,geno.freq=0.99,P3D=TRUE)
              }
              if(ncol(pheno) == 4){
                params <- set.params(fixed=c(colnames(pheno)[3:4]), fixed.type=rep("numeric",2),n.PC=2,MAF=maf,geno.freq=0.99,P3D=TRUE)
              }
              if(ncol(pheno) > 4){
                params <- set.params(fixed=c(colnames(pheno)[3:5]), fixed.type=rep("numeric",3),n.PC=3,MAF=maf,geno.freq=0.99,P3D=TRUE)
              }

              if ( ploidy == "2" ) {
                if (model == "Add"){
                  GWAS.fitted <- GWASpoly(Kinship, models = c("additive"), n.core = cores, quiet = FALSE, params=params)
                  models <- c("additive")
                }
                if (model == "Add_Dom"){
                  GWAS.fitted <- GWASpoly(Kinship, models = c("additive","1-dom"), n.core = cores, quiet = FALSE, params=params)
                  models <- c("additive","1-dom-ref","1-dom-alt")
                }
              }
              if ( ploidy == "4" ) {
                GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom","additive"), n.core = cores, quiet = FALSE, params=params)
                models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","additive")
              }
              if ( ploidy == "6" ) {
                GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom","3-dom","additive"), n.core = cores, quiet = FALSE, params=params)
                models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt","additive")
              }
              if ( ploidy == "8" ) {
                GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom", "3-dom", "4-dom","additive"),
                                        traits=c(traitname), n.core = cores, quiet = FALSE, params=params)
                models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt","4-dom-ref","4-dom-alt","additive")
              }
            }

            # qqplot(GWAS.fitted) + ggtitle(label="Q-Q plot")

            SNP <- GWAS.fitted@map[["Marker"]]
            GWAS_logP<- data.frame(
              SNP = SNP,
              scores  = GWAS.fitted@scores[[colnames(pheno[2])]]
            )
            GWAS_logP$no_missing <- apply(GWAS_logP, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )
            GWAS_logP <- subset(GWAS_logP, no_missing != ncol(GWAS_logP)-2)
            GWAS_logP <- subset(GWAS_logP, select=-c(no_missing))
            colnames(GWAS_logP) <- paste0(colnames(GWAS_logP), "_scores")
            rownames(GWAS_logP) <- GWAS_logP$SNP; GWAS_logP <- GWAS_logP[,-1]
            GWAS_logP <- cbind(SNP = rownames(GWAS_logP), GWAS_logP)
            colnames(GWAS_logP) <- gsub("^scores\\.", "", colnames(GWAS_logP))

            if ( !is.null(GWAS_logP) ) {
              GWAS_effects <- GWAS.fitted@effects[[colnames(pheno[2])]]
              GWAS_effects$no_missing <- apply(GWAS_effects, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )
              GWAS_effects <- subset(GWAS_effects, no_missing != ncol(GWAS_effects)-1)
              GWAS_effects <- subset(GWAS_effects, select=-c(no_missing))
              colnames(GWAS_effects) <- paste0(colnames(GWAS_effects), "_effects")
              GWAS_effects <- cbind(SNP = rownames(GWAS_effects), GWAS_effects)
              GWAS_scores_effects <- merge(GWAS_logP, GWAS_effects, by=c("SNP"))

              geno_mat <- t(as.matrix(geno[ , -c(2:3)]))
              colnames(geno_mat) <- geno_mat[1,]; geno_mat <- geno_mat[-1,]
              geno_mat <- apply(geno_mat, 2, as.numeric)
              geno_maf <- data.frame(
                SNP = colnames(geno_mat),
                MAF = apply(geno_mat, 2, function(x) {
                  x <- as.numeric(x)
                  x <- x[!is.na(x)]
                  if (length(x) == 0) return(NA)
                  p <- mean(x) / as.numeric(ploidy)   # allele frequency adjusted by ploidy
                  return(min(p, 1 - p))               # minor allele frequency
                })
              )
              GWAS_scores_effects <- merge(GWAS_scores_effects, geno_maf, by = "SNP")
              var_y <- var(pheno[,2], na.rm = TRUE)
              colnames(GWAS_scores_effects) <- gsub("^scores\\.", "", colnames(GWAS_scores_effects))
              # Calculate PVE using the variance of the predictor fitted by each model
              geno_mat_imputed <- as.matrix(geno_mat)
              if (is.null(colnames(geno_mat_imputed))) {
                stop("geno_mat must have SNP IDs as column names.")
              }
              if (!is.numeric(var_y) || length(var_y) != 1L ||
                  !is.finite(var_y) || var_y <= 0) {
                stop("var_y must be a single positive finite value.")
              }

              # GWASpoly normally imputes genotype data during import.
              # This fallback reproduces its general imputation behaviour:
              #   - fractional dosage data: marker mean
              #   - integer dosage data: marker mode
              impute_mean <- function(x) {
                if (all(is.na(x))) {
                  stop("A genotype column contains only missing values.")
                }

                x[is.na(x)] <- mean(x, na.rm = TRUE)
                x
              }
              impute_mode <- function(x) {
                if (all(is.na(x))) {
                  stop("A genotype column contains only missing values.")
                }

                observed <- x[!is.na(x)]
                mode_value <- as.numeric(names(which.max(table(observed))))
                x[is.na(x)] <- mode_value
                x
              }
              has_fractional_dosage <- any(
                abs(geno_mat_imputed - round(geno_mat_imputed)) > 1e-8,
                na.rm = TRUE
              )
              if (anyNA(geno_mat_imputed)) {
                if (has_fractional_dosage) {
                  geno_mat_imputed <- apply(
                    geno_mat_imputed,
                    2,
                    impute_mean
                  )
                } else {
                  geno_mat_imputed <- apply(
                    geno_mat_imputed,
                    2,
                    impute_mode
                  )
                }

                geno_mat_imputed <- as.matrix(geno_mat_imputed)
              }

              # Match GWAS results to genotype columns once
              snp_index <- match(
                GWAS_scores_effects$SNP,
                colnames(geno_mat_imputed)
              )
              if (anyNA(snp_index)) {
                warning(
                  sum(is.na(snp_index)),
                  " SNP(s) in GWAS_scores_effects were not found in geno_mat."
                )
              }
              calculate_pve <- function(effect, predictor_variance) {
                pve <- (effect^2 * predictor_variance) / var_y
                pve[!is.finite(pve)] <- NA_real_

                # Retained from the original code.
                # Consider removing pmin() during diagnostics so values >1 are visible.
                pmin(pve, 1)
              }

              # ------------------------------------------------------------------
              # Additive PVE: predictor is alternate-allele dosage
              # ------------------------------------------------------------------
              additive_variance <- apply(
                geno_mat_imputed,
                2,
                var
              )
              v_add <- additive_variance[snp_index]
              GWAS_scores_effects$additive_PVE <- calculate_pve(
                effect = GWAS_scores_effects$additive_effects,
                predictor_variance = v_add
              )

              # ------------------------------------------------------------------
              # Dominance PVE: predictor is the model-specific binary coding
              # ------------------------------------------------------------------
              geno_dom <- round(geno_mat_imputed)
              max_dom <- ploidy %/% 2L
              for (d in seq_len(max_dom)) {

                alt_effect_col <- paste0(d, "-dom-alt_effects")
                ref_effect_col <- paste0(d, "-dom-ref_effects")

                # GWASpoly d-dom-alt coding:
                # 0 when alternate dosage < d
                # 1 when alternate dosage >= d
                dom_alt <- 1L * (geno_dom >= d)

                # GWASpoly d-dom-ref coding:
                # 0 when alternate dosage <= ploidy - d
                # 1 when alternate dosage >  ploidy - d
                dom_ref <- 1L * (geno_dom > (ploidy - d))

                dom_alt_variance <- apply(dom_alt, 2, var)
                dom_ref_variance <- apply(dom_ref, 2, var)

                if (alt_effect_col %in% colnames(GWAS_scores_effects)) {
                  GWAS_scores_effects[[paste0(d, "-dom-alt_PVE")]] <-
                    calculate_pve(
                      effect = GWAS_scores_effects[[alt_effect_col]],
                      predictor_variance = dom_alt_variance[snp_index]
                    )
                }

                if (ref_effect_col %in% colnames(GWAS_scores_effects)) {
                  GWAS_scores_effects[[paste0(d, "-dom-ref_PVE")]] <-
                    calculate_pve(
                      effect = GWAS_scores_effects[[ref_effect_col]],
                      predictor_variance = dom_ref_variance[snp_index]
                    )
                }
              }


              write.table(GWAS_scores_effects, file=paste("./scores_effects/","score_effects_",traitname,".txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
              colnames(GWAS_logP) <- gsub("_scores", "", colnames(GWAS_logP)); GWAS_logP <- GWAS_logP[,-1]
              colnames(GWAS_effects) <- gsub("_effects", "", colnames(GWAS_effects)); GWAS_effects <- GWAS_effects[,-1]

              pvalues <- as.data.frame(10^(-1 * as.data.frame(GWAS_logP)))
              score_models <- colnames(pvalues)
              if(ncol(pvalues) ==1 ){colnames(pvalues) <- score_models}
              pvalues <- as.data.frame(reshape2::melt(pvalues)); colnames(pvalues) <- c("model","pvalue")
              pvalues <- na.omit(pvalues); pvalues$pvalue <- as.numeric(as.character(pvalues$pvalue))
              ps <- pvalues; ci <- 0.95
              ps <- ps[order(ps$model,ps$pvalue),]
              ps$observed <- -log10(ps$pvalue)
              ps$expected <- ps$model; ps$clower <- ps$model; ps$cupper <- ps$model
              ps$expected <- ps$clower <- ps$cupper <- rep(NA_real_, nrow(ps))

              for (jj in score_models) {
                count_df <- as.data.frame(table(ps$model))
                count_sub <- subset(count_df, Var1 == jj)  # use count_sub here

                if (nrow(count_sub) == 0) next  # skip if no rows

                count <- count_sub[, 2]  # get the number of points
                idx <- ps$model == jj

                ps$expected[idx] <- -log10(ppoints(count))

                p1 <- (1 - ci) / 2
                p2 <- (1 + ci) / 2
                shape1 <- 1:count
                shape2 <- count:1

                ps$clower[idx] <- -log10(qbeta(p1, shape1, shape2))
                ps$cupper[idx] <- -log10(qbeta(p2, shape1, shape2))
              }
              ps$expected <- as.numeric(as.character(ps$expected))
              ps$clower <- as.numeric(as.character(ps$clower))
              ps$cupper <- as.numeric(as.character(ps$cupper))
              log10Pe <- expression(paste("Expected -log"[10],plain(P),sep="")); log10Po <- expression(paste("Observed -log"[10], plain(P),sep=""))
              qqplot_metric <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("model","Trait","dev_norm_avGW"))
              for (jjj in c(score_models)) {
                psout <- subset(ps, ps$model == jjj)
                dev_norm_av <- round(mean(psout$observed - psout$expected), digits=5)
                qqplot_metric1 <- as.data.frame(t(c(jjj,colnames(pheno[2]),dev_norm_av))); names(qqplot_metric1) <- c("model","Trait","dev_norm_avGW")
                qqplot_metric <- rbind(qqplot_metric, qqplot_metric1)
              }
              write.table(qqplot_metric, file=paste("./qqplots/",traitname,"_qqplot_metric.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")


              if (ploidy == 2) {
                qqplot <- ggplot(ps, aes(x=expected, y=observed, group=model)) + facet_wrap(~model, ncol=1)+
                  geom_point(shape = 16, size = 2, colour="cornflowerblue") + scale_shape_identity()+
                  geom_abline(intercept = 0, slope = 1, alpha = 0.5, color="red", linetype="dashed", size=1.0)+
                  geom_ribbon(aes(x=expected, ymin=clower, ymax=cupper),alpha = 0.1) + theme_gray()+
                  theme(plot.title = element_text(hjust = 0.5), strip.background=element_rect(fill="grey"),
                        strip.text = element_text(size=10, color="black"))+
                  xlab(log10Pe) + ylab(log10Po) + labs(title=paste(traitname))
                ggsave(file=paste("./qqplots/","QQplot_",traitname,"_fdr0.05.tiff",sep=""), plot=qqplot, width=4.5, height=4, units=("in"), dpi=120, compression = "lzw")
              }
              if (ploidy == 4) {
                qqplot <- ggplot(ps, aes(x=expected, y=observed, group=model)) + facet_wrap(~model, ncol=1)+
                  geom_point(shape = 16, size = 2, colour="cornflowerblue") + scale_shape_identity()+
                  geom_abline(intercept = 0, slope = 1, alpha = 0.5, color="red", linetype="dashed", size=1.0)+
                  geom_ribbon(aes(x=expected, ymin=clower, ymax=cupper),alpha = 0.1) + theme_gray()+
                  theme(plot.title = element_text(hjust = 0.5), strip.background=element_rect(fill="grey"),
                        strip.text = element_text(size=10, color="black"))+
                  xlab(log10Pe) + ylab(log10Po) + labs(title=paste(traitname))
                ggsave(file=paste("./qqplots/","QQplot_",traitname,"_fdr0.05.tiff",sep=""), plot=qqplot, width=4.5, height=8, units=("in"), dpi=120, compression = "lzw")
              }
              if (ploidy == 6) {
                qqplot <- ggplot(ps, aes(x=expected, y=observed, group=model)) + facet_wrap(~model, ncol=1)+
                  geom_point(shape = 16, size = 2, colour="cornflowerblue") + scale_shape_identity()+
                  geom_abline(intercept = 0, slope = 1, alpha = 0.5, color="red", linetype="dashed", size=1.0)+
                  geom_ribbon(aes(x=expected, ymin=clower, ymax=cupper),alpha = 0.1) + theme_gray()+
                  theme(plot.title = element_text(hjust = 0.5), strip.background=element_rect(fill="grey"),
                        strip.text = element_text(size=10, color="black"))+
                  xlab(log10Pe) + ylab(log10Po) + labs(title=paste(traitname))
                ggsave(file=paste("./qqplots/","QQplot_",traitname,"_fdr0.05.tiff",sep=""), plot=qqplot, width=4.5, height=12, units=("in"), dpi=120, compression = "lzw")
              }
              if (ploidy == 8) {
                qqplot <- ggplot(ps, aes(x=expected, y=observed, group=model)) + facet_wrap(~model, ncol=1)+
                  geom_point(shape = 16, size = 2, colour="cornflowerblue") + scale_shape_identity()+
                  geom_abline(intercept = 0, slope = 1, alpha = 0.5, color="red", linetype="dashed", size=1.0)+
                  geom_ribbon(aes(x=expected, ymin=clower, ymax=cupper),alpha = 0.1) + theme_gray()+
                  theme(plot.title = element_text(hjust = 0.5), strip.background=element_rect(fill="grey"),
                        strip.text = element_text(size=10, color="black"))+
                  xlab(log10Pe) + ylab(log10Po) + labs(title=paste(traitname))
                ggsave(file=paste("./qqplots/","QQplot_",traitname,"_fdr0.05.tiff",sep=""), plot=qqplot, width=4.5, height=16, units=("in"), dpi=120, compression = "lzw")
              }

              # data_r2_fdr <- set.threshold(GWAS.fitted,method="FDR",level=0.05,n.core=cores)
              # SigQTL_r2_fdr <- fit.QTL(data_r2_fdr, trait=paste(traitname), qtl=qtl[,c("Marker","Model")])
              # if (is.null(SigQTL_r2_fdr) == "TRUE") {print ("file is empty")} else{
              #   if (is.null(SigQTL_r2_fdr) == "FALSE") {
              #     write.table(SigQTL_r2_fdr, paste("./sigFDR/","Significant_R2_",traitname,"_fdr0.05.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
              #   }}
              # data_r2_Bonferroni <- set.threshold(GWAS.fitted,method="Bonferroni",level=0.05,n.core=cores)
              # SigQT_r2L_Bonferroni <- fit.QTL(data_r2_Bonferroni, trait=paste(traitname), qtl=qtl[,c("Marker","Model")])
              # if (is.null(SigQTL_r2_Bonferroni) == "TRUE") {print ("file is empty")} else{
              #   if (is.null(SigQTL_r2_Bonferroni) == "FALSE") {
              #     write.table(SigQTL_r2_Bonferroni, paste("./sigBonferroni/","Significant_R2_",traitname,"_Bonferroni0.05.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
              #   }}
              # data_r2_Meff <- set.threshold(GWAS.fitted,method="M.eff",level=0.05,n.core=cores)
              # SigQT_r2L_Meff <- fit.QTL(data_r2_Meff, trait=paste(traitname), qtl=qtl[,c("Marker","Model")])
              # if (is.null(SigQTL_r2_Meff) == "TRUE") {print ("file is empty")} else{
              #   if (is.null(SigQTL_r2_Meff) == "FALSE") {
              #     write.table(SigQTL_r2_Meff, paste("./sigBonferroni/","Significant_R2_",traitname,"_Meff.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
              #   }}
              # if (is.null(permutations)) { permutations=0}
              # if (permutations >= 100) {
              #   data_r2_permute <- set.threshold(GWAS.fitted,method="permute",n.permute=permutations,level=0.05,n.core=cores)
              #   SigQTL_r2_permute <- fit.QTL(data_permute, trait=paste(traitname), qtl=qtl[,c("Marker","Model")])
              #   if (is.null(SigQTL_r2_permute) == "TRUE") {print ("file is empty")} else{
              #     if (is.null(SigQTL_r2_permute) == "FALSE") {
              #       write.table(SigQTL_r2_permute, paste("./sigpermute/","Significant_R2_",traitname,"_permute0.05.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
              #     }}
              # }

              # --- Determine dominance terms based on ploidy ---
              get_dom_terms <- function(ploidy) {
                # diploid (2): only 1 level of dominance
                max_dom <- ploidy / 2
                dom_terms <- c()
                for (i in 1:max_dom) {
                  dom_terms <- c(
                    dom_terms,
                    paste0(i, "-dom-alt_PVE"),
                    paste0(i, "-dom-ref_PVE")
                  )
                }
                return(dom_terms)
              }
              dom_terms <- get_dom_terms(ploidy)
              # filter only columns that exist in the dataset
              dom_terms <- dom_terms[dom_terms %in% colnames(GWAS_scores_effects)]
              # --- Build wide table ---
              GWAS_scores_effects_wide <- GWAS_scores_effects %>%
                dplyr::select(SNP,MAF,additive_PVE,dplyr::all_of(dom_terms))

              # --- Pivot longer ---
              GWAS_scores_effects_long <- GWAS_scores_effects_wide %>%
                pivot_longer(
                  cols = -c(SNP, MAF),
                  names_to = "Model",
                  values_to = "PVE"
                ) %>%
                mutate(
                  # clean model labels
                  Model = case_when(
                    Model == "additive_PVE" ~ "additive",
                    str_detect(Model, "dom-alt") ~ str_replace(Model, "_PVE", ""),
                    str_detect(Model, "dom-ref") ~ str_replace(Model, "_PVE", ""),
                    TRUE ~ Model
                  )
                )

              # GWAS_scores_effects_wide <- GWAS_scores_effects %>%
              #   select(SNP, MAF, additive_PVE, `1-dom-alt_PVE`, `1-dom-alt_PVE` = `1-dom-alt_PVE`)
              # GWAS_scores_effects_long <- GWAS_scores_effects_wide %>%
              #   pivot_longer(cols = c(additive_PVE, `1-dom-alt_PVE`, `1-dom-alt_PVE`),names_to = "Model",values_to = "PVE"
              #   ) %>%
              #   mutate(Model = case_when(Model %in% c("additive_PVE", "PVE_additive") ~ "additive",
              #                            Model %in% c("1-dom-ralt_PVE", "1-dom-alt_PVE") ~ "1-dom-alt",Model %in% c("1-dom-ref_PVE") ~ "1-dom-ref",
              #                            TRUE ~ Model)
              #   )
              colnames(GWAS_scores_effects_long)[1] <- "Marker"
              data_fdr <- set.threshold(GWAS.fitted,method="FDR",level=0.05,n.core=cores)
              SigQTL_fdr <- get.QTL(data_fdr)
              colnames(SigQTL_fdr) <- gsub("^scores\\.", "", colnames(SigQTL_fdr))
              if (nrow(SigQTL_fdr) > 0) {
                SigQTL_fdr <- merge(SigQTL_fdr, GWAS_scores_effects_long, by = c("Marker", "Model"))
              }
              if (is.null(SigQTL_fdr) == "TRUE") {print ("file is empty")} else{
                if (is.null(SigQTL_fdr) == "FALSE") {
                  write.table(SigQTL_fdr, paste("./sigFDR/","Significant_effect_",traitname,"_fdr0.05.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
                }}

              data_sugg <- set.threshold(GWAS.fitted,method="FDR",level=0.5,n.core=cores)
              data_sugg <- get.QTL(data_sugg)
              data_sugg <- subset(data_sugg, Score >= threshold_suggestive)
              colnames(data_sugg) <- gsub("^scores\\.", "", colnames(data_sugg))
              if (nrow(data_sugg) > 0) {
                data_sugg <- merge(data_sugg, GWAS_scores_effects_long, by = c("Marker", "Model"))
              }
              write.table(data_sugg, paste("./sigSuggestive/","Significant_effect_",traitname,"_sugg",threshold_suggestive,".txt",sep=""), row.names=F, quote = FALSE, sep = "\t")

              data_Bonferroni <- set.threshold(GWAS.fitted,method="Bonferroni",level=0.05,n.core=cores)
              SigQTL_Bonferroni <- get.QTL(data_Bonferroni)
              colnames(SigQTL_Bonferroni) <- gsub("^scores\\.", "", colnames(SigQTL_Bonferroni))
              if (nrow(SigQTL_Bonferroni) > 0) {
                SigQTL_Bonferroni <- merge(SigQTL_Bonferroni, GWAS_scores_effects_long, by = c("Marker", "Model"))
              }
              if (is.null(SigQTL_Bonferroni) == "TRUE") {print ("file is empty")} else{
                if (is.null(SigQTL_Bonferroni) == "FALSE") {
                  write.table(SigQTL_Bonferroni, paste("./sigBonferroni/","Significant_effect_",traitname,"_Bonferroni0.05.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
                }}

              # data_Meff <- set.threshold(GWAS.fitted,method="M.eff",level=0.05,n.core=cores)
              # SigQTL_Meff <- get.QTL(data_Meff)
              # if (is.null(SigQTL_Meff) == "TRUE") {print ("file is empty")} else{
              #   if (is.null(SigQTL_Meff) == "FALSE") {
              #     write.table(SigQTL_Meff, paste("./sigMeff/","Significant_effect_",traitname,"_Meff.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
              #   }}
              if (is.null(permutations)) { permutations=0}
              if (permutations >= 100) {
                data_permute <- set.threshold(GWAS.fitted,method="permute",n.permute=permutations,level=0.05,n.core=cores)
                SigQTL_permute <- get.QTL(data_permute)
                colnames(SigQTL_permute) <- gsub("^scores\\.", "", colnames(SigQTL_permute))
                if (nrow(SigQTL_permute) > 0) {
                  SigQTL_permute <- merge(SigQTL_permute, GWAS_scores_effects_long, by = c("Marker", "Model"))
                }
                if (is.null(SigQTL_permute) == "TRUE") {print ("file is empty")} else{
                  if (is.null(SigQTL_permute) == "FALSE") {
                    write.table(SigQTL_permute, paste("./sigpermute/","Significant_effect_",traitname,"_permute0.05.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
                  }}
              }

              scores <- GWAS.fitted@scores[[colnames(pheno[2])]]
              colnames(scores) <- gsub("^scores\\.", "", colnames(scores))
              scores <- setDT(scores, keep.rownames = TRUE)
              scores$Chrom <- gsub("_.+$", "", scores$rn); scores$bp <- gsub("^.+_", "", scores$rn); scores$rn <- NULL
              scores <- as.data.frame(reshape2::melt(scores, id=c("Chrom","bp"))); colnames(scores) <- c("Chrom","bp","models","scores")
              scores <-na.omit(scores); scores$scores <- as.numeric(as.character(scores$scores)); scores$bp <- as.numeric(as.character(scores$bp))

              fdr_threshold <- t(data_fdr@threshold)
              fdr_threshold <- setDT(as.data.frame(fdr_threshold), keep.rownames = TRUE)
              colnames(fdr_threshold) <- c("model","threshold")
              hline.fdr <- data.frame(z1 = c(fdr_threshold$threshold), models = c(fdr_threshold$model))

              Bonferroni_threshold <- t(data_Bonferroni@threshold)
              Bonferroni_threshold <- setDT(as.data.frame(Bonferroni_threshold), keep.rownames = TRUE)
              colnames(Bonferroni_threshold) <- c("model","threshold")
              hline.Bonferroni <- data.frame(z2 = c(Bonferroni_threshold$threshold), models = c(Bonferroni_threshold$model))

              suggestive_threshold <- fdr_threshold
              suggestive_threshold$threshold <- rep(threshold_suggestive, nrow(suggestive_threshold))
              hline.suggestive <- data.frame(z1 = c(suggestive_threshold$threshold), models = c(suggestive_threshold$model))


              if (permutations >= 100) {
                permute_threshold <- t(data_permute@threshold)
                permute_threshold <- setDT(as.data.frame(permute_threshold), keep.rownames = TRUE)
                colnames(permute_threshold) <- c("model","threshold")
                hline.permute <- data.frame(z3 = c(permute_threshold$threshold), models = c(permute_threshold$model))
              } else { z3 <- -1 ; hline.permute <- NULL; threshold_permute <- FALSE }
              if (!exists("threshold_suggestive")){z4 <- -1; hline.suggestive_threshold <- NULL }
              if (is.null(threshold_suggestive)){z4 <- -1; hline.suggestive_threshold <- NULL }
              if (exists("threshold_suggestive")){
                hline.suggestive_threshold <- hline.Bonferroni
                hline.suggestive_threshold[,1] <- threshold_suggestive
                colnames(hline.suggestive_threshold)[colnames(hline.suggestive_threshold) == 'z2'] <- 'z4'
              }
              if (threshold_FDR == FALSE) {z1 <- -1 ; hline.fdr <- NULL}
              if (threshold_Bonferroni == FALSE) {z2 <- -1 ; hline.Bonferroni <- NULL}
              scores$Chrom <- gsub("CHR", "", scores$Chrom, ignore.case = TRUE)
              scores$Chrom <- gsub("CHROM", "", scores$Chrom, ignore.case = TRUE)
              scores$Chrom <- as.numeric(scores$Chrom)
              scores <- scores %>%
                arrange(Chrom, bp) %>%
                mutate(Chrom = factor(Chrom, levels = unique(Chrom)),
                       Chrom_bp = interaction(Chrom, bp, sep = "_"))
              scores$Chrom <- factor(scores$Chrom, levels = sort(as.numeric(levels(scores$Chrom))))

              # Create a dummy data frame for legend
              legend_df <- data.frame(x = 1:length(c("FDR", "Bonferroni", "Permutation", "Suggestive")),
                                      y = 1,  # dummy y, won't affect plot
                                      label = c("FDR", "Bonferroni", "Permutation", "Suggestive"),
                                      color = c("grey10", "tomato", "green4", "cornflowerblue"))
              if(is.null(threshold_suggestive)) {threshold_suggested <- FALSE} else {threshold_suggested <- TRUE}
              if(permutations <= 1) {threshold_permuted <- FALSE}
              legend_df$value <- c(threshold_FDR, threshold_Bonferroni, threshold_permuted, threshold_suggested)
              legend_df <- subset(legend_df, value == TRUE)
              legend_df$threshold <- factor(legend_df$label, levels = legend_df$label)



              if (ploidy == 2) {
                if(model == "Add"){
                  dim=10
                  hline.suggestive_threshold <- hline.suggestive_threshold[1,]
                }
                if(model == "Add_Dom"){
                  dim=20
                }
                if (diversity == TRUE) {
                  manplot <- ggplot(scores, aes(x = Chrom_bp, y=scores, group=Chrom)) +
                    geom_line(aes(group=1), lwd=1.0,color="grey20") + geom_point(size=10, pch=20, aes(colour=factor(Chrom)), show.legend = FALSE) + theme(text = element_text(size=40)) +
                    scale_color_manual(values=c('cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange')) +
                    geom_hline(aes(yintercept = z1), hline.fdr, color="grey10", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z2), hline.Bonferroni, color="tomato", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z3), hline.permute, color="green4", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z4), hline.suggestive_threshold, color="cornflowerblue", size=1.0, linetype="dashed")+
                    facet_grid(models ~ Chrom, scales="free_x", space="free_x") + xlab("\nChromosomes and Physical Map Position (bp)")+ ylab("-log10(P)\n") +
                    theme(panel.spacing = unit(0.5, "lines"), panel.grid.major = element_line(colour ="grey95"), panel.background = element_blank(),
                          panel.border = element_blank(), strip.background = element_blank(), axis.text.x=element_blank(),
                          strip.text.x = element_text(size=30,color="black"), strip.text.y = element_text(size=40,color="black"),
                          axis.line=element_line(colour="white")) +
                    theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0,NA)) +
                    # Invisible points used only to construct the legend
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape = 22, alpha = 0, show.legend = TRUE) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(colour="none", fill = guide_legend(title = "Thresholds: ", nrow = 1, keywidth = 3, keyheight = 3, override.aes = list(alpha = 1, shape = 22, size = 10))) +
                    theme(legend.position = "top", legend.justification = "right",  legend.box.just = "right", legend.background = element_rect(fill = "white", color = "white")) +
                    coord_cartesian(clip = "off") +
                    labs(title= paste(traitname,"\n",sep="")) + theme(plot.title = element_text(hjust = 0, vjust=-10))
                  ggsave(file=paste("./manplots/","manplot_",traitname,".tiff",sep=""), plot=manplot, width=30, height=dim, units=("in"), dpi=600, compression = "lzw")
                }
                if (biparental == TRUE) {
                  manplot <- ggplot(scores, aes(x = Chrom_bp, y=scores, group=Chrom)) +
                    geom_line(aes(colour=factor(Chrom), y=rollmean(scores, 2, na.pad=TRUE)), size=3, show.legend = FALSE) + theme(text = element_text(size=40)) +
                    geom_point(size=1, pch=".", aes(colour=factor(Chrom)), show.legend = FALSE)+
                    scale_color_manual(values=c('cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange')) +
                    geom_hline(aes(yintercept = z1), hline.fdr, color="grey50", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z2), hline.Bonferroni, color="tomato", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z3), hline.permute, color="green4", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z4), hline.suggestive_threshold, color="cornflowerblue", size=1.0, linetype="dashed")+
                    facet_grid(models ~ Chrom, scales="free_x", space="free_x") + xlab("\nChromosomes and Physical Map Position (bp)")+ ylab("-log10(P)\n") +
                    theme(panel.spacing = unit(0.5, "lines"), panel.grid.major = element_line(colour ="grey95"), panel.background = element_blank(),
                          panel.border = element_blank(), strip.background = element_blank(), axis.text.x=element_blank(),
                          strip.text.x = element_text(size=30,color="black"), strip.text.y = element_text(size=40,color="black"),
                          axis.line=element_line(colour="white")) +
                    theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0,NA)) +
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape = 22, alpha = 0, show.legend = TRUE) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(colour="none", fill = guide_legend(title = "Thresholds: ", nrow = 1, keywidth = 3, keyheight = 3, override.aes = list(alpha = 1, shape = 22, size = 10))) +
                    theme(legend.position = "top", legend.justification = "right",  legend.box.just = "right", legend.background = element_rect(fill = "white", color = "white")) +
                    coord_cartesian(clip = "off") +
                    labs(title= paste(traitname,"\n",sep="")) + theme(plot.title = element_text(hjust = 0, vjust=-10))
                  ggsave(file=paste("./manplots/","QTLprofile_",traitname,".tiff",sep=""), plot=manplot, width=30, height=dim, units=("in"), dpi=600, compression = "lzw")
                }
              }
              if (ploidy == 4) {
                if (diversity == TRUE) {
                  manplot <- ggplot(scores, aes(x = Chrom_bp, y=scores, group=Chrom)) +
                    geom_line(lwd=1.0,color="grey20") + geom_point(size=10, pch=20, aes(colour=factor(Chrom)), show.legend = FALSE)+ theme(text = element_text(size=40)) +
                    scale_color_manual(values=c('cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange')) +
                    geom_hline(aes(yintercept = z1), hline.fdr, color="grey10", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z2), hline.Bonferroni, color="tomato", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z3), hline.permute, color="green4", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z4), hline.suggestive_threshold, color="cornflowerblue", size=1.0, linetype="dashed")+
                    facet_grid(models ~ Chrom, scales="free_x", space="free_x") + xlab("\nChromosomes and Physical Map Position (bp)")+ ylab("-log10(P)\n") +
                    theme(panel.spacing = unit(0.5, "lines"), panel.grid.major = element_line(colour ="grey95"), panel.background = element_blank(),
                          panel.border = element_blank(), strip.background = element_blank(), axis.text.x=element_blank(),
                          strip.text.x = element_text(size=30,color="black"), strip.text.y = element_text(size=40,color="black"),
                          axis.line=element_line(colour="white")) +
                    theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0,NA)) +
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape = 22, alpha = 0, show.legend = TRUE) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(colour="none", fill = guide_legend(title = "Thresholds: ", nrow = 1, keywidth = 3, keyheight = 3, override.aes = list(alpha = 1, shape = 22, size = 10))) +
                    theme(legend.position = "top", legend.justification = "right",  legend.box.just = "right", legend.background = element_rect(fill = "white", color = "white")) +
                    coord_cartesian(clip = "off") +
                    labs(title= paste(traitname,"\n",sep="")) + theme(plot.title = element_text(hjust = 0, vjust=-10))
                  ggsave(file=paste("./manplots/","manplot_",traitname,".tiff",sep=""), plot=manplot, width=30, height=30, units=("in"), dpi=300, compression = "lzw")
                }
                if (biparental == TRUE) {
                  manplot <- ggplot(scores, aes(x = Chrom_bp, y=scores, group=Chrom)) +
                    geom_line(aes(colour=factor(Chrom), y=rollmean(scores, 2, na.pad=TRUE)), size=3, show.legend = FALSE) + theme(text = element_text(size=40)) +
                    geom_point(size=1, pch=".", aes(colour=factor(Chrom)), show.legend = FALSE)+
                    scale_color_manual(values=c('cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange')) +
                    geom_hline(aes(yintercept = z1), hline.fdr, color="grey50", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z2), hline.Bonferroni, color="tomato", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z3), hline.permute, color="green4", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z4), hline.suggestive_threshold, color="cornflowerblue", size=1.0, linetype="dashed")+
                    facet_grid(models ~ Chrom, scales="free_x", space="free_x") + xlab("\nChromosomes and Physical Map Position (bp)")+ ylab("-log10(P)\n") +
                    theme(panel.spacing = unit(0.5, "lines"), panel.grid.major = element_line(colour ="grey95"), panel.background = element_blank(),
                          panel.border = element_blank(), strip.background = element_blank(), axis.text.x=element_blank(),
                          strip.text.x = element_text(size=30,color="black"), strip.text.y = element_text(size=40,color="black"),
                          axis.line=element_line(colour="white")) +
                    theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0,NA)) +
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape = 22, alpha = 0, show.legend = TRUE) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(colour="none", fill = guide_legend(title = "Thresholds: ", nrow = 1, keywidth = 3, keyheight = 3, override.aes = list(alpha = 1, shape = 22, size = 10))) +
                    theme(legend.position = "top", legend.justification = "right",  legend.box.just = "right", legend.background = element_rect(fill = "white", color = "white")) +
                    coord_cartesian(clip = "off") +
                    labs(title= paste(traitname,"\n",sep="")) + theme(plot.title = element_text(hjust = 0, vjust=-10))
                  ggsave(file=paste("./manplots/","QTLprofile_",traitname,".tiff",sep=""), plot=manplot, width=30, height=30, units=("in"), dpi=300, compression = "lzw")
                }
              }
              if (ploidy == 6) {
                if (diversity == TRUE) {
                  manplot <- ggplot(scores, aes(x = Chrom_bp, y=scores, group=Chrom)) +
                    geom_line(lwd=1.0,color="grey20") + geom_point(size=10, pch=20, aes(colour=factor(Chrom)), show.legend = FALSE)+ theme(text = element_text(size=40)) +
                    scale_color_manual(values=c('cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange')) +
                    geom_hline(aes(yintercept = z1), hline.fdr, color="grey10", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z2), hline.Bonferroni, color="tomato", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z3), hline.permute, color="green4", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z4), hline.suggestive_threshold, color="cornflowerblue", size=1.0, linetype="dashed")+
                    facet_grid(models ~ Chrom, scales="free_x", space="free_x") + xlab("\nChromosomes and Physical Map Position (bp)")+ ylab("-log10(P)\n") +
                    theme(panel.spacing = unit(0.5, "lines"), panel.grid.major = element_line(colour ="grey95"), panel.background = element_blank(),
                          panel.border = element_blank(), strip.background = element_blank(), axis.text.x=element_blank(),
                          strip.text.x = element_text(size=30,color="black"), strip.text.y = element_text(size=40,color="black"),
                          axis.line=element_line(colour="white")) +
                    theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0,NA)) +
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape = 22, alpha = 0, show.legend = TRUE) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(colour="none", fill = guide_legend(title = "Thresholds: ", nrow = 1, keywidth = 3, keyheight = 3, override.aes = list(alpha = 1, shape = 22, size = 10))) +
                    theme(legend.position = "top", legend.justification = "right",  legend.box.just = "right", legend.background = element_rect(fill = "white", color = "white")) +
                    coord_cartesian(clip = "off") +
                    labs(title= paste(traitname,"\n",sep="")) + theme(plot.title = element_text(hjust = 0, vjust=-10))
                  ggsave(file=paste("./manplots/","manplot_",traitname,".tiff",sep=""), plot=manplot, width=30, height=40, units=("in"), dpi=300, compression = "lzw")
                }
                if (biparental == TRUE) {
                  manplot <- ggplot(scores, aes(x = Chrom_bp, y=scores, group=Chrom)) +
                    geom_line(aes(colour=factor(Chrom), y=rollmean(scores, 2, na.pad=TRUE)), size=3, show.legend = FALSE) + theme(text = element_text(size=40)) +
                    geom_point(size=1, pch=".", aes(colour=factor(Chrom)), show.legend = FALSE)+
                    scale_color_manual(values=c('cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange')) +
                    geom_hline(aes(yintercept = z1), hline.fdr, color="grey50", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z2), hline.Bonferroni, color="tomato", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z3), hline.permute, color="green4", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z4), hline.suggestive_threshold, color="cornflowerblue", size=1.0, linetype="dashed")+
                    facet_grid(models ~ Chrom, scales="free_x", space="free_x") + xlab("\nChromosomes and Physical Map Position (bp)")+ ylab("-log10(P)\n") +
                    theme(panel.spacing = unit(0.5, "lines"), panel.grid.major = element_line(colour ="grey95"), panel.background = element_blank(),
                          panel.border = element_blank(), strip.background = element_blank(), axis.text.x=element_blank(),
                          strip.text.x = element_text(size=30,color="black"), strip.text.y = element_text(size=40,color="black"),
                          axis.line=element_line(colour="white")) +
                    theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0,NA)) +
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape = 22, alpha = 0, show.legend = TRUE) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(colour="none", fill = guide_legend(title = "Thresholds: ", nrow = 1, keywidth = 3, keyheight = 3, override.aes = list(alpha = 1, shape = 22, size = 10))) +
                    theme(legend.position = "top", legend.justification = "right",  legend.box.just = "right", legend.background = element_rect(fill = "white", color = "white")) +
                    coord_cartesian(clip = "off") +
                    labs(title= paste(traitname,"\n",sep="")) + theme(plot.title = element_text(hjust = 0, vjust=-10))
                  ggsave(file=paste("./manplots/","QTLprofile_",traitname,".tiff",sep=""), plot=manplot, width=30, height=40, units=("in"), dpi=300, compression = "lzw")
                }
              }
              if (ploidy == 8) {
                if (diversity == TRUE) {
                  manplot <- ggplot(scores, aes(x = Chrom_bp, y=scores, group=Chrom)) +
                    geom_line(lwd=1.0,color="grey20") + geom_point(size=10, pch=20, aes(colour=factor(Chrom)), show.legend = FALSE)+ theme(text = element_text(size=40)) +
                    scale_color_manual(values=c('cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange')) +
                    geom_hline(aes(yintercept = z1), hline.fdr, color="grey10", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z2), hline.Bonferroni, color="tomato", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z3), hline.permute, color="green4", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z4), hline.suggestive_threshold, color="cornflowerblue", size=1.0, linetype="dashed")+
                    facet_grid(models ~ Chrom, scales="free_x", space="free_x") + xlab("\nChromosomes and Physical Map Position (bp)")+ ylab("-log10(P)\n") +
                    theme(panel.spacing = unit(0.5, "lines"), panel.grid.major = element_line(colour ="grey95"), panel.background = element_blank(),
                          panel.border = element_blank(), strip.background = element_blank(), axis.text.x=element_blank(),
                          strip.text.x = element_text(size=30,color="black"), strip.text.y = element_text(size=40,color="black"),
                          axis.line=element_line(colour="white")) +
                    theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0,NA)) +
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape = 22, alpha = 0, show.legend = TRUE) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(colour="none", fill = guide_legend(title = "Thresholds: ", nrow = 1, keywidth = 3, keyheight = 3, override.aes = list(alpha = 1, shape = 22, size = 10))) +
                    theme(legend.position = "top", legend.justification = "right",  legend.box.just = "right", legend.background = element_rect(fill = "white", color = "white")) +
                    coord_cartesian(clip = "off") +
                    labs(title= paste(traitname,"\n",sep="")) + theme(plot.title = element_text(hjust = 0, vjust=-10))
                  ggsave(file=paste("./manplots/","manplot_",traitname,".tiff",sep=""), plot=manplot, width=30, height=50, units=("in"), dpi=300, compression = "lzw")
                }
                if (biparental == TRUE) {
                  manplot <- ggplot(scores, aes(x = Chrom_bp, y=scores, group=Chrom)) +
                    geom_line(aes(colour=factor(Chrom), y=rollmean(scores, 2, na.pad=TRUE)), size=3, show.legend = FALSE) + theme(text = element_text(size=40)) +
                    geom_point(size=1, pch=".", aes(colour=factor(Chrom)), show.legend = FALSE)+
                    scale_color_manual(values=c('cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                                'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange')) +
                    geom_hline(aes(yintercept = z1), hline.fdr, color="grey50", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z2), hline.Bonferroni, color="tomato", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z3), hline.permute, color="green4", size=1.0, linetype="dashed")+
                    geom_hline(aes(yintercept = z4), hline.suggestive_threshold, color="cornflowerblue", size=1.0, linetype="dashed")+
                    facet_grid(models ~ Chrom, scales="free_x", space="free_x") + xlab("\nChromosomes and Physical Map Position (bp)")+ ylab("-log10(P)\n") +
                    theme(panel.spacing = unit(0.5, "lines"), panel.grid.major = element_line(colour ="grey95"), panel.background = element_blank(),
                          panel.border = element_blank(), strip.background = element_blank(), axis.text.x=element_blank(),
                          strip.text.x = element_text(size=30,color="black"), strip.text.y = element_text(size=40,color="black"),
                          axis.line=element_line(colour="white")) +
                    theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0,NA)) +
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape = 22, alpha = 0, show.legend = TRUE) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(colour="none", fill = guide_legend(title = "Thresholds: ", nrow = 1, keywidth = 3, keyheight = 3, override.aes = list(alpha = 1, shape = 22, size = 10))) +
                    theme(legend.position = "top", legend.justification = "right",  legend.box.just = "right", legend.background = element_rect(fill = "white", color = "white")) +
                    coord_cartesian(clip = "off") +
                    labs(title= paste(traitname,"\n",sep="")) + theme(plot.title = element_text(hjust = 0, vjust=-10))
                  ggsave(file=paste("./manplots/","QTLprofile_",traitname,".tiff",sep=""), plot=manplot, width=30, height=50, units=("in"), dpi=300, compression = "lzw")
                }
              }

              hline.fdr <- NULL
              hline.Bonferroni <- NULL
              hline.permute <- NULL
              z1 <- NULL; z2 <- NULL; z3 <- NULL
            }
            tcorr <- NULL
          }
          delete_empty_dirs <- function(path = "./") {
            dirs <- list.dirs(path, recursive = TRUE, full.names = TRUE)
            empty_dirs <- dirs[sapply(dirs, function(d) length(list.files(d, all.files = TRUE, no.. = TRUE)) == 0)]
            invisible(lapply(empty_dirs, unlink, recursive = TRUE, force = TRUE))
            return(empty_dirs)
          }
          deleted <- delete_empty_dirs("./")
          print(deleted)
          setwd("../")
        }
      }
    }
  }


  # p <- read.table("./FER_GWAS_pheno_nocov_2x/pvalues/logPNI_2021.txt", sep="\t", head = TRUE)
  # for (i in c(1:ncol(p))){
  #   q <- qvalue(10^-(p[,i]))
  #   q_value <- as.data.frame(q$qvalues)
  #   p[,paste(i)] <- q_value
  # }
}


# Resolve input paths once, before the GWAS engine changes directories.
.sgwas_path <- function(path, base) {
  if(is.null(path) || identical(path,"NULL")) return(NULL)
  if(length(path)!=1L || !is.character(path) || is.na(path) || !nzchar(path))
    stop("Expected one nonempty input file path")
  if(!grepl("^(/|[A-Za-z]:[/\\\\]|~)",path)) path <- file.path(base,path)
  normalizePath(path.expand(path),mustWork=TRUE)
}

# Integrated association models: additive dosage, complete-case cohorts,
# REML null covariance estimates and P3D/GLS marker tests. Base R only.
.sgwas_write_tsv <- function(x, file) {
  write.table(x, file, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
}

.sgwas_integer <- function(x, name, minimum=1L) {
  if(length(x)!=1L || !is.numeric(x) || !is.finite(x) || x<minimum || x!=floor(x))
    stop(name," must be an integer >= ",minimum)
  as.integer(x)
}

.sgwas_read_dosage <- function(file, ploidy) {
  d <- read.delim(file, check.names=FALSE, stringsAsFactors=FALSE, colClasses="character")
  if(ncol(d)<7L || !identical(names(d)[1:5], c("SNP","CHROM","POS","REF","ALT")))
    stop("Integrated GWAS needs SNP, CHROM, POS, REF, ALT, then at least two sample columns")
  if(anyDuplicated(names(d)) || any(!nzchar(names(d)))) stop("Genotype columns must be unique and nonempty")
  if(!nrow(d) || anyNA(d$SNP) || any(!nzchar(d$SNP)) || anyDuplicated(d$SNP))
    stop("SNP IDs must be nonempty and unique")
  if(anyNA(d$CHROM) || any(!nzchar(d$CHROM))) stop("Every marker needs a chromosome")
  d$POS <- .sgwas_numeric(d$POS,"POS")
  if(any(!is.finite(d$POS)) || any(d$POS<0)) stop("Marker positions must be finite and nonnegative")
  G <- vapply(d[-(1:5)], .sgwas_numeric, numeric(nrow(d)), label="marker dosage")
  G <- t(matrix(G,nrow=nrow(d),dimnames=list(d$SNP,names(d)[-(1:5)])))
  if(any(!is.na(G) & (!is.finite(G) | G<0 | G>ploidy)))
    stop("Observed dosages must be finite and between 0 and ploidy")
  list(G=G, map=d[1:5])
}

.sgwas_dosage_qc <- function(raw, ids, ploidy, maf, max_missing) {
  G <- raw$G[ids,,drop=FALSE]
  missing <- colMeans(is.na(G)); af <- colMeans(G,na.rm=TRUE)/ploidy
  marker_maf <- pmin(af,1-af)
  variable <- vapply(seq_len(ncol(G)),function(j) {
    v <- G[,j]; v <- v[is.finite(v)]; length(v)>1L && diff(range(v))>1e-10
  },logical(1))
  keep <- is.finite(marker_maf) & marker_maf>=maf & missing<=max_missing & variable
  ledger <- cbind(raw$map, MAF=marker_maf, Missing_Rate=missing, Variable=variable, Retained=keep)
  if(!any(keep)) stop("No variable markers pass MAF and missingness filters on the analysed cohort")
  G <- G[,keep,drop=FALSE]; mu <- colMeans(G,na.rm=TRUE)
  for(j in which(colSums(is.na(G))>0)) G[is.na(G[,j]),j] <- mu[j]
  list(G=G, centered=sweep(G,2,mu,"-"), map=ledger[keep,,drop=FALSE], ledger=ledger)
}

.sgwas_kinship_eigen <- function(crossproduct) {
  K <- (crossproduct+t(crossproduct))/2
  scale <- mean(diag(K))
  if(!is.finite(scale) || scale<=1e-12) stop("No variable background markers are available for kinship")
  K <- K/scale
  eg <- eigen(K,symmetric=TRUE)
  if(min(eg$values)< -1e-7*max(eg$values)) stop("Kinship is not positive semidefinite")
  eg$values <- pmax(eg$values,0)
  eg
}

.sgwas_inverse <- function(A) {
  ch <- chol((A+t(A))/2)
  list(inverse=chol2inv(ch), logdet=2*sum(log(diag(ch))))
}

# Common normal-equation engine; observations are rotated by eigenvectors of K.
# For two traits, the observation order is all trait 1, then all trait 2.
.sgwas_gls_null <- function(Y, X, weights) {
  Y <- as.matrix(Y); q <- ncol(X); d <- ncol(Y)
  if(d==1L) {
    A <- crossprod(X,X*weights[[1]])
    r <- crossprod(X,Y[,1]*weights[[1]])
    yy <- sum(Y[,1]^2*weights[[1]])
  } else {
    a <- weights[[1]]; b <- weights[[2]]; c <- weights[[3]]
    A12 <- crossprod(X,X*b)
    A <- rbind(cbind(crossprod(X,X*a),A12),cbind(A12,crossprod(X,X*c)))
    r <- c(crossprod(X,a*Y[,1]+b*Y[,2]),crossprod(X,b*Y[,1]+c*Y[,2]))
    yy <- sum(a*Y[,1]^2+2*b*Y[,1]*Y[,2]+c*Y[,2]^2)
  }
  inv <- .sgwas_inverse(A); beta <- inv$inverse%*%r
  residual <- Y-X%*%matrix(beta,nrow=q,ncol=d)
  Py <- if(d==1L) matrix(weights[[1]]*residual[,1],ncol=1) else
    cbind(weights[[1]]*residual[,1]+weights[[2]]*residual[,2],
          weights[[2]]*residual[,1]+weights[[3]]*residual[,2])
  rss <- sum(residual*Py)
  if(!is.finite(rss) || rss<=0) stop("Nonpositive residual sum of squares")
  list(X=X,Y=Y,weights=weights,Ainv=inv$inverse,beta=beta,Py=Py,rss=rss,
       logdetA=inv$logdet,nobs=length(Y),nfixed=q*d)
}

.sgwas_fit_univariate <- function(y, X, eig) {
  y <- as.numeric(y); U <- eig$vectors; lambda <- eig$values
  yr <- crossprod(U,y); xr <- crossprod(U,X); df <- length(y)-ncol(X)
  profile <- function(h, details=FALSE) {
    v <- (1-h)+h*lambda
    fit <- .sgwas_gls_null(yr,xr,list(1/v))
    s <- fit$rss/df
    value <- .5*(df*(log(2*pi)+1+log(s))+sum(log(v))+fit$logdetA)
    if(details) return(list(fit=fit,variance=s,h=h,objective=value))
    value
  }
  upper <- 1-1e-6
  opt <- optimize(profile,c(0,upper),tol=1e-8)
  candidates <- c(0,opt$minimum,upper)
  values <- vapply(candidates,profile,numeric(1))
  best <- profile(candidates[which.min(values)],TRUE)
  best$G <- matrix(best$variance*best$h,1,1)
  best$R <- matrix(best$variance*(1-best$h),1,1)
  best$boundary <- best$h<1e-5 || best$h>upper-1e-5
  best$converged <- TRUE; best$optimizer <- "profile_REML"
  best
}

.sgwas_bivariate_weights <- function(lambda,G,R) {
  a <- lambda*G[1,1]+R[1,1]; b <- lambda*G[1,2]+R[1,2]
  c <- lambda*G[2,2]+R[2,2]; det <- a*c-b*b
  if(any(!is.finite(det)) || any(det<=0)) stop("Nonpositive bivariate covariance determinant")
  list(weights=list(c/det,-b/det,a/det),logdet=sum(log(det)))
}

.sgwas_cov2 <- function(par) {
  v <- exp(par[1:2]); cv <- tanh(par[3])*sqrt(prod(v))
  matrix(c(v[1],cv,cv,v[2]),2,2)
}

.sgwas_fit_bivariate <- function(Y,X,eig,maxit=300L) {
  centers <- colMeans(Y); scales <- apply(Y,2,sd)
  Ys <- sweep(sweep(Y,2,centers,"-"),2,scales,"/")
  yr <- crossprod(eig$vectors,Ys); xr <- crossprod(eig$vectors,X)
  objective <- function(par,details=FALSE) {
    G <- .sgwas_cov2(par[1:3]); R <- .sgwas_cov2(par[4:6])
    ans <- tryCatch({
      w <- .sgwas_bivariate_weights(eig$values,G,R)
      fit <- .sgwas_gls_null(yr,xr,w$weights)
      value <- .5*((fit$nobs-fit$nfixed)*log(2*pi)+w$logdet+fit$logdetA+fit$rss)
      if(details) list(fit=fit,G=G,R=R,objective=value) else value
    },error=function(e) if(details) stop(e) else 1e100)
    ans
  }
  # K proportional to I cannot distinguish genetic from residual covariance.
  # Fail explicitly rather than report arbitrary genetic correlations.
  if(diff(range(eig$values))<1e-8)
    stop("Bivariate genetic/residual covariance is not identifiable: kinship is proportional to identity")
  residual <- qr.resid(qr(X),Ys)
  vr <- crossprod(residual)/(nrow(Y)-ncol(X))
  if(min(eigen(vr,symmetric=TRUE,only.values=TRUE)$values)<1e-8)
    stop("The two outcomes are linearly dependent after covariate adjustment")
  corr <- max(-.95,min(.95,cov2cor(vr)[1,2]))
  fractions <- c(.2,.5,.8)
  starts <- lapply(fractions,function(f) c(log(diag(vr)*f),atanh(corr),
                                          log(diag(vr)*(1-f)),atanh(corr)))
  fits <- lapply(starts,function(s) tryCatch(nlminb(s,objective,
    lower=c(-18,-18,-6,-18,-18,-6),upper=c(8,8,6,8,8,6),
    control=list(iter.max=maxit,eval.max=maxit*3L,rel.tol=1e-9)),error=identity))
  good <- vapply(fits,function(f) !inherits(f,"error") && f$convergence==0L &&
                  is.finite(f$objective) && f$objective<1e90,logical(1))
  attempts <- do.call(rbind,lapply(seq_along(fits),function(i) {
    f <- fits[[i]]
    if(inherits(f,"error")) data.frame(Start=i,Convergence=NA,Objective=NA,Message=conditionMessage(f)) else
      data.frame(Start=i,Convergence=f$convergence,Objective=f$objective,Message=f$message)
  }))
  if(!any(good)) stop("Bivariate REML did not converge; increase association_maxit. ",
                       paste(attempts$Message,collapse="; "))
  best <- fits[[which(good)[which.min(vapply(fits[good],`[[`,numeric(1),"objective"))]]]
  ans <- objective(best$par,TRUE)
  ans$centers <- centers; ans$scales <- scales; ans$converged <- TRUE
  ans$optimizer <- "nlminb_REML_3_starts"; ans$attempts <- attempts
  ans$boundary <- any(best$par[c(1,2,4,5)]< -15) || any(abs(best$par[c(3,6)])>5.5)
  ans$G_original <- ans$G*outer(scales,scales); ans$R_original <- ans$R*outer(scales,scales)
  ans
}

# Test one or two marker columns after projecting out the null fixed effects.
# The covariance shape is held fixed; the scalar residual variance is estimated
# in the full marker model. F tests are approximate because shape was estimated.
.sgwas_marker_test <- function(fit,g,h=NULL,bivariate=FALSE) {
  g <- as.matrix(g); X <- fit$X; w <- fit$weights; m <- ncol(g)
  two <- bivariate || !is.null(h)
  if(bivariate) {
    B1 <- rbind(crossprod(X,g*w[[1]]),crossprod(X,g*w[[2]]))
    B2 <- rbind(crossprod(X,g*w[[2]]),crossprod(X,g*w[[3]]))
    S11 <- colSums(g*g*w[[1]]); S12 <- colSums(g*g*w[[2]]); S22 <- colSums(g*g*w[[3]])
    t1 <- as.numeric(crossprod(g,fit$Py[,1])); t2 <- as.numeric(crossprod(g,fit$Py[,2]))
  } else {
    B1 <- crossprod(X,g*w[[1]]); S11 <- colSums(g*g*w[[1]])
    t1 <- as.numeric(crossprod(g,fit$Py[,1]))
    if(two) {
      h <- as.matrix(h); B2 <- crossprod(X,h*w[[1]])
      S12 <- colSums(g*h*w[[1]]); S22 <- colSums(h*h*w[[1]])
      t2 <- as.numeric(crossprod(h,fit$Py[,1]))
    }
  }
  S11 <- S11-colSums(B1*(fit$Ainv%*%B1))
  df <- fit$nobs-fit$nfixed-if(two) 2L else 1L
  if(df<=0L) stop("Insufficient residual degrees of freedom")
  if(two) {
    S12 <- S12-colSums(B1*(fit$Ainv%*%B2))
    S22 <- S22-colSums(B2*(fit$Ainv%*%B2))
    det <- S11*S22-S12*S12
    valid <- is.finite(det) & S11>1e-10 & S22>1e-10 & det>1e-10*S11*S22
    det[!valid] <- NA_real_
    b1 <- (S22*t1-S12*t2)/det; b2 <- (S11*t2-S12*t1)/det
    explained <- pmax(0,b1*t1+b2*t2)
    residual <- pmax(0,fit$rss-explained)
    variance <- residual/df
    se1 <- sqrt(variance*S22/det); se2 <- sqrt(variance*S11/det)
    pjoint <- pf(explained/2/variance,2,df,lower.tail=FALSE)
    p1 <- pf((b1/se1)^2,1,df,lower.tail=FALSE); p2 <- pf((b2/se2)^2,1,df,lower.tail=FALSE)
  } else {
    valid <- is.finite(S11) & S11>1e-10; S11[!valid] <- NA_real_
    b1 <- t1/S11; explained <- pmax(0,b1*t1); variance <- pmax(0,fit$rss-explained)/df
    se1 <- sqrt(variance/S11); p1 <- pf((b1/se1)^2,1,df,lower.tail=FALSE)
    b2 <- se2 <- p2 <- pjoint <- rep(NA_real_,m)
  }
  valid <- valid & is.finite(variance) & variance>.Machine$double.eps*fit$rss/df
  out <- data.frame(Beta1=b1,SE1=se1,P1=p1,Beta2=b2,SE2=se2,P2=p2,PJoint=pjoint,
                    DF_Residual=df,Status=ifelse(valid,"ok","rank_or_residual_failure"))
  out[!valid,1:7] <- NA_real_
  out
}

.sgwas_association_pair <- function(Y,X,G,map,eig,mode,chunk_size,maxit) {
  fits <- list(Trait=.sgwas_fit_univariate(Y[,1],X,eig),
               Microbial=.sgwas_fit_univariate(Y[,2],X,eig))
  if(mode %in% c("joint","both")) fits$Bivariate <- .sgwas_fit_bivariate(Y,X,eig,maxit)
  M <- as.numeric(scale(Y[,2])); XX <- cbind(X,Microbial_score=M)
  if(mode %in% c("interaction","both")) {
    if(qr(XX)$rank<ncol(XX)) stop("Microbial score is collinear with a phenotype covariate or genotype PC")
    fits$Interaction <- .sgwas_fit_univariate(Y[,1],XX,eig)
  }
  res <- vector("list",ceiling(ncol(G)/chunk_size))
  for(k in seq_along(res)) {
    js <- seq.int((k-1L)*chunk_size+1L,min(k*chunk_size,ncol(G)))
    g <- G[,js,drop=FALSE]; gr <- crossprod(eig$vectors,g)
    a <- .sgwas_marker_test(fits$Trait$fit,gr); b <- .sgwas_marker_test(fits$Microbial$fit,gr)
    r <- cbind(map[js,,drop=FALSE],
      Beta_Trait=a$Beta1,SE_Trait=a$SE1,P_Trait=a$P1,Status_FER=a$Status,
      Beta_Microbial=b$Beta1,SE_Microbial=b$SE1,P_Microbial=b$P1,Status_Microbial=b$Status)
    if(!is.null(fits$Bivariate)) {
      z <- .sgwas_marker_test(fits$Bivariate$fit,gr,bivariate=TRUE); s <- fits$Bivariate$scales
      r <- cbind(r,Beta_Trait_joint=z$Beta1*s[1],SE_Trait_joint=z$SE1*s[1],P_Trait_joint=z$P1,
                 Beta_Microbial_joint=z$Beta2*s[2],SE_Microbial_joint=z$SE2*s[2],P_Microbial_joint=z$P2,
                 P_Bivariate=z$PJoint,Status_Bivariate=z$Status)
    }
    if(!is.null(fits$Interaction)) {
      z <- .sgwas_marker_test(fits$Interaction$fit,gr,crossprod(eig$vectors,g*M))
      r <- cbind(r,Beta_SNP_at_mean=z$Beta1,SE_SNP_at_mean=z$SE1,P_SNP_at_mean=z$P1,
                 Beta_SNPxM=z$Beta2,SE_SNPxM=z$SE2,P_SNPxM=z$P2,
                 P_SNP_joint=z$PJoint,Status_Interaction=z$Status)
    }
    res[[k]] <- r
  }
  list(results=do.call(rbind,res),fits=fits,microbial_center=mean(Y[,2]),microbial_scale=sd(Y[,2]))
}

.sgwas_association_plots <- function(results,file,fdr,bonferroni,suggestive,total_tests) {
  tests <- names(results)[grepl("^P_",names(results)) & !grepl("_(BH|Bonferroni)$",names(results))]
  pdf(file,width=11,height=5); on.exit(dev.off(),add=TRUE)
  for(test in tests) {
    p <- results[[test]]; use <- is.finite(p) & p>=0 & p<=1
    if(!any(use)) next
    d <- results[use,,drop=FALSE]; p <- p[use]
    chr <- unique(d$CHROM); ix <- order(match(d$CHROM,chr),d$POS)
    ypos <- -log10(pmax(p[ix],.Machine$double.xmin))
    colors <- c("#276D83","#84A98C")[1+(match(d$CHROM[ix],chr)%%2)]
    par(mfrow=c(1,2),mar=c(5,4,3,1))
    plot(seq_along(ix),ypos,pch=20,cex=.55,col=colors,xaxt="n",xlab="Chromosome",ylab="-log10(p)",main=test)
    ticks <- vapply(chr,function(ch) mean(which(d$CHROM[ix]==ch)),numeric(1))
    axis(1,at=ticks,labels=chr,las=2,cex.axis=.7)
    if(bonferroni && total_tests[[test]]>0) abline(h=-log10(.05/total_tests[[test]]),col="#C14953",lty=2)
    sig <- results[[paste0(test,"_BH")]][use]<=.05
    if(fdr && any(sig,na.rm=TRUE)) abline(h=-log10(max(p[sig],na.rm=TRUE)),col="#6E44AA",lty=3)
    if(is.finite(suggestive)) abline(h=suggestive,col="grey60",lty=3)
    observed <- sort(-log10(pmax(p,.Machine$double.xmin)))
    expected <- sort(-log10(ppoints(length(p))))
    plot(expected,observed,pch=20,cex=.6,xlab="Expected -log10(p)",ylab="Observed -log10(p)",main="QQ plot")
    abline(0,1,col="grey50")
  }
}

.sgwas_integrated_gwas <- function(args,prepared,mode,n_pcs,max_missing,chunk_size,maxit) {
  if(!identical(.sgwas_names(args$method),"MLM") || !identical(.sgwas_names(args$model_effect),"Add"))
    stop("joint/interaction/both currently require method = MLM and model_effect = Add")
  if(length(args$perm)!=1L || !is.numeric(args$perm) || !is.finite(args$perm) || args$perm<0 || args$perm!=floor(args$perm))
    stop("perm must be a nonnegative integer")
  if(args$perm>=100) stop("Permutation thresholds are not implemented for integrated GWAS; use perm = 1")
  for(nm in c("LOCO","fdr","bonferroni"))
    if(length(args[[nm]])!=1L || !is.logical(args[[nm]]) || is.na(args[[nm]])) stop(nm," must be TRUE or FALSE")
  if(length(args$suggestive)!=1L || !is.numeric(args$suggestive) || !is.finite(args$suggestive) || args$suggestive<0)
    stop("suggestive must be a finite nonnegative number")
  if(length(args$maf)!=1L || !is.numeric(args$maf) || !is.finite(args$maf) || args$maf<0 || args$maf>.5)
    stop("maf must be between 0 and 0.5")
  traits <- .sgwas_names(args$trait_names); covs <- .sgwas_names(args$covariate_pheno)
  if(!length(traits)) stop("Supply trait_names, e.g. FER, for integrated GWAS")
  if(any(traits %in% covs)) stop("An outcome cannot be its own covariate")
  ph <- .sgwas_phenotypes(args$phenofile)
  if(is.null(ph)) stop("Integrated GWAS requires phenofile with the measured outcome")
  if(any(!c(traits,covs) %in% names(ph))) stop("Some requested traits/covariates are absent from phenofile")
  if(prepared$summary$Trait %in% c(traits,covs)) stop("Use a microbial_trait_name distinct from outcomes and covariates")
  if(args$cores>1) message("[systemsGWAS] Integrated scans use serial vectorised chunks; cores controls the legacy GWASpoly engine only.")
  provenance <- tolower(prepared$summary$Source_Method)
  supervised <- isTRUE(provenance == "spls")
  inference <- if(supervised) "exploratory_supervised_score" else
    if(provenance %in% c("single","sum","pc1","pc2","pc3")) "unsupervised_score" else "unknown_weight_provenance"
  if(supervised) warning("Supervised microbial weights: nominal association p-values do not account for learning the score. Treat as exploratory unless trained independently.",call.=FALSE)
  if(inference=="unknown_weight_provenance") warning("Saved weights have unknown training provenance; verify independence from the measured outcome before confirmatory interpretation.",call.=FALSE)
  outdir <- file.path(args$wdir,paste0(args$projname,"_integrated_gwas"))
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  complete <- FALSE
  run_status <- function(status) .sgwas_write_tsv(data.frame(Status=status,Mode=mode,Time=as.character(Sys.time())),
                                                file.path(outdir,"RUN_STATUS.tsv"))
  run_status("running")
  on.exit(if(!complete) run_status("failed"),add=TRUE)
  results <- list(); summaries <- list(); audit <- list(); allfits <- list(); outputs <- list(); idx <- 0L
  score <- prepared$phenotype[,1:2,drop=FALSE]; micro_name <- names(score)[2]
  for(ploidy in unique(as.numeric(args$ploidy_levels))) {
    raw <- .sgwas_read_dosage(args[[paste0("genofile_",ploidy,"x")]],ploidy)
    for(trait in traits) {
      idx <- idx+1L; label <- paste0("pair",idx,"_",ploidy,"x")
      folder <- file.path(outdir,label); dir.create(folder,showWarnings=FALSE)
      ids <- score$Plant_ID
      y <- .sgwas_id_values(ph,trait,ids)
      C <- matrix(numeric(),nrow=length(ids),ncol=0L,dimnames=list(NULL,character()))
      if(length(covs)) C <- vapply(covs,function(nm) .sgwas_id_values(ph,nm,ids),numeric(length(ids)))
      colnames(C) <- covs
      keep <- is.finite(y) & is.finite(score[[2]]) & rowSums(!is.finite(C))==0L & ids %in% rownames(raw$G)
      audit[[idx]] <- data.frame(Pair=label,Trait=trait,Microbial_Trait=micro_name,Ploidy=ploidy,
        Plant_ID=ids,Observed_Trait=is.finite(y),Complete_Covariates=rowSums(!is.finite(C))==0L,
        Has_Genotype=ids %in% rownames(raw$G),Included=keep)
      ids <- ids[keep]; Y <- cbind(Trait=y[keep],Microbial=score[[2]][keep]); C <- C[keep,,drop=FALSE]
      if(length(ids)<max(10L,ncol(C)+n_pcs+5L)) stop("Too few complete samples for ",trait)
      if(any(!is.finite(apply(Y,2,sd))) || any(apply(Y,2,sd)<=0)) stop("Both outcomes must vary on the analysed cohort")
      if(ncol(C)) {
        cs <- apply(C,2,sd)
        if(any(!is.finite(cs)) || any(cs<=0)) stop("Constant or invalid phenotype covariate on the analysed cohort")
        C <- scale(C)
      }
      qc <- .sgwas_dosage_qc(raw,ids,ploidy,args$maf,max_missing)
      .sgwas_write_tsv(qc$ledger,file.path(folder,"marker_qc.tsv"))
      cp <- tcrossprod(qc$centered); full_eig <- .sgwas_kinship_eigen(cp)
      X <- cbind(Intercept=1,C)
      if(n_pcs) {
        if(n_pcs>sum(full_eig$values>max(full_eig$values)*1e-8)) stop("Too few nonzero genomic PCs")
        pcs <- full_eig$vectors[,seq_len(n_pcs),drop=FALSE]
        colnames(pcs) <- paste0("Genotype_PC",seq_len(n_pcs)); X <- cbind(X,pcs)
      }
      if(qr(X)$rank<ncol(X)) stop("Collinear covariates/PCs in integrated GWAS for ",trait)
      .sgwas_write_tsv(data.frame(Plant_ID=ids,Y,Microbial_standardised=as.numeric(scale(Y[,2])),X,check.names=FALSE),
                       file.path(folder,"analysed_samples.tsv"))
      groups <- if(args$LOCO) unique(qc$map$CHROM) else "All"
      if(args$LOCO && length(groups)<2L) stop("LOCO requires retained markers on at least two chromosomes")
      result_chr <- list(); fit_chr <- list()
      for(ch in groups) {
        selected <- if(args$LOCO) which(qc$map$CHROM==ch) else seq_len(ncol(qc$G))
        eg <- if(args$LOCO) .sgwas_kinship_eigen(cp-tcrossprod(qc$centered[,selected,drop=FALSE])) else full_eig
        message("[systemsGWAS] ",trait," + ",micro_name,"; ",ploidy,"x; ",length(ids)," samples; ",
                length(selected)," markers; ",if(args$LOCO) paste0("LOCO ",ch) else "whole-genome kinship")
        fit <- tryCatch(.sgwas_association_pair(Y,X,qc$centered[,selected,drop=FALSE],qc$map[selected,,drop=FALSE],
                        eg,mode,chunk_size,maxit),error=function(e) {
          .sgwas_write_tsv(data.frame(Trait=trait,Chromosome=ch,Error=conditionMessage(e)),file.path(folder,"FAILED.tsv"))
          stop(e)
        })
        r <- cbind(Pair=label,Trait=trait,Microbial_Trait=micro_name,Ploidy=ploidy,N=length(ids),
                   Inference=inference,fit$results)
        result_chr[[ch]] <- r; fit_chr[[ch]] <- fit$fits
        for(nm in names(fit$fits)) {
          f <- fit$fits[[nm]]; gg <- if(nm=="Bivariate") f$G_original else f$G
          rr <- if(nm=="Bivariate") f$R_original else f$R
          summaries[[length(summaries)+1L]] <- data.frame(Pair=label,Chromosome=ch,Model=nm,
            N=length(ids),Converged=f$converged,Boundary=f$boundary,Optimizer=f$optimizer,
            Objective=f$objective,Genetic_Var1=gg[1,1],Residual_Var1=rr[1,1],
            Genetic_Var2=if(nrow(gg)==2) gg[2,2] else NA_real_,
            Residual_Var2=if(nrow(rr)==2) rr[2,2] else NA_real_,
            Genetic_Cov=if(nrow(gg)==2) gg[1,2] else NA_real_,
            Residual_Cov=if(nrow(rr)==2) rr[1,2] else NA_real_,
            Genetic_Correlation=if(nrow(gg)==2 && !f$boundary) cov2cor(gg)[1,2] else NA_real_,
            Residual_Correlation=if(nrow(rr)==2) cov2cor(rr)[1,2] else NA_real_)
        }
      }
      results[[idx]] <- do.call(rbind,result_chr); allfits[[label]] <- fit_chr
      outputs[[idx]] <- folder
      if(file.exists(file.path(folder,"FAILED.tsv"))) unlink(file.path(folder,"FAILED.tsv"))
    }
  }
  results <- do.call(rbind,results); rownames(results) <- NULL
  pcols <- names(results)[startsWith(names(results),"P_")]
  total_tests <- setNames(integer(length(pcols)),pcols)
  for(p in pcols) {
    good <- is.finite(results[[p]]); total_tests[[p]] <- sum(good)
    for(adj in c("BH","bonferroni")) {
      value <- rep(NA_real_,nrow(results)); value[good] <- p.adjust(results[[p]][good],method=adj)
      results[[paste0(p,if(adj=="BH") "_BH" else "_Bonferroni")]] <- value
    }
  }
  .sgwas_write_tsv(results,file.path(outdir,"ASSOCIATIONS.tsv"))
  summary <- do.call(rbind,summaries); .sgwas_write_tsv(summary,file.path(outdir,"NULL_MODELS.tsv"))
  cohort <- do.call(rbind,audit); .sgwas_write_tsv(cohort,file.path(outdir,"SAMPLE_AUDIT.tsv"))
  for(i in seq_along(outputs)) {
    r <- results[results$Pair==unique(results$Pair)[i],,drop=FALSE]
    .sgwas_write_tsv(r,file.path(outputs[[i]],"associations.tsv"))
    .sgwas_association_plots(r,file.path(outputs[[i]],"Manhattan_QQ.pdf"),args$fdr,args$bonferroni,args$suggestive,total_tests)
  }
  metadata <- list(mode=mode,arguments=args,n_pcs=n_pcs,max_missing=max_missing,chunk_size=chunk_size,maxit=maxit,
    inference=inference,microbiome=prepared$model,models=allfits,tests=total_tests,session=sessionInfo(),
    method="Additive REML null/P3D covariance shape; per-marker GLS residual scale and approximate F tests",
    multiple_testing="Across markers, trait pairs and ploidies in this call, separately for each P column")
  saveRDS(metadata,file.path(outdir,"MODEL_DETAILS.rds"))
  complete <- TRUE; run_status("complete")
  message("[systemsGWAS] Integrated association results: ",file.path(outdir,"ASSOCIATIONS.tsv"))
  invisible(list(results=results,null_models=summary,sample_audit=cohort,metadata=metadata,outdir=outdir))
}

# Public entry point. Legacy named arguments remain accepted; prefer
# microbial_features + microbial_method for new microbial-phenotype analyses.
systemsGWAS <- function(
    wdir="./", projname="GWAS", ploidy_levels=c(2,4,6,8),
    trait_names=c("trait1","trait2"), alternate_trait=NULL,
    trait_microbial_proxy=c("NULL"), model_effect=c("Add","Dom"),
    fdr=TRUE, bonferroni=TRUE, suggestive=5, perm=1, cores=1,
    genofile_2x=NULL, genofile_4x=NULL, genofile_6x=NULL, genofile_8x=NULL,
    phenofile=NULL, method=c("GLM","MLM"), covariate_pheno=NULL,
    covariate_metag=FALSE, maf=0.02, LOCO=FALSE,
    pheno_taxa_strain=NULL, pheno_taxa_species=NULL,
    metag_data_strains=NULL, metag_data_species=NULL,
    microbial_features=NULL, microbial_method="auto", microbial_response=NULL,
    microbial_weights=NULL, microbial_trait_name="Microbial_composite",
    microbial_transform="auto", microbial_prevalence=0.05,
    microbial_bootstraps=100L, microbial_stability=0.60,
    microbial_seed=123L, microbial_pseudocount=1e-6, run_gwas=TRUE,
    association_mode="standard", association_n_pcs=0L,
    association_max_missing=0.20, association_chunk_size=256L,
    association_maxit=300L) {
  association_mode <- match.arg(tolower(association_mode),c("standard","joint","interaction","both"))
  integrated <- association_mode!="standard"
  if(integrated) {
    if(missing(method)) method <- "MLM"
    if(missing(model_effect)) model_effect <- "Add"
    association_n_pcs <- .sgwas_integer(association_n_pcs,"association_n_pcs",0L)
    association_chunk_size <- .sgwas_integer(association_chunk_size,"association_chunk_size")
    association_maxit <- .sgwas_integer(association_maxit,"association_maxit")
    if(length(association_max_missing)!=1L || !is.numeric(association_max_missing) ||
       !is.finite(association_max_missing) || association_max_missing<0 || association_max_missing>1)
      stop("association_max_missing must be between 0 and 1")
  }
  if(length(run_gwas)!=1L || !is.logical(run_gwas) || is.na(run_gwas))
    stop("run_gwas must be TRUE or FALSE")
  wdir <- normalizePath(wdir,mustWork=TRUE)
  if(length(projname)!=1L || is.na(projname) || !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$",projname))
    stop("projname must be a simple project prefix (letters, digits, _, . or -)")
  numeric_cores <- suppressWarnings(as.numeric(cores))
  if(length(numeric_cores)!=1L || !is.finite(numeric_cores) || numeric_cores<1 || numeric_cores!=floor(numeric_cores))
    stop("cores must be a positive integer")
  new_mode <- length(.sgwas_names(microbial_features))>0L || length(.sgwas_names(microbial_weights))>0L
  if(integrated && !new_mode) stop("Integrated GWAS requires microbial_features or microbial_weights")
  legacy_names <- names(formals(.systemsGWAS_engine))
  args <- mget(legacy_names,envir=environment(),inherits=FALSE)
  args$cores <- as.integer(numeric_cores)
  prepared <- NULL
  if(new_mode) {
    if(length(.sgwas_names(trait_microbial_proxy)) || length(.sgwas_names(alternate_trait)) ||
       length(.sgwas_names(pheno_taxa_species)) || length(.sgwas_names(pheno_taxa_strain)))
      stop("Use microbial_features/microbial_method without the legacy proxy or pheno_taxa selectors")
    if(!identical(covariate_metag,FALSE))
      stop("Use covariate_metag = FALSE when mapping a constructed microbial phenotype")
    if(is.null(.sgwas_names(metag_data_species))) stop("Supply metag_data_species for microbial phenotype construction")
    # FER is read only for sPLS, or if a phenotype-table covariate is requested.
    needs_pheno <- identical(tolower(microbial_method),"spls") || length(.sgwas_names(covariate_pheno))>0L
    outdir <- file.path(wdir,paste0(projname,"_microbial_traits"))
    prepared <- systemsGWAS_prepare_microbiome(
      microbiome_file=.sgwas_path(metag_data_species,wdir), features=microbial_features,
      method=microbial_method, response=microbial_response,
      phenofile=if(needs_pheno) .sgwas_path(phenofile,wdir) else NULL,
      weights_file=.sgwas_path(microbial_weights,wdir), trait_name=microbial_trait_name,
      transform=microbial_transform, prevalence=microbial_prevalence,
      bootstraps=microbial_bootstraps, stability=microbial_stability,
      seed=microbial_seed, pseudocount=microbial_pseudocount,
      outdir=outdir, covariates=covariate_pheno)
    if(!run_gwas) return(invisible(prepared))
    if(!integrated) {
      args$trait_names <- microbial_trait_name
      args$phenofile <- normalizePath(prepared$paths[['phenotype']],mustWork=TRUE)
      args[c('alternate_trait','trait_microbial_proxy','pheno_taxa_strain',
             'pheno_taxa_species','metag_data_strains','metag_data_species')] <- rep(list(NULL),6L)
      args$covariate_metag <- FALSE
    }
  } else {
    if(!run_gwas) stop("run_gwas = FALSE requires microbial_features or microbial_weights")
    if(!identical(microbial_method,"auto") || length(.sgwas_names(microbial_response)))
      stop("Supply microbial_features (or microbial_weights) to activate microbial phenotype construction")
  }
  # Check genotype/phenotype inputs before the engine loads its GWAS packages.
  ploidies <- suppressWarnings(as.numeric(ploidy_levels))
  if(!length(ploidies) || anyNA(ploidies) || any(!ploidies %in% c(2,4,6,8)))
    stop("ploidy_levels must contain 2, 4, 6 or 8")
  for(p in ploidies) {
    nm <- paste0('genofile_',p,'x')
    if(is.null(.sgwas_names(args[[nm]]))) stop("Supply ",nm," for ploidy ",p)
    args[[nm]] <- .sgwas_path(args[[nm]],wdir)
  }
  for(nm in c('phenofile','metag_data_species','metag_data_strains'))
    args[[nm]] <- .sgwas_path(args[[nm]],wdir)
  gwas_result <- if(integrated) .sgwas_integrated_gwas(args,prepared,association_mode,
    association_n_pcs,association_max_missing,association_chunk_size,association_maxit) else
    do.call(.systemsGWAS_engine,args)
  if(!is.null(prepared)) {prepared$gwas_result <- gwas_result;return(invisible(prepared))}
  invisible(gwas_result)
}
