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
  if (any(!is.finite(X)))
    stop("Microbiome values must be finite numbers (use 0 or numeric LSMeans, not NA/Inf)")
  min_val <- min(X)
  if (min_val < 0) {
    message("[systemsGWAS] Negative microbiome values detected (min = ", round(min_val, 4),
            "); shifting all values up by ", round(abs(min_val), 4), " so minimum is 0.")
    X <- X - min_val
  }
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

.sgwas_transform <- function(X, method, pseudocount) {
  if (method == "none") return(X)
  if (method == "log1p") return(log1p(X))
  total <- rowSums(X)
  if (any(total <= 0)) {
    # Shift by pseudocount if any sample has zero total after min-shifting
    X <- X + pseudocount
    total <- rowSums(X)
  }
  relative <- X / total
  if (method == "relative") return(relative)
  if (method == "clr") {
    # Close to proportions before adding the pseudocount. Reference columns
    # are fixed in the exported weights, including zero-weight reference taxa.
    logx <- log(relative + pseudocount)
    return(logx - rowMeans(logx))
  }
  stop("Unknown microbial transformation: ", method)
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
    prevalence=0.05, bootstraps=100L, stability=0.30, seed=123L,
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

.sgwas_marginal_scan <- function(y, X, G, map, eig, chunk_size,
                                 model_effect = "Add", ploidy = 2L, dosage = NULL, max_geno_freq = 0.99) {
  models <- .sgwas_marker_models(model_effect, ploidy)
  if (any(models != "additive") && (is.null(dosage) || !identical(dim(dosage), dim(G))))
    stop("Dominance scans require dosages matching G")
  
  var_y <- var(y, na.rm = TRUE)
  calc_pve <- function(beta, var_x) {
    if (!is.finite(var_y) || var_y <= 0) return(rep(NA_real_, length(beta)))
    pve <- (beta^2 * var_x) / var_y
    pve[!is.finite(pve)] <- NA_real_
    pmin(pve, 1)
  }
  
  fits <- list(Trait = .sgwas_fit_univariate(y, X, eig))
  chunks <- ceiling(ncol(G) / chunk_size)
  res <- vector("list", chunks * length(models)); result_index <- 0L
  
  for (model in models) for (k in seq_len(chunks)) {
    js <- seq.int((k - 1L) * chunk_size + 1L, min(k * chunk_size, ncol(G)))
    coding <- if (model == "additive")
      list(G = G[, js, drop = FALSE], n0 = rep(NA_integer_, length(js)), n1 = rep(NA_integer_, length(js)),
           status = rep("ok", length(js))) else
             .sgwas_encode_dominance(dosage[, js, drop = FALSE], model, ploidy, max_geno_freq)
    g <- coding$G
    var_x <- apply(g, 2, var, na.rm = TRUE)
    gr <- crossprod(eig$vectors, g)
    a <- .sgwas_marker_test(fits$Trait$fit, gr)
    r <- cbind(map[js, , drop = FALSE], Marker_Model = model, N_Group0 = coding$n0, N_Group1 = coding$n1,
               Marker_Coding_Status = coding$status,
               Beta_Trait = a$Beta1, SE_Trait = a$SE1, PVE_Trait = calc_pve(a$Beta1, var_x),
               P_Trait = a$P1, Status_FER = a$Status)
    bad <- coding$status != "ok"
    if (any(bad)) {
      r[bad, grepl("^(Beta_|SE_|PVE_|P_)", names(r))] <- NA_real_
      r[bad, startsWith(names(r), "Status_")] <- coding$status[bad]
    }
    result_index <- result_index + 1L; res[[result_index]] <- r
  }
  list(results = do.call(rbind, res), fits = fits)
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
    fdr_level = 0.05,
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
    print("Please provide <trait names> or <taxa names>. You can't provide both")
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
            
            data <- .sgwas_read_gwaspoly(ploidy,pheno,geno)
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
              
              
              ps$Model <- factor(ps$model, levels = score_models)
              qqplot <- ggplot(ps, aes(x = expected, y = observed)) +
                geom_ribbon(aes(ymin = clower, ymax = cupper), alpha = 0.12, fill = "grey60") +
                geom_point(size = 2.5, colour = "cornflowerblue") +
                geom_abline(intercept = 0, slope = 1, colour = "red", linewidth = 0.8, linetype = "dashed") +
                facet_wrap(~Model, ncol = 1) +
                theme_bw(base_size = 16) +
                theme(panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                      legend.box.spacing = ggplot2::unit(0.2, "cm"), 
                      legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                      strip.text = ggplot2::element_text(size = 18, colour = "black"),
                      plot.title = ggplot2::element_text(size = 15, hjust = 0.5),
                      axis.title = ggplot2::element_text(size = 18),
                      axis.text = ggplot2::element_text(size = 14)) +
                xlab(log10Pe) + ylab(log10Po) + labs(title = traitname)
              ggsave(file = paste("./qqplots/", "QQplot_", traitname, "_fdr0.05.tiff", sep = ""),
                     plot = qqplot, width = 5, height = max(4, 2 * length(score_models) + 1),
                     units = "in", dpi = 300, compression = "lzw", limitsize = FALSE)
              
              
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
              scores_mat <- as.matrix(GWAS.fitted@scores[[traitname]])
              colnames(scores_mat) <- gsub("^scores\\.", "", colnames(scores_mat))
              marker_ids <- if (!is.null(rownames(scores_mat))) rownames(scores_mat) else GWAS.fitted@map[["Marker"]]
              fdr_thresh_vec <- setNames(rep(Inf, ncol(scores_mat)), colnames(scores_mat))
              sig_fdr_list <- list()
              for (m in colnames(scores_mat)) {
                sc <- as.numeric(scores_mat[, m])
                qv <- .sgwas_fdr(10^(-sc))
                sig <- which(is.finite(qv) & qv <= fdr_level)
                if (length(sig)) {
                  fdr_thresh_vec[m] <- min(sc[sig], na.rm = TRUE)
                  sig_fdr_list[[m]] <- data.frame(
                    Marker = marker_ids[sig],
                    Model = m,
                    Score = sc[sig],
                    FDR = qv[sig],
                    stringsAsFactors = FALSE
                  )
                }
              }
              
              SigQTL_fdr <- if (length(sig_fdr_list)) {
                out_fdr <- do.call(rbind, sig_fdr_list)
                merge(out_fdr, GWAS_scores_effects_long, by = c("Marker", "Model"))
              } else {
                data.frame()
              }
              write.table(SigQTL_fdr, paste("./sigFDR/", "Significant_effect_", traitname, "_fdr", fdr_level, ".txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
              
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
              
              hline.fdr <- data.frame(z1 = as.numeric(fdr_thresh_vec), models = names(fdr_thresh_vec))
              hline.fdr <- hline.fdr[is.finite(hline.fdr$z1), , drop = FALSE]
              
              Bonferroni_threshold <- t(data_Bonferroni@threshold)
              Bonferroni_threshold <- setDT(as.data.frame(Bonferroni_threshold), keep.rownames = TRUE)
              colnames(Bonferroni_threshold) <- c("model","threshold")
              hline.Bonferroni <- data.frame(z2 = c(Bonferroni_threshold$threshold), models = c(Bonferroni_threshold$model))
              
              scores$Chrom <- gsub("^(chr|chrom)_?", "", scores$Chrom, ignore.case = TRUE)
              chr_levels <- unique(as.character(scores$Chrom))
              chr_num <- suppressWarnings(as.numeric(chr_levels))
              chr_levels <- chr_levels[order(is.na(chr_num), chr_num, chr_levels)]
              scores$Chrom <- factor(as.character(scores$Chrom), levels = chr_levels, ordered = TRUE)
              
              scores$models <- gsub("\\.", "-", as.character(scores$models))
              score_models_clean <- gsub("\\.", "-", score_models)
              scores$models <- factor(scores$models, levels = score_models_clean)
              scores <- scores[!is.na(scores$models), ]
              scores <- scores[order(scores$models, scores$Chrom, scores$bp), ]
              h_list <- list()
              if (threshold_FDR && nrow(hline.fdr) > 0)
                h_list[[length(h_list) + 1L]] <- data.frame(models = gsub("\\.", "-", hline.fdr$models), Threshold = "FDR", Value = hline.fdr$z1)
              if (threshold_Bonferroni && nrow(hline.Bonferroni) > 0)
                h_list[[length(h_list) + 1L]] <- data.frame(models = gsub("\\.", "-", hline.Bonferroni$models), Threshold = "Bonferroni", Value = hline.Bonferroni$z2)
              if (!is.null(threshold_suggestive) && is.finite(threshold_suggestive))
                h_list[[length(h_list) + 1L]] <- data.frame(models = score_models_clean, Threshold = "Suggestive", Value = threshold_suggestive)
              h <- if (length(h_list)) do.call(rbind, h_list) else data.frame(models = character(), Threshold = character(), Value = numeric())
              h <- h[is.finite(h$Value), , drop = FALSE]
              
              chrom_cols <- rep(c("cornflowerblue", "orange"), length.out = length(chr_levels))
              names(chrom_cols) <- chr_levels
              palette <- c(chrom_cols, FDR = "grey10", Bonferroni = "tomato", Suggestive = "cornflowerblue")
              
              base <- ggplot(scores, aes(x = bp, y = scores, colour = Chrom)) +
                facet_grid(models ~ Chrom, scales = "free_x", space = "free_x") +
                scale_colour_manual(values = palette, breaks = h$Threshold) +
                scale_linetype_manual(values = c(FDR = "dashed", Bonferroni = "dashed", Suggestive = "dashed"), breaks = h$Threshold) +
                scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
                theme_bw(base_size = 18) +
                theme(panel.spacing.x = grid::unit(0, "lines"), panel.spacing.y = grid::unit(0.5, "lines"),
                      panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_blank(),
                      panel.border = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                      axis.ticks.x = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "white"),
                      axis.title.x = ggplot2::element_text(size = 15, margin = margin(t = 12)),
                      axis.title.y = ggplot2::element_text(size = 15, margin = margin(r = 12)),
                      axis.text.y = ggplot2::element_text(size = 15, colour = "black"), strip.background = ggplot2::element_blank(),
                      strip.text.x = ggplot2::element_text(size = 15, colour = "black"), strip.text.y = ggplot2::element_text(size = 12, colour = "black"),
                      legend.position = "top", legend.justification = "right",
                      legend.box.just = "right", legend.background = ggplot2::element_rect(fill = "white", colour = "white"),
                      legend.box.spacing = ggplot2::unit(0.2, "cm"), 
                      legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                      legend.title = ggplot2::element_text(size = 17), legend.text = ggplot2::element_text(size = 15),
                      legend.key.width = grid::unit(1.5, "cm"), legend.key.height = grid::unit(0.7, "cm"),
                      plot.title = ggplot2::element_text(size = 5, hjust = 0.5, face = "plain", margin = margin(b = 15)),
                      plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
                labs(title = traitname, x = "Chromosomes and Physical Map Position (bp)",
                     y = "-log10(P)", colour = "Threshold", linetype = "Threshold") +
                coord_cartesian(clip = "off")
              
              plot_height <- max(4.5, 1.8 * length(score_models_clean) + 1.5)
              
              if (diversity == TRUE) {
                manplot <- base + geom_point(size = 2, shape = 20, show.legend = FALSE) +
                  geom_hline(data = h, aes(yintercept = Value, colour = Threshold, linetype = Threshold), linewidth = 0.35, inherit.aes = FALSE) +
                  scale_colour_manual(values = palette, breaks = h$Threshold) +
                  scale_linetype_manual(values = c(FDR = "dashed", Bonferroni = "dashed", Suggestive = "dashed"), breaks = h$Threshold)
                ggsave(file = paste("./manplots/", "manplot_", traitname, ".tiff", sep = ""),
                       plot = manplot, width = 12, height = plot_height, units = "in", dpi = 300, compression = "lzw", limitsize = FALSE)
              }
              
              if (biparental == TRUE) {
                scores$Smoothed <- scores$scores
                for (ix in split(seq_len(nrow(scores)), interaction(scores$models, scores$Chrom, drop = TRUE))) {
                  if (length(ix) > 1L) {
                    scores$Smoothed[ix] <- c(NA_real_, (head(scores$scores[ix], -1) + tail(scores$scores[ix], -1)) / 2)
                  }
                }
                qtlplot <- base +
                  geom_line(data = scores, aes(y = Smoothed, group = interaction(models, Chrom)),
                            linewidth = 1.0, show.legend = FALSE, na.rm = TRUE) +
                  geom_hline(data = h, aes(yintercept = Value, colour = Threshold, linetype = Threshold), linewidth = 0.35, inherit.aes = FALSE) +
                  scale_colour_manual(values = palette, breaks = h$Threshold) +
                  scale_linetype_manual(values = c(FDR = "dashed", Bonferroni = "dashed", Suggestive = "dashed"), breaks = h$Threshold)
                ggsave(file = paste("./manplots/", "QTLprofile_", traitname, ".tiff", sep = ""),
                       plot = qtlplot, width = 12, height = plot_height, units = "in", dpi = 300, compression = "lzw", limitsize = FALSE)
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

# Keep response-training input separate from the constructed GWAS phenotype.
.sgwas_training_phenofile <- function(path, base, spls=FALSE) {
  tryCatch(.sgwas_path(path,base),error=function(e) {
    hint <- if(spls) paste0("For microbial_method = 'spls', phenofile must be the original measured-trait TSV ",
                            "containing microbial_response (e.g. FER). Do not supply the generated composite phenotype as its training file. ") else ""
    stop(hint,"Cannot open phenofile: ",conditionMessage(e),call.=FALSE)
  })
}

.sgwas_standard_proxy_args <- function(args, prepared) {
  file <- normalizePath(prepared$paths[['phenotype']],mustWork=TRUE)
  trait <- as.character(prepared$summary$Trait)
  ph <- .sgwas_phenotypes(file)
  required <- c(trait,.sgwas_names(args$covariate_pheno))
  if(any(!required %in% names(ph)))
    stop("Generated GWAS phenotype file is missing: ",paste(setdiff(required,names(ph)),collapse=", "))
  if(trait %in% .sgwas_names(args$covariate_pheno)) stop("A constructed GWAS outcome cannot be its own covariate")
  training_file <- args$phenofile
  args$trait_names <- trait
  args$phenofile <- file
  args[c('alternate_trait','trait_microbial_proxy','pheno_taxa_strain',
         'pheno_taxa_species','metag_data_strains','metag_data_species')] <- rep(list(NULL),6L)
  args$covariate_metag <- FALSE
  .sgwas_write_tsv(data.frame(Association_Mode="gwaspoly-marginal",GWAS_Trait=trait,GWAS_Phenotype_File=file,
                              Original_Phenofile=if(is.null(training_file)) NA_character_ else training_file,
                              N_Scores=sum(is.finite(ph[[trait]]))),
                   file.path(dirname(file),paste0(trait,"_gwas_input.tsv")))
  message("[systemsGWAS] Standard proxy GWAS outcome: ",trait,"; phenotype file: ",file)
  args
}

.sgwas_read_gwaspoly <- function(ploidy, phenotypes, genotypes) {
  if(!requireNamespace("GWASpoly",quietly=TRUE)) stop("Install GWASpoly for gwaspoly-marginal GWAS")
  map_columns <- if(ncol(genotypes)>=5L && identical(toupper(names(genotypes)[4:5]),c("REF","ALT"))) 5L else 3L
  columns <- seq.int(map_columns+1L,ncol(genotypes)); ids <- names(genotypes)[columns]
  if(anyDuplicated(ids) || any(!nzchar(ids))) stop("Genotype sample IDs must be unique and nonempty")
  ph_ids <- as.character(phenotypes[[1]])
  selected <- match(ph_ids,ids)
  if(anyNA(selected)) stop("Some phenotype IDs are absent from the genotype table passed to GWASpoly")
  # GWASpoly read.table() type-converts its phenotype IDs (0001 becomes 1).
  # Alias both inputs during import and restore exact original IDs immediately.
  aliases <- paste0("sgwas_sample_",seq_along(ids))
  names(genotypes)[columns] <- aliases
  phenotypes[[1]] <- aliases[selected]
  directory <- tempfile("systemsGWAS_native_import_"); dir.create(directory)
  on.exit(unlink(directory,recursive=TRUE),add=TRUE)
  phfile <- file.path(directory,"phenotypes.tsv"); genofile <- file.path(directory,"genotypes.tsv")
  write.table(phenotypes,phfile,sep="\t",row.names=FALSE,quote=TRUE,na="NA")
  write.table(genotypes,genofile,sep="\t",row.names=FALSE,quote=TRUE,na="NA")
  data <- GWASpoly::read.GWASpoly(ploidy,phfile,genofile,"numeric",1,delim="\t")
  data@pheno[[1]] <- ids[match(as.character(data@pheno[[1]]),aliases)]
  rownames(data@geno) <- ids[match(rownames(data@geno),aliases)]
  if(anyNA(data@pheno[[1]]) || anyNA(rownames(data@geno))) stop("GWASpoly import lost sample ID alignment")
  data
}

# Integrated association models: additive/dominance, complete-case cohorts,
# REML null covariance estimates and P3D/GLS marker tests. Base R only.
.sgwas_write_tsv <- function(x, file) {
  write.table(x, file, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
}

.sgwas_fdr <- function(p) {
  q <- rep(NA_real_, length(p))
  good <- is.finite(p) & p >= 0 & p <= 1
  if (!any(good)) return(q)
  pv <- p[good]
  if (requireNamespace("qvalue", quietly = TRUE) && length(pv) >= 20L) {
    res <- tryCatch(qvalue::qvalue(pv)$qvalues, error = function(e) NULL)
    if (!is.null(res) && length(res) == length(pv) && all(is.finite(res))) {
      q[good] <- pmin(1, res)
      return(q)
    }
  }
  q[good] <- p.adjust(pv, method = "BH")
  q
}

.sgwas_integer <- function(x, name, minimum=1L) {
  if(length(x)!=1L || !is.numeric(x) || !is.finite(x) || x<minimum || x!=floor(x))
    stop(name," must be an integer >= ",minimum)
  as.integer(x)
}

.sgwas_validate_integrated_models <- function(method, model_effect) {
  if(!identical(.sgwas_names(method),"MLM"))
    stop('The sgwas engine requires method = "MLM"')
  effects <- .sgwas_names(model_effect)
  if(!length(effects) || any(!effects %in% c("Add","Dom")))
    stop('Integrated model_effect must contain "Add", "Dom", or both')
  invisible(effects)
}

.sgwas_marker_models <- function(model_effect, ploidy) {
  effects <- .sgwas_validate_integrated_models("MLM",model_effect)
  models <- if("Add" %in% effects) "additive" else character()
  if("Dom" %in% effects) {
    # Same allele-threshold families as the legacy GWASpoly scans.
    for(d in seq_len(ploidy %/% 2L))
      models <- c(models,paste0(d,"-dom-alt"),paste0(d,"-dom-ref"))
  }
  models
}

.sgwas_encode_dominance <- function(dosage, model, ploidy, max_geno_freq = 0.99) {
  parts <- strsplit(model, "-", fixed = TRUE)[[1]]
  d <- as.integer(parts[1])
  calls <- round(dosage)
  B <- if (parts[3] == "alt") 1 * (calls >= d) else 1 * (calls > ploidy - d)
  
  n0 <- colSums(B == 0, na.rm = TRUE)
  n1 <- colSums(B == 1, na.rm = TRUE)
  n_obs <- n0 + n1
  max_dom_gf <- pmax(n0, n1) / pmax(1L, n_obs)
  
  mu <- colMeans(B, na.rm = TRUE); mu[!is.finite(mu)] <- 0
  for (j in which(colSums(is.na(B)) > 0)) B[is.na(B[, j]), j] <- mu[j]
  
  list(
    G = sweep(B, 2, mu, "-"),
    n0 = n0,
    n1 = n1,
    status = ifelse(max_dom_gf > max_geno_freq | pmin(n0, n1) == 0L, "insufficient_dominance_group", "ok")
  )
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

.sgwas_dosage_qc <- function(raw, ids, ploidy, maf, max_missing, max_geno_freq = 0.99) {
  G <- raw$G[ids, , drop = FALSE]
  missing <- colMeans(is.na(G))
  af <- colMeans(G, na.rm = TRUE) / ploidy
  marker_maf <- pmin(af, 1 - af)
  
  # Calculate maximum genotype class frequency (GWASpoly geno.freq filter)
  calls <- round(G)
  max_gf <- vapply(seq_len(ncol(calls)), function(j) {
    v <- calls[, j]
    v <- v[is.finite(v)]
    if (!length(v)) return(1)
    max(table(v)) / length(v)
  }, numeric(1))
  
  variable <- vapply(seq_len(ncol(G)), function(j) {
    v <- G[, j]; v <- v[is.finite(v)]; length(v) > 1L && diff(range(v)) > 1e-10
  }, logical(1))
  
  keep <- is.finite(marker_maf) & marker_maf >= maf & missing <= max_missing &
    max_gf <= max_geno_freq & variable
  
  ledger <- cbind(raw$map, MAF = marker_maf, Max_Geno_Freq = max_gf, Missing_Rate = missing,
                  N_Observed = colSums(!is.na(G)),
                  N_Het_calls = colSums(calls > 0 & calls < ploidy, na.rm = TRUE),
                  N_Fractional = colSums(abs(G - calls) > 1e-8, na.rm = TRUE),
                  Variable = variable, Retained = keep)
  if (!any(keep)) stop("No variable markers pass MAF, geno.freq, and missingness filters")
  
  G <- G[, keep, drop = FALSE]
  # Mode imputation to match GWASpoly on integer dosage data
  for (j in which(colSums(is.na(G)) > 0)) {
    obs <- G[!is.na(G[, j]), j]
    G[is.na(G[, j]), j] <- as.numeric(names(which.max(table(obs))))
  }
  mu <- colMeans(G, na.rm = TRUE)
  list(G = G, centered = sweep(G, 2, mu, "-"), map = ledger[keep, , drop = FALSE], ledger = ledger)
}

.sgwas_kinship_eigen <- function(K_raw, denom) {
  if (!is.finite(denom) || denom <= 1e-12) stop("No variable background markers are available for kinship")
  K <- ((K_raw + t(K_raw)) / 2) / denom
  
  # Match legacy GWASpoly pipeline normalize_kinmat()
  K <- K - min(K)
  max_k <- max(K)
  if (!is.finite(max_k) || max_k <= 1e-12) stop("Kinship matrix is constant after normalization")
  K <- K / max_k
  
  eg_init <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  diag(K) <- diag(K) - min(eg_init)
  K <- K / mean(diag(K))   # <-- Matches GWASpoly::set.K() internal rescaling!
  
  eg <- eigen(K, symmetric = TRUE)
  if (min(eg$values) < -1e-7 * max(eg$values)) stop("Kinship is not positive semidefinite")
  eg$values <- pmax(eg$values, 0)
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
  yr <- crossprod(U, y); xr <- crossprod(U, X); df <- length(y) - ncol(X)
  profile <- function(h, details = FALSE) {
    v <- (1 - h) + h * lambda
    fit <- .sgwas_gls_null(yr, xr, list(1 / v))
    s <- fit$rss / df
    value <- .5 * (df * (log(2 * pi) + 1 + log(s)) + sum(log(v)) + fit$logdetA)
    if (details) return(list(fit = fit, variance = s, h = h, objective = value))
    value
  }
  upper <- 1 - 1e-6
  grid_h <- seq(0, upper, length.out = 25L)
  grid_val <- vapply(grid_h, profile, numeric(1))
  best_i <- which.min(grid_val)
  lo <- grid_h[max(1L, best_i - 1L)]
  hi <- grid_h[min(length(grid_h), best_i + 1L)]
  
  opt <- optimize(profile, c(lo, hi), tol = 1e-8)
  candidates <- unique(c(0, grid_h[best_i], opt$minimum, upper))
  values <- vapply(candidates, profile, numeric(1))
  best <- profile(candidates[which.min(values)], TRUE)
  
  best$G <- matrix(best$variance * best$h, 1, 1)
  best$R <- matrix(best$variance * (1 - best$h), 1, 1)
  best$boundary <- best$h < 1e-5 || best$h > upper - 1e-5
  best$converged <- TRUE; best$optimizer <- "grid_bracketed_Brent_REML"
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
# Under EMMAX / P3D, both the covariance structure and the scalar residual
# variance are held fixed at their null REML estimates. Approximate F tests
# are computed using the residual degrees of freedom of the full marker model.
.sgwas_marker_test <- function(fit, g, h = NULL, bivariate = FALSE) {
  g <- as.matrix(g); X <- fit$X; w <- fit$weights; m <- ncol(g)
  two <- bivariate || !is.null(h)
  
  if (bivariate) {
    B1 <- rbind(crossprod(X, g * w[[1]]), crossprod(X, g * w[[2]]))
    B2 <- rbind(crossprod(X, g * w[[2]]), crossprod(X, g * w[[3]]))
    S11 <- colSums(g * g * w[[1]]); S12 <- colSums(g * g * w[[2]]); S22 <- colSums(g * g * w[[3]])
    t1 <- as.numeric(crossprod(g, fit$Py[, 1])); t2 <- as.numeric(crossprod(g, fit$Py[, 2]))
  } else {
    B1 <- crossprod(X, g * w[[1]]); S11 <- colSums(g * g * w[[1]])
    t1 <- as.numeric(crossprod(g, fit$Py[, 1]))
    if (two) {
      h <- as.matrix(h); B2 <- crossprod(X, h * w[[1]])
      S12 <- colSums(g * h * w[[1]]); S22 <- colSums(h * h * w[[1]])
      t2 <- as.numeric(crossprod(h, fit$Py[, 1]))
    }
  }
  
  S11 <- S11 - colSums(B1 * (fit$Ainv %*% B1))
  df_null <- fit$nobs - fit$nfixed
  df <- df_null - if (two) 2L else 1L
  if (df <= 0L) stop("Insufficient residual degrees of freedom")
  
  # GWASpoly / EMMAX P3D: hold null REML variance fixed rather than re-estimating per marker
  null_variance <- fit$rss / df_null
  
  if (two) {
    S12 <- S12 - colSums(B1 * (fit$Ainv %*% B2))
    S22 <- S22 - colSums(B2 * (fit$Ainv %*% B2))
    det <- S11 * S22 - S12 * S12
    valid <- is.finite(det) & S11 > 1e-10 & S22 > 1e-10 & det > 1e-10 * S11 * S22
    det[!valid] <- NA_real_
    
    b1 <- (S22 * t1 - S12 * t2) / det
    b2 <- (S11 * t2 - S12 * t1) / det
    explained <- pmax(0, b1 * t1 + b2 * t2)
    
    se1 <- sqrt(null_variance * S22 / det)
    se2 <- sqrt(null_variance * S11 / det)
    pjoint <- pf((explained / 2) / null_variance, 2, df, lower.tail = FALSE)
    p1 <- pf((b1 / se1)^2, 1, df, lower.tail = FALSE)
    p2 <- pf((b2 / se2)^2, 1, df, lower.tail = FALSE)
  } else {
    valid <- is.finite(S11) & S11 > 1e-10; S11[!valid] <- NA_real_
    b1 <- t1 / S11
    se1 <- sqrt(null_variance / S11)
    p1 <- pf((b1 / se1)^2, 1, df, lower.tail = FALSE)
    b2 <- se2 <- p2 <- pjoint <- rep(NA_real_, m)
  }
  
  valid <- valid & is.finite(null_variance) & null_variance > .Machine$double.eps
  out <- data.frame(Beta1 = b1, SE1 = se1, P1 = p1, Beta2 = b2, SE2 = se2, P2 = p2, PJoint = pjoint,
                    DF_Residual = df, Status = ifelse(valid, "ok", "rank_or_residual_failure"))
  out[!valid, 1:7] <- NA_real_
  out
}

.sgwas_association_pair <- function(Y,X,G,map,eig,mode,chunk_size,maxit,
                                    model_effect="Add",ploidy=2L,dosage=NULL,max_geno_freq=0.99) {
  models <- .sgwas_marker_models(model_effect,ploidy)
  if(any(models!="additive") && (is.null(dosage) || !identical(dim(dosage),dim(G))))
    stop("Dominance scans require original dosages, with missing calls retained as NA")
  
  var_y1 <- var(Y[,1], na.rm = TRUE)
  var_y2 <- var(Y[,2], na.rm = TRUE)
  calc_pve <- function(beta, var_x, var_y) {
    if (!is.finite(var_y) || var_y <= 0) return(rep(NA_real_, length(beta)))
    pve <- (beta^2 * var_x) / var_y
    pve[!is.finite(pve)] <- NA_real_
    pmin(pve, 1)
  }
  
  fits <- list(Trait=.sgwas_fit_univariate(Y[,1],X,eig),
               Microbial=.sgwas_fit_univariate(Y[,2],X,eig))
  if(mode %in% c("sgwas-joint","sgwas-all")) fits$Bivariate <- .sgwas_fit_bivariate(Y,X,eig,maxit)
  M <- as.numeric(scale(Y[,2])); XX <- cbind(X,Microbial_score=M)
  if(mode %in% c("sgwas-interaction","sgwas-all")) {
    if(qr(XX)$rank<ncol(XX)) stop("Microbial score is collinear with a phenotype covariate or genotype PC")
    fits$Interaction <- .sgwas_fit_univariate(Y[,1],XX,eig)
  }
  chunks <- ceiling(ncol(G)/chunk_size)
  res <- vector("list",chunks*length(models)); result_index <- 0L
  for(model in models) for(k in seq_len(chunks)) {
    js <- seq.int((k-1L)*chunk_size+1L,min(k*chunk_size,ncol(G)))
    coding <- if(model=="additive")
      list(G=G[,js,drop=FALSE],n0=rep(NA_integer_,length(js)),n1=rep(NA_integer_,length(js)),
           status=rep("ok",length(js))) else
             .sgwas_encode_dominance(dosage[,js,drop=FALSE],model,ploidy,max_geno_freq)
    g <- coding$G
    var_x <- apply(g, 2, var, na.rm = TRUE)
    gr <- crossprod(eig$vectors,g)
    a <- .sgwas_marker_test(fits$Trait$fit,gr); b <- .sgwas_marker_test(fits$Microbial$fit,gr)
    r <- cbind(map[js,,drop=FALSE],Marker_Model=model,N_Group0=coding$n0,N_Group1=coding$n1,
               Marker_Coding_Status=coding$status,
               Beta_Trait=a$Beta1,SE_Trait=a$SE1,PVE_Trait=calc_pve(a$Beta1, var_x, var_y1),P_Trait=a$P1,Status_FER=a$Status,
               Beta_Microbial=b$Beta1,SE_Microbial=b$SE1,PVE_Microbial=calc_pve(b$Beta1, var_x, var_y2),P_Microbial=b$P1,Status_Microbial=b$Status)
    if(!is.null(fits$Bivariate)) {
      z <- .sgwas_marker_test(fits$Bivariate$fit,gr,bivariate=TRUE); s <- fits$Bivariate$scales
      b1_orig <- z$Beta1*s[1]; b2_orig <- z$Beta2*s[2]
      r <- cbind(r,Beta_Trait_joint=b1_orig,SE_Trait_joint=z$SE1*s[1],PVE_Trait_joint=calc_pve(b1_orig, var_x, var_y1),P_Trait_joint=z$P1,
                 Beta_Microbial_joint=b2_orig,SE_Microbial_joint=z$SE2*s[2],PVE_Microbial_joint=calc_pve(b2_orig, var_x, var_y2),P_Microbial_joint=z$P2,
                 P_Bivariate=z$PJoint,Status_Bivariate=z$Status)
    }
    if(!is.null(fits$Interaction)) {
      z <- .sgwas_marker_test(fits$Interaction$fit,gr,crossprod(eig$vectors,g*M))
      r <- cbind(r,Beta_SNP_at_mean=z$Beta1,SE_SNP_at_mean=z$SE1,PVE_SNP_at_mean=calc_pve(z$Beta1, var_x, var_y1),P_SNP_at_mean=z$P1,
                 Beta_SNPxM=z$Beta2,SE_SNPxM=z$SE2,P_SNPxM=z$P2,
                 P_SNP_joint=z$PJoint,Status_Interaction=z$Status)
    }
    bad <- coding$status!="ok"
    if(any(bad)) {
      r[bad,grepl("^(Beta_|SE_|PVE_|P_)",names(r))] <- NA_real_
      r[bad,startsWith(names(r),"Status_")] <- coding$status[bad]
    }
    result_index <- result_index+1L; res[[result_index]] <- r
  }
  list(results=do.call(rbind,res),fits=fits,microbial_center=mean(Y[,2]),microbial_scale=sd(Y[,2]))
}


.sgwas_association_plots <- function(results,file,fdr,bonferroni,suggestive,total_tests,format="pdf",fdr_level = 0.05) {
  tests <- names(results)[grepl("^P_",names(results)) & !grepl("_(BH|Bonferroni)$",names(results))]
  if(format=="pdf") pdf(file,width=11,height=5) else
    tiff(file,width=11,height=5,units="in",res=200,compression="lzw")
  on.exit(dev.off(),add=TRUE)
  if(!"Marker_Model" %in% names(results)) results$Marker_Model <- "additive"
  for(model in unique(results$Marker_Model)) for(test in tests) {
    rr <- results[results$Marker_Model==model,,drop=FALSE]
    p <- rr[[test]]; use <- is.finite(p) & p>=0 & p<=1
    if(!any(use)) next
    d <- rr[use,,drop=FALSE]; p <- p[use]
    chr <- unique(d$CHROM); ix <- order(match(d$CHROM,chr),d$POS)
    ypos <- -log10(pmax(p[ix],.Machine$double.xmin))
    colors <- c("#276D83","#84A98C")[1+(match(d$CHROM[ix],chr)%%2)]
    par(mfrow=c(1,2),mar=c(5,4,3,1))
    plot(seq_along(ix),ypos,pch=20,cex=.55,col=colors,xaxt="n",xlab="Chromosome",ylab="-log10(p)",
         main=paste(model,test,sep=" / "))
    ticks <- vapply(chr,function(ch) mean(which(d$CHROM[ix]==ch)),numeric(1))
    axis(1,at=ticks,labels=chr,las=2,cex.axis=.7)
    
    if(bonferroni && total_tests[[test]]>0) abline(h=-log10(.05/total_tests[[test]]),col="#C14953",lty=2)
    sig <- rr[[paste0(test,"_BH")]][use]<=fdr_level
    if(fdr && any(sig,na.rm=TRUE)) abline(h=-log10(max(p[sig],na.rm=TRUE)),col="#6E44AA",lty=3)
    if(is.finite(suggestive)) abline(h=suggestive,col="grey60",lty=3)
    observed <- sort(-log10(pmax(p,.Machine$double.xmin)))
    expected <- sort(-log10(ppoints(length(p))))
    plot(expected,observed,pch=20,cex=.6,xlab="Expected -log10(p)",ylab="Observed -log10(p)",main="QQ plot")
    abline(0,1,col="grey50")
  }
}

.sgwas_report_specs <- function(results) {
  specs <- data.frame(
    Test=c("Trait_marginal","Microbial_marginal","Trait_joint","Microbial_joint","Bivariate_joint",
           "SNP_at_mean","Interaction","SNP_joint"),
    P=c("P_Trait","P_Microbial","P_Trait_joint","P_Microbial_joint","P_Bivariate",
        "P_SNP_at_mean","P_SNPxM","P_SNP_joint"),
    Beta=c("Beta_Trait","Beta_Microbial","Beta_Trait_joint","Beta_Microbial_joint","",
           "Beta_SNP_at_mean","Beta_SNPxM",""),
    SE=c("SE_Trait","SE_Microbial","SE_Trait_joint","SE_Microbial_joint","",
         "SE_SNP_at_mean","SE_SNPxM",""),
    PVE=c("PVE_Trait","PVE_Microbial","PVE_Trait_joint","PVE_Microbial_joint","",
          "PVE_SNP_at_mean","",""),
    Status=c("Status_FER","Status_Microbial","Status_Bivariate","Status_Bivariate","Status_Bivariate",
             "Status_Interaction","Status_Interaction","Status_Interaction"),
    DF=c(1L,1L,1L,1L,2L,1L,1L,2L),stringsAsFactors=FALSE)
  if("P_Multitrait" %in% names(results)) {
    d <- unique(results$N_Traits)
    if(length(d)!=1L || !is.finite(d) || d<2L) stop("Invalid multi-trait report dimension")
    specs <- rbind(specs,data.frame(Test="Multitrait_joint",P="P_Multitrait",Beta="",SE="",PVE="",
                                    Status="Status_Multitrait",DF=as.integer(d),stringsAsFactors=FALSE))
  }
  specs[specs$P %in% names(results),,drop=FALSE]
}

# Reporting adapter only: fitted joint/interaction tests are never relabelled
# as native GWASpoly fits, nor are joint p-values turned into separate-trait scans.
.sgwas_gwas_results_reports <- function(results,outdir,all_results,fdr,bonferroni,suggestive,
                                        formats=c("pdf","tiff"),fdr_level=0.05) {
  specs <- .sgwas_report_specs(results)
  models <- unique(results$Marker_Model)
  root <- outdir; dir.create(root,showWarnings=FALSE)
  manifest <- list()
  for(i in seq_len(nrow(specs))) {
    spec <- specs[i,]; pcol <- spec$P
    folder <- file.path(root,spec$Test)
    for(sub in c("scores_effects","sigFDR","sigBonferroni","sigSuggestive","qqplots","manplots"))
      dir.create(file.path(folder,sub),recursive=TRUE,showWarnings=FALSE)
    label <- gsub("[^A-Za-z0-9_.-]","_",paste(results$Trait[1],results$Microbial_Trait[1],spec$Test,sep="_"))
    p <- results[[pcol]]; valid <- is.finite(p) & p>=0 & p<=1
    score <- rep(NA_real_,length(p)); score[valid] <- -log10(pmax(p[valid],.Machine$double.xmin))
    beta <- if(nzchar(spec$Beta)) results[[spec$Beta]] else rep(NA_real_,length(p))
    se <- if(nzchar(spec$SE)) results[[spec$SE]] else rep(NA_real_,length(p))
    pve <- if(nzchar(spec$PVE)) results[[spec$PVE]] else rep(NA_real_,length(p))
    bh <- results[[paste0(pcol,"_BH")]]; bon <- results[[paste0(pcol,"_Bonferroni")]]
    good <- is.finite(all_results[[pcol]])
    significant <- valid & is.finite(bh) & bh <= fdr_level
    fdr_cut <- if(any(significant)) -log10(max(p[significant],.Machine$double.xmin)) else Inf
    bon_cut <- if(sum(good)) -log10(.05/sum(good)) else Inf
    map_names <- intersect(c("SNP","CHROM","POS","REF","ALT","MAF"),names(results))
    wide <- results[!duplicated(results$SNP),map_names,drop=FALSE]
    long <- data.frame(Trait=results$Trait,Microbial_Trait=results$Microbial_Trait,Test=spec$Test,
                       Model=results$Marker_Model,Marker=results$SNP,Chrom=results$CHROM,Position=results$POS,
                       Ref=results$REF,Alt=results$ALT,MAF=results$MAF,Score=score,Effect=beta,SE=se,PVE=pve,
                       Pvalue=p,FDR=bh,Bonferroni=bon,DF_Test=spec$DF,Status=results[[spec$Status]],
                       Inference=if(paste0(pcol,"_Inference") %in% names(results)) results[[paste0(pcol,"_Inference")]] else results$Inference,
                       stringsAsFactors=FALSE)
    # Multi-trait imputation mode has different observed counts for marginal tests.
    for(field in c("N","N_Cohort","N_Observed_Trait","N_Observed_Phenotypes"))
      if(paste0(pcol,"_",field) %in% names(results)) long[[field]] <- results[[paste0(pcol,"_",field)]]
    for(model in models) {
      ix <- match(wide$SNP,results$SNP[results$Marker_Model==model])
      rows <- which(results$Marker_Model==model)[ix]
      for(nm in c("scores","effects","SE","PVE","pvalue","FDR","Bonferroni","status")) {
        v <- switch(nm,scores=score,effects=beta,SE=se,PVE=pve,pvalue=p,FDR=bh,Bonferroni=bon,status=long$Status)
        wide[[paste0(model,"_",nm)]] <- v[rows]
      }
    }
    .sgwas_write_tsv(wide,file.path(folder,"scores_effects",paste0("score_effects_",label,".txt")))
    .sgwas_write_tsv(long,file.path(folder,"scores_effects",paste0("score_effects_",label,"_long.tsv")))
    cuts <- data.frame(Test=spec$Test,P_Column=pcol,N_Tests=sum(good), FDR_Cutoff = fdr_cut, FDR_Level = fdr_level,
                       Bonferroni_0.05=bon_cut,Suggestive=suggestive,
                       Family="All eligible markers, marker models, trait pairs and ploidies in this call; one P column")
    .sgwas_write_tsv(cuts,file.path(folder,"thresholds.tsv"))
    for(kind in c("FDR","Bonferroni","Suggestive")) {
      keep <- switch(kind,
                     FDR = is.finite(bh) & bh <= fdr_level,
                     Bonferroni = is.finite(bon) & bon <= .05,
                     Suggestive = valid & score >= suggestive)
      threshold <- switch(kind,FDR=fdr_cut,Bonferroni=bon_cut,Suggestive=suggestive)
      selected <- long[keep,,drop=FALSE]; selected$Threshold <- rep(threshold,nrow(selected))
      suffix <- switch(kind,FDR=paste0("fdr",fdr_level),Bonferroni="Bonferroni0.05",Suggestive=paste0("sugg",suggestive))
      .sgwas_write_tsv(selected,file.path(folder,paste0("sig",kind),paste0("Significant_effect_",label,"_",suffix,".txt")))
    }
    
    qq <- list(); metric <- list()
    for(model in models) {
      pv <- sort(p[valid & results$Marker_Model==model]); nn <- length(pv)
      if(nn) {
        rank <- seq_len(nn)
        q <- data.frame(Model=model,Expected=-log10(ppoints(nn)),Observed=-log10(pmax(pv,.Machine$double.xmin)),
                        Lower=-log10(qbeta(.975,rank,nn+1-rank)),Upper=-log10(qbeta(.025,rank,nn+1-rank)))
        qq[[model]] <- q
      }
      metric[[model]] <- data.frame(model=model,Trait=results$Trait[1],Test=spec$Test,N_Tests=nn,
                                    dev_norm_avGW=if(nn) mean(q$Observed-q$Expected) else NA_real_)
    }
    .sgwas_write_tsv(do.call(rbind,metric),file.path(folder,"qqplots",paste0(label,"_qqplot_metric.txt")))
    manifest[[i]] <- data.frame(Test=spec$Test,P_Column=pcol,DF_Test=spec$DF,Directory=spec$Test,
                                Effect=if(nzchar(spec$Beta)) spec$Beta else "NA: omnibus test has multiple coefficients",
                                Engine="systemsGWAS REML/P3D GLS",Report_Style="GWASpoly-style tables and systemsGWAS facets")
    if(!length(formats)) next
    save_plot <- function(plot,sub,name,width,height) {
      for(format in formats) {
        opts <- list(filename=file.path(folder,sub,paste0(name,".",format)),plot=plot,
                     width=width,height=height,units="in",device=format,dpi=300,limitsize=FALSE)
        if(format=="tiff") opts$compression <- "lzw"
        do.call(ggplot2::ggsave,opts)
      }
    }
    if(!any(valid)) {
      blank <- ggplot2::ggplot()+ggplot2::annotate("text",x=0,y=0,label="No testable markers")+
        ggplot2::theme_void()+ggplot2::ggtitle(label)
      save_plot(blank,"manplots",paste0("manplot_",label),12,4)
      save_plot(blank,"manplots",paste0("QTLprofile_",label),12,4)
      save_plot(blank,"qqplots",paste0("QQplot_",label),5,4)
      next
    }
    d <- long[valid, , drop = FALSE]
    d$Model <- factor(d$Model, levels = models)
    chromosomes <- unique(as.character(results$CHROM))
    chr_num <- suppressWarnings(as.numeric(sub("^(chr|chrom)_?", "", chromosomes, ignore.case = TRUE)))
    chromosomes <- chromosomes[order(is.na(chr_num), chr_num, chromosomes)]
    d$Chrom <- factor(as.character(d$Chrom), levels = chromosomes, ordered = TRUE)
    d <- d[order(d$Model, d$Chrom, d$Position), ]
    palette <- setNames(rep(c("cornflowerblue", "orange"), length.out = length(chromosomes)), chromosomes)
    h <- data.frame(Threshold = c("FDR", "Bonferroni", "Suggestive"), Value = c(fdr_cut, bon_cut, suggestive))
    h <- h[c(fdr, bonferroni, TRUE) & is.finite(h$Value), , drop = FALSE]
    
    # Threshold colour palette
    chrom_cols <- rep(c("cornflowerblue", "orange"), length.out = length(unique(d$Chrom)))
    names(chrom_cols) <- unique(d$Chrom)
    palette <- c(chrom_cols, FDR = "grey10", Bonferroni = "tomato", Suggestive = "cornflowerblue")
    # Base Manhattan plot
    base <- ggplot2::ggplot(d, ggplot2::aes(x = Position, y = Score, colour = Chrom)) +
      ggplot2::facet_grid(Model ~ Chrom, scales = "free_x", space = "free_x") +
      ggplot2::scale_colour_manual(values = palette, breaks = h$Threshold) +
      ggplot2::scale_linetype_manual(values = c(FDR = "dashed", Bonferroni = "dashed", Suggestive = "dashed"), breaks = h$Threshold) +
      ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::theme_bw(base_size = 18) +
      ggplot2::theme(panel.spacing.x = grid::unit(0, "lines"), panel.spacing.y = grid::unit(0.5, "lines"),
                     panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "white"),
                     axis.title.x = ggplot2::element_text(size = 15, margin = ggplot2::margin(t = 12)),
                     axis.title.y = ggplot2::element_text(size = 15, margin = ggplot2::margin(r = 12)),
                     axis.text.y = ggplot2::element_text(size = 15, colour = "black"), strip.background = ggplot2::element_blank(),
                     strip.text.x = ggplot2::element_text(size = 15, colour = "black"), strip.text.y = ggplot2::element_text(size = 15, colour = "black"),
                     legend.position = "top", legend.justification = "right",
                     legend.box.just = "right", legend.background = ggplot2::element_rect(fill = "white", colour = "white"), legend.title = ggplot2::element_text(size = 17),
                     legend.text = ggplot2::element_text(size = 15), legend.key.width = grid::unit(1.5, "cm"),
                     legend.key.height = grid::unit(0.7, "cm"), plot.title = ggplot2::element_text(size = 15, hjust = 0.5, face = "plain", margin = ggplot2::margin(b = 15)),
                     plot.margin = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)) +
      ggplot2::labs(title = label, x = "Chromosomes and Physical Map Position (bp)", y = "-log10(P)", colour = "Threshold", linetype = "Threshold") +
      ggplot2::coord_cartesian(clip = "off")
    
    # Manhattan plot
    man <- base + ggplot2::geom_point(size = 2, shape = 20, show.legend = FALSE) +
      ggplot2::geom_hline(data = h, ggplot2::aes(yintercept = Value, colour = Threshold, linetype = Threshold), linewidth = 0.35, inherit.aes = FALSE) +
      ggplot2::scale_colour_manual(values = palette, breaks = h$Threshold) +
      ggplot2::scale_linetype_manual(values = c(FDR = "dashed", Bonferroni = "dashed", Suggestive = "dashed"),breaks = h$Threshold)
    save_plot(man, "manplots", paste0("manplot_", label), 12, max(4, 1.4 * length(models) + 1.5))
    
    # QTL profile plot
    d$Smoothed <- d$Score
    for (ix in split(seq_len(nrow(d)), interaction(d$Model, d$Chrom, drop = TRUE))) {
      if (length(ix) > 1L) {
        d$Smoothed[ix] <- c(NA_real_, (head(d$Score[ix], -1) + tail(d$Score[ix], -1)) / 2 )
      }
    }
    profile <- base +
      ggplot2::geom_line(data = d, ggplot2::aes(y = Smoothed, group = interaction(Model, Chrom)),
                         linewidth = 1.0, show.legend = FALSE, na.rm = TRUE)
    save_plot(profile, "manplots", paste0("QTLprofile_", label), 12, max(4, 1.4 * length(models) + 1.5))
    
    # QQ plots
    q <- do.call(rbind, qq)
    q$Model <- factor( q$Model, levels = models)
    qp <- ggplot2::ggplot(q, ggplot2::aes(x = Expected, y = Observed)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper), alpha = 0.12, fill = "grey60") +
      ggplot2::geom_point(size = 2.5, colour = "cornflowerblue") +
      ggplot2::geom_abline(intercept = 0, slope = 1, colour = "red", linewidth = 0.8, linetype = "dashed") +
      ggplot2::facet_wrap(~Model, ncol = 1) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(strip.background = ggplot2::element_blank(), legend.box.spacing = grid::unit(0.2, "cm"),  panel.grid = ggplot2::element_blank(),
                     legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                     strip.text = ggplot2::element_text(size = 18, colour = "black"),
                     plot.title = ggplot2::element_text(size = 5, hjust = 0.5),
                     axis.title = ggplot2::element_text(size = 18), axis.text = ggplot2::element_text(size = 14)) +
      ggplot2::labs(title = label, x = "Expected -log10(P)", y = "Observed -log10(P)")
    save_plot(qp, "qqplots", paste0("QQplot_", label), 5, max(4, 2 * length(models) + 1))
  }
  .sgwas_write_tsv(do.call(rbind,manifest),file.path(root,"REPORT_MANIFEST.tsv"))
}

.sgwas_integrated_gwas <- function(args,prepared,mode,n_pcs,max_missing,chunk_size,maxit,max_geno_freq=0.99,
                                   output_style="verbose",plot_formats=c("pdf","tiff")) {
  effects <- .sgwas_validate_integrated_models(args$method,args$model_effect)
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
  marginal_only <- identical(mode, "sgwas-marginal")
  if(!marginal_only && prepared$summary$Trait %in% c(traits,covs))
    stop("Use a microbial_trait_name distinct from outcomes and covariates")
  if(args$cores>1) message("[systemsGWAS] Base-R REML/GLS scans use serial vectorised chunks; cores controls the legacy GWASpoly engine only.")
  if(marginal_only) {
    inference <- "standard_univariate"
    outdir <- file.path(args$wdir,paste0(args$projname,"_base_gwas"))
  } else {
    provenance <- tolower(prepared$summary$Source_Method)
    supervised <- isTRUE(provenance == "spls")
    inference <- if(supervised) "exploratory_supervised_score" else
      if(provenance %in% c("single","sum","pc1","pc2","pc3")) "unsupervised_score" else "unknown_weight_provenance"
    if(supervised) warning("Supervised microbial weights: nominal association p-values do not account for learning the score. Treat as exploratory unless trained independently.",call.=FALSE)
    if(inference=="unknown_weight_provenance") warning("Saved weights have unknown training provenance; verify independence from the measured outcome before confirmatory interpretation.",call.=FALSE)
    outdir <- file.path(args$wdir,paste0(args$projname,"_integrated_gwas"))
  }
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  complete <- FALSE
  run_status <- function(status) .sgwas_write_tsv(data.frame(Status=status,Mode=mode,Time=as.character(Sys.time())),
                                                  file.path(outdir,"RUN_STATUS.tsv"))
  run_status("running")
  on.exit(if(!complete) run_status("failed"),add=TRUE)
  results <- list(); summaries <- list(); audit <- list(); allfits <- list(); outputs <- list(); idx <- 0L
  score <- if(marginal_only) NULL else prepared$phenotype[,1:2,drop=FALSE]
  micro_name <- if(marginal_only) "None" else names(score)[2]
  for(ploidy in unique(as.numeric(args$ploidy_levels))) {
    raw <- .sgwas_read_dosage(args[[paste0("genofile_",ploidy,"x")]],ploidy)
    ploidy_folder <- file.path(outdir, paste0(ploidy, "x"))
    dir.create(ploidy_folder, showWarnings = FALSE)
    for(trait in traits) {
      idx <- idx+1L; label <- paste0("pair",idx,"_",ploidy,"x")
      folder <- ploidy_folder
      
      ids <- if(marginal_only) unique(ph$Plant_ID) else score$Plant_ID
      y <- .sgwas_id_values(ph,trait,ids)
      C <- matrix(numeric(),nrow=length(ids),ncol=0L,dimnames=list(NULL,character()))
      if(length(covs)) C <- vapply(covs,function(nm) .sgwas_id_values(ph,nm,ids),numeric(length(ids)))
      colnames(C) <- covs
      has_micro <- if(marginal_only) rep(TRUE,length(ids)) else is.finite(score[[2]])
      keep <- is.finite(y) & has_micro & rowSums(!is.finite(C))==0L & ids %in% rownames(raw$G)
      audit[[idx]] <- data.frame(Pair=label,Trait=trait,Microbial_Trait=micro_name,Ploidy=ploidy,
                                 Plant_ID=ids,Observed_Trait=is.finite(y),Complete_Covariates=rowSums(!is.finite(C))==0L,
                                 Has_Genotype=ids %in% rownames(raw$G),Included=keep)
      ids <- ids[keep]; C <- C[keep,,drop=FALSE]
      Y <- if(marginal_only) cbind(Trait=y[keep]) else cbind(Trait=y[keep],Microbial=score[[2]][keep])
      if(length(ids)<max(10L,ncol(C)+n_pcs+5L)) stop("Too few complete samples for ",trait)
      if(any(!is.finite(apply(Y,2,sd))) || any(apply(Y,2,sd)<=0)) stop("Outcome(s) must vary on the analysed cohort")
      if(ncol(C)) {
        cs <- apply(C,2,sd)
        if(any(!is.finite(cs)) || any(cs<=0)) stop("Constant or invalid phenotype covariate on the analysed cohort")
        C <- scale(C)
      }
      qc <- .sgwas_dosage_qc(raw,ids,ploidy,args$maf,max_missing,max_geno_freq=max_geno_freq)
      tr_prefix <- gsub("[^A-Za-z0-9_.-]", "_", trait)
      .sgwas_write_tsv(qc$ledger, file.path(folder, paste0(tr_prefix, "_marker_qc.tsv")))
      if("Dom" %in% effects) {
        if(any(qc$map$N_Fractional>0)) warning(
          "Dominance encodings round observed fractional dosages to integer calls, as GWASpoly does; genotype probabilities are not propagated.",call.=FALSE)
        if(all(qc$map$N_Het_calls==0) && !any(qc$map$N_Fractional>0)) message(
          "[systemsGWAS] No heterozygous calls in this cohort: dominance encodings duplicate the additive contrasts; they cannot establish dominance.")
      }
      af <- colMeans(qc$G, na.rm = TRUE) / ploidy
      marker_denom <- ploidy * af * (1 - af)
      total_denom <- sum(marker_denom)
      cp <- tcrossprod(qc$centered)
      full_eig <- .sgwas_kinship_eigen(cp, total_denom)
      X <- cbind(Intercept=1,C)
      if(n_pcs) {
        if(n_pcs>sum(full_eig$values>max(full_eig$values)*1e-8)) stop("Too few nonzero genomic PCs")
        pcs <- full_eig$vectors[,seq_len(n_pcs),drop=FALSE]
        colnames(pcs) <- paste0("Genotype_PC",seq_len(n_pcs)); X <- cbind(X,pcs)
      }
      if(qr(X)$rank<ncol(X)) stop("Collinear covariates/PCs in integrated GWAS for ",trait)
      sample_df <- if(marginal_only) {
        data.frame(Plant_ID=ids,Y,X,check.names=FALSE)
      } else {
        data.frame(Plant_ID=ids,Y,Microbial_standardised=as.numeric(scale(Y[,2])),X,check.names=FALSE)
      }
      .sgwas_write_tsv(sample_df, file.path(folder, paste0(tr_prefix, "_analysed_samples.tsv")))
      
      groups <- if(args$LOCO) unique(qc$map$CHROM) else "All"
      if(args$LOCO && length(groups)<2L) stop("LOCO requires retained markers on at least two chromosomes")
      result_chr <- list(); fit_chr <- list()
      for(ch in groups) {
        selected <- if(args$LOCO) which(qc$map$CHROM==ch) else seq_len(ncol(qc$G))
        eg <- if(args$LOCO) {
          .sgwas_kinship_eigen(
            cp - tcrossprod(qc$centered[, selected, drop = FALSE]),
            total_denom - sum(marker_denom[selected]))
        } else {full_eig}
        message("[systemsGWAS] ",trait," + ",micro_name,"; ",ploidy,"x; ",length(ids)," samples; ",
                length(selected)," markers; ",if(args$LOCO) paste0("LOCO ",ch) else "whole-genome kinship")
        
        
        fit <- tryCatch({
          if(marginal_only) {
            .sgwas_marginal_scan(Y[,1],X,qc$centered[,selected,drop=FALSE],qc$map[selected,,drop=FALSE],
                                 eg,chunk_size,model_effect=effects,ploidy=ploidy,
                                 dosage=if("Dom" %in% effects) qc$G[,selected,drop=FALSE] else NULL,
                                 max_geno_freq=max_geno_freq)
          } else {
            .sgwas_association_pair(Y,X,qc$centered[,selected,drop=FALSE],qc$map[selected,,drop=FALSE],
                                    eg,mode,chunk_size,maxit,model_effect=effects,ploidy=ploidy,
                                    dosage=if("Dom" %in% effects) qc$G[,selected,drop=FALSE] else NULL,
                                    max_geno_freq=max_geno_freq)
          }
        },error=function(e) {
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
  grp <- interaction(results$Pair, results$Marker_Model, drop = TRUE)
  for(p in pcols) {
    good <- is.finite(results[[p]]); total_tests[[p]] <- sum(good)
    fdr_val <- rep(NA_real_, nrow(results))
    for (ix in split(seq_len(nrow(results)), grp)) {
      fdr_val[ix] <- .sgwas_fdr(results[[p]][ix])
    }
    results[[paste0(p, "_BH")]] <- fdr_val
    bon_val <- rep(NA_real_, nrow(results)); bon_val[good] <- p.adjust(results[[p]][good], method = "bonferroni")
    results[[paste0(p, "_Bonferroni")]] <- bon_val
  }
  .sgwas_write_tsv(results,file.path(outdir,"ASSOCIATIONS.tsv"))
  summary <- do.call(rbind,summaries); .sgwas_write_tsv(summary,file.path(outdir,"NULL_MODELS.tsv"))
  cohort <- do.call(rbind,audit); .sgwas_write_tsv(cohort,file.path(outdir,"SAMPLE_AUDIT.tsv"))
  for(i in seq_along(outputs)) {
    r <- results[results$Pair==unique(results$Pair)[i],,drop=FALSE]
    tr <- gsub("[^A-Za-z0-9_.-]", "_", r$Trait[1])
    .sgwas_write_tsv(r,file.path(outputs[[i]],paste0("associations_",tr,".tsv")))
    # Around line 947 inside .sgwas_integrated_gwas():
    if(output_style=="verbose") .sgwas_gwas_results_reports(r,outputs[[i]],results,args$fdr,
                                                            args$bonferroni,args$suggestive,plot_formats,fdr_level=args$fdr_level)
    
    # Preserve the compact PDF as a quick overview, alongside the faceted reports.
    if("pdf" %in% plot_formats)
      .sgwas_association_plots(r,file.path(outputs[[i]],paste0("Manhattan_QQ_",tr,".pdf")),args$fdr,args$bonferroni,args$suggestive,total_tests,fdr_level=args$fdr_level)
    if(output_style=="compact" && "tiff" %in% plot_formats)
      .sgwas_association_plots(r,file.path(outputs[[i]],paste0("Manhattan_QQ_",tr,"_%03d.tiff")),args$fdr,args$bonferroni,args$suggestive,total_tests,"tiff",fdr_level=args$fdr_level)
  }
  metadata <- list(mode=mode,arguments=args,n_pcs=n_pcs,max_missing=max_missing,chunk_size=chunk_size,maxit=maxit,
                   max_geno_freq=max_geno_freq,marker_effects=effects,output_style=output_style,plot_formats=plot_formats,
                   inference=inference,microbiome=if(!is.null(prepared)) prepared$model else NULL,models=allfits,tests=total_tests,session=sessionInfo(),
                   method="VanRaden additive kinship REML null/P3D covariance shape; separate additive/dominance marker encodings; null REML residual scale and approximate F tests",
                   dominance_coding="round observed ALT dosage; d-dom-alt: dosage >= d; d-dom-ref: dosage > ploidy-d; missing calls imputed with encoded mean; centre each encoding",
                   multiple_testing="Per-trait and per-marker-model Storey qvalue FDR; global Bonferroni across all eligible tests")
  saveRDS(metadata,file.path(outdir,"MODEL_DETAILS.rds"))
  complete <- TRUE; run_status("complete")
  message("[systemsGWAS] Integrated association results: ",file.path(outdir,"ASSOCIATIONS.tsv"))
  invisible(list(results=results,null_models=summary,sample_audit=cohort,metadata=metadata,outdir=outdir))
}

# Multivariate phenotype GWAS. Trait-major observation order throughout:
# Var(vec(Y)) = G_trait (x) K + R_trait (x) I. No microbial score is required.
.sgwas_mt_cov <- function(par, d) {
  L <- matrix(0, d, d); ix <- lower.tri(L, diag=TRUE)
  L[ix] <- par; diag(L) <- exp(diag(L))
  list(cov=tcrossprod(L), L=L)
}

.sgwas_mt_pack <- function(V) {
  L <- t(chol(V)); diag(L) <- log(diag(L))
  L[lower.tri(L, diag=TRUE)]
}

.sgwas_mt_weights <- function(lambda, G, R) {
  d <- nrow(G); cr <- chol(R); ri <- backsolve(cr, diag(d))
  gg <- crossprod(ri, G %*% ri)
  eg <- eigen((gg+t(gg))/2, symmetric=TRUE)
  v <- 1 + outer(lambda, pmax(eg$values, 0))
  if(any(!is.finite(v)) || any(v<=0)) stop("Invalid multivariate covariance")
  B <- ri %*% eg$vectors; iv <- 1/v
  W <- array(0, c(length(lambda),d,d))
  for(a in seq_len(d)) for(b in seq_len(a)) {
    W[,a,b] <- as.numeric(iv %*% (B[a,]*B[b,]))
    W[,b,a] <- W[,a,b]
  }
  list(weights=W, logdet=length(lambda)*2*sum(log(diag(cr)))+sum(log(v)))
}

.sgwas_mt_gls <- function(Y, X, W) {
  n <- nrow(Y); d <- ncol(Y); p <- ncol(X)
  A <- matrix(0,p*d,p*d); WY <- matrix(0,n,d)
  blocks <- lapply(seq_len(d),function(a) seq.int((a-1L)*p+1L,a*p))
  for(a in seq_len(d)) for(b in seq_len(d)) {
    A[blocks[[a]],blocks[[b]]] <- crossprod(X,X*W[,a,b])
    WY[,a] <- WY[,a]+W[,a,b]*Y[,b]
  }
  ai <- .sgwas_inverse(A); r <- as.vector(crossprod(X,WY))
  beta <- ai$inverse %*% r
  residual <- Y-X%*%matrix(beta,p,d)
  Py <- matrix(0,n,d)
  for(a in seq_len(d)) for(b in seq_len(d)) Py[,a] <- Py[,a]+W[,a,b]*residual[,b]
  rss <- sum(residual*Py)
  if(!is.finite(rss) || rss<=0) stop("Nonpositive multivariate residual sum of squares")
  list(X=X,Y=Y,weights=W,Ainv=ai$inverse,beta=beta,Py=Py,rss=rss,
       logdetA=ai$logdet,nobs=n*d,nfixed=p*d)
}

# Missing outcomes are integrated out exactly, rather than treated as observed.
# A nuisance indicator for each missing cell gives the same REML objective and
# marker information as V[observed, observed], while retaining spectral blocks.
.sgwas_mt_missing_design <- function(Xr, U, missing) {
  n <- nrow(missing); d <- ncol(missing); ix <- which(missing,arr.ind=TRUE)
  Z <- matrix(0,n*d,nrow(ix))
  for(a in seq_len(d)) {
    js <- which(ix[,2]==a)
    if(length(js)) Z[seq.int((a-1L)*n+1L,a*n),js] <- t(U[ix[js,1],,drop=FALSE])
  }
  cbind(kronecker(diag(d),Xr),Z)
}

.sgwas_mt_gls_missing <- function(Y,X,W,design) {
  n <- nrow(Y); d <- ncol(Y)
  blocks <- lapply(seq_len(d),function(a) seq.int((a-1L)*n+1L,a*n))
  WD <- matrix(0,n*d,ncol(design)); WY <- matrix(0,n,d)
  for(a in seq_len(d)) for(b in seq_len(d)) {
    WD[blocks[[a]],] <- WD[blocks[[a]],]+design[blocks[[b]],,drop=FALSE]*W[,a,b]
    WY[,a] <- WY[,a]+Y[,b]*W[,a,b]
  }
  ai <- .sgwas_inverse(crossprod(design,WD))
  beta <- ai$inverse%*%crossprod(design,as.vector(WY))
  residual <- Y-matrix(design%*%beta,n,d)
  Py <- matrix(0,n,d)
  for(a in seq_len(d)) for(b in seq_len(d)) Py[,a] <- Py[,a]+W[,a,b]*residual[,b]
  rss <- sum(residual*Py)
  if(!is.finite(rss) || rss<=0) stop("Nonpositive observed-data multivariate RSS")
  list(X=X,Y=Y,weights=W,Ainv=ai$inverse,beta=beta,Py=Py,rss=rss,
       logdetA=ai$logdet,nobs=n*d,nfixed=ncol(design),WD=WD)
}

# Exact derivative of the REML objective, including the fixed-effect projection.
.sgwas_mt_gradient <- function(fit, lambda, LG, LR) {
  d <- ncol(fit$Y); W <- fit$weights; X <- fit$X
  BX <- if(is.null(fit$WD)) lapply(seq_len(d),function(a)
    do.call(cbind,lapply(seq_len(d),function(b) X*W[,a,b]))) else
      lapply(seq_len(d),function(a) fit$WD[seq.int((a-1L)*nrow(X)+1L,a*nrow(X)),,drop=FALSE])
  BA <- lapply(BX,function(x) x%*%fit$Ainv)
  dg <- dr <- matrix(0,d,d)
  for(a in seq_len(d)) for(b in seq_len(a)) {
    z <- W[,a,b]-rowSums(BX[[a]]*BA[[b]])-fit$Py[,a]*fit$Py[,b]
    dg[a,b] <- dg[b,a] <- .5*sum(lambda*z)
    dr[a,b] <- dr[b,a] <- .5*sum(z)
  }
  pack_gradient <- function(score,L) {
    grad <- 2*score%*%L; diag(grad) <- diag(grad)*diag(L)
    grad[lower.tri(grad,diag=TRUE)]
  }
  c(pack_gradient(dg,LG),pack_gradient(dr,LR))
}

.sgwas_fit_multitrait <- function(Y,X,eig,maxit=300L) {
  Y <- as.matrix(Y); X <- as.matrix(X); n <- nrow(Y); d <- ncol(Y)
  if(d<2L) stop("Multi-trait GWAS requires at least two traits")
  if(nrow(X)!=n || any(!is.finite(X))) stop("Invalid multi-trait trait/design matrix")
  missing <- !is.finite(Y); Y[missing] <- NA_real_
  if(any(rowSums(!missing)==0L)) stop("Samples with all selected traits missing carry no outcome information; exclude them")
  for(a in seq_len(d)) {
    observed <- !missing[,a]
    if(sum(observed)<=ncol(X)+1L || qr(X[observed,,drop=FALSE])$rank<ncol(X))
      stop("Too few observed values or rank-deficient covariates for trait ",colnames(Y)[a])
  }
  if(qr(X)$rank<ncol(X)) stop("Collinear multi-trait covariates")
  if(n-ncol(X)<=d) stop("Too few residual observations for the selected traits")
  scales <- apply(Y,2,sd,na.rm=TRUE); centers <- colMeans(Y,na.rm=TRUE)
  if(any(!is.finite(scales)) || any(scales<=0)) stop("All selected traits must vary on the analysed cohort")
  Ys <- sweep(sweep(Y,2,centers,"-"),2,scales,"/")
  complete <- rowSums(missing)==0L
  if(sum(complete)>ncol(X)+d && qr(X[complete,,drop=FALSE])$rank==ncol(X)) {
    residual <- qr.resid(qr(X[complete,,drop=FALSE]),Ys[complete,,drop=FALSE])
    ev <- eigen(crossprod(residual),symmetric=TRUE,only.values=TRUE)$values
    if(min(ev)<=max(ev)*1e-8) stop("Selected traits are linearly dependent or nearly collinear after covariate adjustment; remove redundant traits/composites")
  }
  if(!any(missing)) {
    residual <- qr.resid(qr(X),Ys); vr <- crossprod(residual)/(n-ncol(X))
  } else {
    residual <- matrix(NA_real_,n,d)
    for(a in seq_len(d)) {
      observed <- !missing[,a]
      residual[observed,a] <- qr.resid(qr(X[observed,,drop=FALSE]),Ys[observed,a])
    }
    # Pairwise covariance is only a starting value, never a completed dataset.
    vr <- cov(residual,use="pairwise.complete.obs")
    if(any(!is.finite(vr))) stop("Insufficient overlapping observations to initialize trait covariance")
    ev <- eigen(vr,symmetric=TRUE)
    vr <- tcrossprod(sweep(ev$vectors,2,sqrt(pmax(ev$values,1e-4)),"*"))
  }
  if(diff(range(eig$values))<1e-8) stop("Multi-trait genetic/residual covariance is not identifiable: kinship is proportional to identity")
  Ys[missing] <- 0 # Only a placeholder; missing-cell indicators remove its information.
  yr <- crossprod(eig$vectors,Ys); xr <- crossprod(eig$vectors,X)
  design <- if(any(missing)) .sgwas_mt_missing_design(xr,eig$vectors,missing) else NULL
  # Identifiability must hold in the residual subspace, not only in the full K.
  # For example, K = a*I + b*11' is unidentifiable with a fitted intercept.
  Q <- qr.Q(qr(xr)); QDQ <- crossprod(Q,Q*eig$values); rdf <- n-ncol(X)
  mean_k <- (sum(eig$values)-sum(diag(QDQ)))/rdf
  mean_k2 <- (sum(eig$values^2)-2*sum((Q*eig$values)^2)+sum(QDQ^2))/rdf
  if(mean_k2-mean_k^2<1e-10*max(1,mean_k2))
    stop("Multi-trait genetic/residual covariance is not identifiable after fixed-effect projection")
  k <- d*(d+1L)/2L; cached_par <- NULL; cached <- NULL
  evaluate <- function(par) {
    if(identical(par,cached_par)) return(cached)
    cached_par <<- par
    cached <<- tryCatch({
      G <- .sgwas_mt_cov(par[seq_len(k)],d); R <- .sgwas_mt_cov(par[k+seq_len(k)],d)
      w <- .sgwas_mt_weights(eig$values,G$cov,R$cov)
      fit <- if(is.null(design)) .sgwas_mt_gls(yr,xr,w$weights) else
        .sgwas_mt_gls_missing(yr,xr,w$weights,design)
      value <- .5*((fit$nobs-fit$nfixed)*log(2*pi)+w$logdet+fit$logdetA+fit$rss)
      gradient <- .sgwas_mt_gradient(fit,eig$values,G$L,R$L)
      if(!is.finite(value) || any(!is.finite(gradient))) stop("Nonfinite REML evaluation")
      list(objective=value,gradient=gradient,fit=fit,G=G$cov,R=R$cov)
    },error=function(e) list(objective=1e100,gradient=rep(0,length(par))))
    cached
  }
  ix <- lower.tri(matrix(0,d,d),diag=TRUE)
  diagonal <- row(matrix(0,d,d))[ix]==col(matrix(0,d,d))[ix]
  lower <- rep(ifelse(diagonal,-12,-10),2L); upper <- rep(ifelse(diagonal,4,10),2L)
  fractions <- c(.2,.5,.8)
  opts <- lapply(fractions,function(f) {
    s <- c(.sgwas_mt_pack(vr*f),.sgwas_mt_pack(vr*(1-f)))
    tryCatch(nlminb(s,function(par) evaluate(par)$objective,
                    gradient=function(par) evaluate(par)$gradient,lower=lower,upper=upper,
                    control=list(iter.max=maxit,eval.max=maxit*4L,rel.tol=1e-9)),error=identity)
  })
  attempts <- do.call(rbind,lapply(seq_along(opts),function(i) {
    z <- opts[[i]]
    if(inherits(z,"error")) data.frame(Start=i,Initial_G_Fraction=fractions[i],Convergence=NA_integer_,
                                       Objective=NA_real_,Iterations=NA_integer_,Message=conditionMessage(z)) else
                                         data.frame(Start=i,Initial_G_Fraction=fractions[i],Convergence=z$convergence,Objective=z$objective,
                                                    Iterations=z$iterations,Message=if(is.null(z$message)) "" else z$message)
  }))
  good <- vapply(opts,function(z) !inherits(z,"error") && z$convergence==0L &&
                   is.finite(z$objective) && z$objective<1e90,logical(1))
  if(!any(good)) stop(structure(list(message=paste0("Multi-trait REML did not converge; increase association_maxit or reduce redundant traits. ",
                                                    paste(attempts$Message,collapse="; ")),call=NULL,attempts=attempts),
                                class=c("sgwas_multitrait_convergence_error","error","condition")))
  best <- opts[[which(good)[which.min(vapply(opts[good],`[[`,numeric(1),"objective"))]]]
  ans <- evaluate(best$par); ans$gradient <- NULL
  ans$centers <- centers; ans$scales <- scales; ans$converged <- TRUE
  ans$optimizer <- "nlminb_REML_analytic_gradient_3_starts"; ans$attempts <- attempts
  # Flag near-singular covariance estimates even when optimization converges.
  scale_total <- max(diag(ans$G+ans$R))
  ans$boundary <- min(eigen(ans$G,symmetric=TRUE,only.values=TRUE)$values)<1e-5*scale_total ||
    min(eigen(ans$R,symmetric=TRUE,only.values=TRUE)$values)<1e-5*scale_total ||
    any(best$par-lower<1e-4 | upper-best$par<1e-4)
  ans$G_original <- ans$G*outer(scales,scales); ans$R_original <- ans$R*outer(scales,scales)
  dimnames(ans$G_original) <- dimnames(ans$R_original) <- list(colnames(Y),colnames(Y))
  ans$n_samples <- n; ans$n_observed <- sum(!missing); ans$n_missing <- sum(missing)
  ans$observed_per_trait <- colSums(!missing)
  ans$inference <- if(any(missing)) "observed_data_REML_with_conditional_imputation" else "complete_cases"
  ans$completed <- Y
  ans$imputation <- data.frame(Sample_Index=integer(),Trait_Index=integer(),Trait=character(),
                               Imputed_Value=double(),Imputation_SE=double())
  if(any(missing)) {
    cells <- which(missing,arr.ind=TRUE); js <- ncol(X)*d+seq_len(nrow(cells))
    # Negative nuisance coefficients are conditional means in standardized units.
    values <- -as.numeric(ans$fit$beta[js])*scales[cells[,2]]+centers[cells[,2]]
    variance <- ans$fit$rss/(ans$fit$nobs-ans$fit$nfixed)
    se <- sqrt(variance*diag(ans$fit$Ainv)[js])*scales[cells[,2]]
    ans$completed[missing] <- values
    ans$imputation <- data.frame(Sample_Index=cells[,1],Trait_Index=cells[,2],Trait=colnames(Y)[cells[,2]],
                                 Imputed_Value=values,Imputation_SE=se)
  }
  ans
}

# Block cross-products followed by a d-by-d solve per SNP. No n*d dense V.
.sgwas_mt_marker_test <- function(fit,g) {
  g <- as.matrix(g); d <- ncol(fit$Y); m <- ncol(g); W <- fit$weights; X <- fit$X
  B <- if(is.null(fit$WD)) lapply(seq_len(d),function(a)
    do.call(rbind,lapply(seq_len(d),function(b) crossprod(X,g*W[,b,a])))) else
      lapply(seq_len(d),function(a) crossprod(fit$WD[seq.int((a-1L)*nrow(X)+1L,a*nrow(X)),,drop=FALSE],g))
  AB <- lapply(B,function(z) fit$Ainv%*%z)
  S <- array(0,c(m,d,d)); score <- crossprod(g,fit$Py)
  for(a in seq_len(d)) for(b in seq_len(a)) {
    S[,a,b] <- colSums(g*g*W[,a,b])-colSums(B[[a]]*AB[[b]])
    S[,b,a] <- S[,a,b]
  }
  df_null <- fit$nobs-fit$nfixed; df <- df_null-d; variance <- fit$rss/df_null
  if(df<=0 || !is.finite(variance) || variance<=0) stop("Invalid multi-trait residual degrees of freedom/scale")
  beta <- se <- p <- matrix(NA_real_,m,d); joint <- rep(NA_real_,m)
  status <- rep("rank_or_residual_failure",m)
  for(j in seq_len(m)) {
    s <- matrix(S[j,,],d,d); dd <- diag(s)
    if(any(!is.finite(s)) || any(dd<=1e-10)) next
    # Scale-free rank check, followed by the original-scale covariance inverse.
    if(rcond(s/sqrt(outer(dd,dd)))<1e-10) next
    si <- tryCatch(.sgwas_inverse(s)$inverse,error=function(e) NULL)
    if(is.null(si)) next
    beta[j,] <- as.numeric(si%*%score[j,]); se[j,] <- sqrt(variance*diag(si))
    p[j,] <- pf((beta[j,]/se[j,])^2,1,df,lower.tail=FALSE)
    joint[j] <- pf(max(0,sum(beta[j,]*score[j,]))/(d*variance),d,df,lower.tail=FALSE)
    status[j] <- "ok"
  }
  list(beta=beta,se=se,p=p,pjoint=joint,df=df,status=status)
}

.sgwas_multitrait_scan <- function(Y,X,G,map,eig,chunk_size,maxit,model_effect,ploidy,dosage,max_geno_freq) {
  traits <- colnames(Y); d <- ncol(Y); models <- .sgwas_marker_models(model_effect,ploidy)
  multi <- .sgwas_fit_multitrait(Y,X,eig,maxit)
  observed <- lapply(seq_len(d),function(a) which(is.finite(Y[,a])))
  K <- if(any(!is.finite(Y))) tcrossprod(sweep(eig$vectors,2,eig$values,"*"),eig$vectors) else NULL
  marginal_eig <- lapply(observed,function(ix) if(length(ix)==nrow(Y)) eig else {
    eg <- eigen(K[ix,ix,drop=FALSE],symmetric=TRUE); eg$values <- pmax(eg$values,0); eg
  })
  marginal <- lapply(seq_len(d),function(a) {
    ix <- observed[[a]]; z <- .sgwas_fit_univariate(Y[ix,a],X[ix,,drop=FALSE],marginal_eig[[a]])
    z$n_samples <- length(ix); z$n_observed <- length(ix); z
  })
  names(marginal) <- traits; rows <- list()
  for(model in models) for(start in seq.int(1L,ncol(G),by=chunk_size)) {
    js <- seq.int(start,min(start+chunk_size-1L,ncol(G)))
    code <- if(model=="additive") list(G=G[,js,drop=FALSE],n0=rep(NA_integer_,length(js)),
                                       n1=rep(NA_integer_,length(js)),status=rep("ok",length(js))) else
                                         .sgwas_encode_dominance(dosage[,js,drop=FALSE],model,ploidy,max_geno_freq)
    gr <- crossprod(eig$vectors,code$G); z <- .sgwas_mt_marker_test(multi$fit,gr)
    base <- cbind(map[js,,drop=FALSE],Marker_Model=model,N_Group0=code$n0,N_Group1=code$n1,
                  Marker_Coding_Status=code$status)
    make_row <- function(test,trait,trait_index,beta,se,p,df_test,df,status,nobs_trait=NA_integer_) {
      bad <- code$status!="ok"; beta[bad] <- se[bad] <- p[bad] <- NA_real_
      status[bad] <- code$status[bad]
      cbind(base,Trait=trait,Trait_Index=trait_index,N_Traits=d,Test=test,Beta=beta,SE=se,
            Pvalue=p,DF_Test=df_test,DF_Residual=df,Status=status,
            N=if(test=="Trait_marginal") nobs_trait else nrow(Y),N_Cohort=nrow(Y),
            N_Observed_Trait=nobs_trait,N_Observed_Phenotypes=if(test=="Trait_marginal") nobs_trait else sum(is.finite(Y)),
            Inference=if(test=="Trait_marginal") "observed_trait_only" else multi$inference)
    }
    for(a in seq_len(d)) {
      gm <- if(length(observed[[a]])==nrow(Y)) gr else
        crossprod(marginal_eig[[a]]$vectors,code$G[observed[[a]],,drop=FALSE])
      single <- .sgwas_marker_test(marginal[[a]]$fit,gm)
      rows[[length(rows)+1L]] <- make_row("Trait_marginal",traits[a],a,single$Beta1,single$SE1,
                                          single$P1,1L,single$DF_Residual,single$Status,length(observed[[a]]))
      rows[[length(rows)+1L]] <- make_row("Trait_joint",traits[a],a,z$beta[,a]*multi$scales[a],
                                          z$se[,a]*multi$scales[a],z$p[,a],1L,z$df,z$status,length(observed[[a]]))
    }
    rows[[length(rows)+1L]] <- make_row("Multitrait_joint","All_selected_traits",0L,
                                        rep(NA_real_,length(js)),rep(NA_real_,length(js)),z$pjoint,d,z$df,z$status)
  }
  list(results=do.call(rbind,rows),fits=list(Multitrait=multi,Marginal=marginal))
}

# Adapt long-form multivariate results to the existing hypothesis-specific plots.
.sgwas_mt_report_data <- function(r) {
  key <- paste(r$Ploidy,r$Trait_Index,r$Marker_Model,r$SNP,sep="\r")
  take <- !duplicated(key); out <- r[take,,drop=FALSE]; out$Microbial_Trait <- "None"
  if(!"Inference" %in% names(out)) out$Inference <- "multitrait_phenotypes"
  for(test in unique(r$Test)) {
    tag <- switch(test,Trait_marginal="Trait",Trait_joint="Trait_joint",Multitrait_joint="Multitrait")
    source <- which(r$Test==test); ix <- source[match(key[take],key[source])]
    for(pair in list(c("P_","Pvalue"),c("Beta_","Beta"),c("SE_","SE")))
      out[[paste0(pair[1],tag)]] <- r[[pair[2]]][ix]
    out[[paste0("P_",tag,"_BH")]] <- r$FDR_BH[ix]
    out[[paste0("P_",tag,"_Bonferroni")]] <- r$Bonferroni[ix]
    for(field in intersect(c("Inference","N","N_Cohort","N_Observed_Trait","N_Observed_Phenotypes"),names(r)))
      out[[paste0("P_",tag,"_",field)]] <- r[[field]][ix]
    out[[paste0("PVE_",tag)]] <- NA_real_
    status_name <- switch(test,Trait_marginal="Status_FER",Trait_joint="Status_Bivariate",Multitrait_joint="Status_Multitrait")
    out[[status_name]] <- r$Status[ix]
  }
  out
}

# Thresholds use the genotype-matched, covariate-complete phenotype cohort BEFORE
# any phenotype exclusions. All-missing rows are counted but cannot be imputed.
.sgwas_mt_missing_policy <- function(Y,eligible,strategy,threshold) {
  missing <- !is.finite(Y); n <- sum(eligible); d <- ncol(Y)
  cells <- colSums(missing[eligible,,drop=FALSE]); complete <- rowSums(missing)==0L
  tab <- data.frame(Scope=c("Overall",rep("Trait",d)),Trait=c("All_selected_traits",colnames(Y)),
                    N_Eligible_Samples=n,N_Cells=c(n*d,rep(n,d)),N_Missing=c(sum(cells),cells),
                    N_Samples_With_Missing=c(sum(rowSums(missing[eligible,,drop=FALSE])>0L),cells))
  tab$Missing_Rate <- if(n) tab$N_Missing/tab$N_Cells else NA_real_
  tab$Missing_Percent <- 100*tab$Missing_Rate; tab$Threshold <- threshold
  tab$Sample_Missing_Rate <- if(n) tab$N_Samples_With_Missing/n else NA_real_
  tab$Sample_Missing_Percent <- 100*tab$Sample_Missing_Rate
  tab$Exceeds_Threshold <- is.finite(tab$Missing_Rate) & tab$Missing_Rate>threshold
  fallback <- strategy=="impute" && any(tab$Exceeds_Threshold)
  used <- if(fallback) "complete_case" else strategy
  keep <- eligible & if(used=="complete_case") complete else rowSums(!missing)>0L
  tab$Requested_Strategy <- strategy; tab$Applied_Strategy <- used; tab$Fallback <- fallback
  tab$N_Complete_Case_Samples <- sum(eligible & complete)
  tab$N_All_Traits_Missing <- sum(eligible & rowSums(!missing)==0L)
  tab$N_Analysed_Samples <- sum(keep)
  list(table=tab,keep=keep,strategy=used,fallback=fallback)
}

.sgwas_multitrait_gwas <- function(args,n_pcs,max_missing,chunk_size,maxit,max_geno_freq,output_style,plot_formats,
                                   missing_strategy="impute",phenotype_missing_threshold=.30) {
  traits <- .sgwas_names(args$trait_names); covs <- .sgwas_names(args$covariate_pheno)
  if(length(traits)<2L) stop("sgwas-multitrait requires at least two distinct trait_names from phenofile")
  if(any(traits %in% covs)) stop("A selected multi-trait outcome cannot also be a covariate")
  if(length(args$perm)!=1L || !is.numeric(args$perm) || !is.finite(args$perm) || !args$perm %in% c(0,1))
    stop("Multi-trait permutation thresholds are not implemented; use perm = 0 or 1")
  for(nm in c("LOCO","fdr","bonferroni"))
    if(length(args[[nm]])!=1L || !is.logical(args[[nm]]) || is.na(args[[nm]])) stop(nm," must be TRUE or FALSE")
  if(length(args$maf)!=1L || !is.numeric(args$maf) || !is.finite(args$maf) || args$maf<0 || args$maf>.5)
    stop("maf must be between 0 and 0.5")
  if(length(args$suggestive)!=1L || !is.numeric(args$suggestive) || !is.finite(args$suggestive) || args$suggestive<0)
    stop("suggestive must be a finite nonnegative number")
  ph <- .sgwas_phenotypes(args$phenofile)
  if(is.null(ph)) stop("sgwas-multitrait requires phenofile")
  absent <- setdiff(c(traits,covs),names(ph))
  if(length(absent)) stop("Columns absent from phenofile: ",paste(absent,collapse=", "))
  if("Plant_ID" %in% c(traits,covs)) stop("The sample identifier cannot be a trait or covariate")
  if(length(traits)>10L) warning("Unstructured multi-trait REML estimates t*(t+1) covariance parameters; consider a smaller prespecified trait set.",call.=FALSE)
  if(args$cores>1L) message("[systemsGWAS] Multi-trait scans use serial vectorised chunks; cores controls the legacy engine only.")
  outdir <- file.path(args$wdir,paste0(args$projname,"_multitrait_gwas")); dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  complete <- FALSE
  status <- function(value) .sgwas_write_tsv(data.frame(Status=value,Mode="sgwas-multitrait",Time=as.character(Sys.time())),file.path(outdir,"RUN_STATUS.tsv"))
  status("running"); on.exit(if(!complete) status("failed"),add=TRUE)
  ids_all <- unique(ph$Plant_ID)
  if(length(ids_all)<10L) stop("Too few phenotype samples for multi-trait GWAS")
  Y_all <- vapply(traits,function(t) .sgwas_id_values(ph,t,ids_all),numeric(length(ids_all)))
  colnames(Y_all) <- traits
  C_all <- matrix(numeric(),length(ids_all),0L,dimnames=list(NULL,character()))
  if(length(covs)) {C_all <- vapply(covs,function(t) .sgwas_id_values(ph,t,ids_all),numeric(length(ids_all))); colnames(C_all) <- covs}
  results <- audits <- nulls <- covariance <- attempts <- allfits <- missingness <- list()
  for(ploidy in unique(as.numeric(args$ploidy_levels))) {
    pkey <- paste0(ploidy,"x"); folder <- file.path(outdir,pkey); dir.create(folder,showWarnings=FALSE)
    raw <- .sgwas_read_dosage(args[[paste0("genofile_",ploidy,"x")]],ploidy)
    complete_y <- rowSums(!is.finite(Y_all))==0L; complete_c <- rowSums(!is.finite(C_all))==0L
    has_g <- ids_all %in% rownames(raw$G)
    policy <- .sgwas_mt_missing_policy(Y_all,complete_c & has_g,missing_strategy,phenotype_missing_threshold)
    keep <- policy$keep
    missingness[[pkey]] <- cbind(Ploidy=ploidy,policy$table)
    .sgwas_write_tsv(do.call(rbind,missingness),file.path(outdir,"MULTITRAIT_MISSINGNESS.tsv"))
    message("[systemsGWAS] Phenotype missingness (",pkey,"; ",sum(complete_c & has_g)," eligible samples): ",
            paste(paste0(policy$table$Trait,"=",round(policy$table$Missing_Percent,2),"%"),collapse="; "),
            "; applied strategy=",policy$strategy)
    if(policy$fallback) {
      flagged <- policy$table[policy$table$Exceeds_Threshold,,drop=FALSE]
      warning("[systemsGWAS] Phenotype missingness exceeds ",100*phenotype_missing_threshold,"% at ",pkey," (",
              paste(paste0(flagged$Trait,"=",round(flagged$Missing_Percent,2),"%"),collapse=", "),
              "); reverting this ploidy to complete-case multi-trait analysis.",call.=FALSE)
    }
    missing_traits <- apply(!is.finite(Y_all),1,function(z) paste(traits[z],collapse=";"))
    observed_count <- rowSums(is.finite(Y_all))
    reason <- ifelse(!has_g,"no_genotype",ifelse(!complete_c,"missing_covariate",
                                                 ifelse(observed_count==0L,"all_traits_missing",ifelse(!keep,"incomplete_traits","included"))))
    audits[[pkey]] <- data.frame(Ploidy=ploidy,Plant_ID=ids_all,Complete_Traits=complete_y,
                                 Missing_Traits=missing_traits,N_Observed_Traits=observed_count,
                                 Missing_Rate=1-observed_count/length(traits),Complete_Covariates=complete_c,
                                 Has_Genotype=has_g,Included=keep,Exclusion_Reason=reason,
                                 Applied_Strategy=policy$strategy,
                                 Imputed_Traits=ifelse(keep & policy$strategy=="impute",missing_traits,""))
    .sgwas_write_tsv(do.call(rbind,audits),file.path(outdir,"SAMPLE_AUDIT.tsv"))
    ids <- ids_all[keep]; Y <- Y_all[keep,,drop=FALSE]; C <- C_all[keep,,drop=FALSE]
    if(length(ids)<max(10L,ncol(C)+n_pcs+length(traits)+2L))
      stop("Too few eligible samples after ",policy$strategy," filtering for multi-trait GWAS at ",pkey)
    if(ncol(C)) {
      cs <- apply(C,2,sd)
      if(any(!is.finite(cs)) || any(cs<=0)) stop("Constant or invalid phenotype covariate on the multi-trait cohort")
      C <- scale(C)
    }
    qc <- .sgwas_dosage_qc(raw,ids,ploidy,args$maf,max_missing,max_geno_freq)
    .sgwas_write_tsv(qc$ledger,file.path(folder,"MARKER_QC.tsv"))
    if("Dom" %in% args$model_effect && any(qc$map$N_Fractional>0))
      warning("Dominance encodings round fractional dosages; genotype uncertainty is not propagated.",call.=FALSE)
    if("Dom" %in% args$model_effect && all(qc$map$N_Het_calls==0) && !any(qc$map$N_Fractional>0))
      message("[systemsGWAS] Homozygous cohort: dominance contrasts duplicate additive evidence.")
    af <- colMeans(qc$G)/ploidy; denom <- ploidy*af*(1-af); cp <- tcrossprod(qc$centered)
    full <- .sgwas_kinship_eigen(cp,sum(denom)); X <- cbind(Intercept=1,C)
    if(n_pcs) {
      # The legacy shift of K creates a constant eigenvector. Genotype PCs must
      # instead come from the centred marker Gram matrix, before that shift.
      pc_eig <- eigen(cp,symmetric=TRUE)
      if(n_pcs>sum(pc_eig$values>max(pc_eig$values)*1e-8)) stop("Too few nonzero genomic PCs")
      pcs <- pc_eig$vectors[,seq_len(n_pcs),drop=FALSE]; colnames(pcs) <- paste0("Genotype_PC",seq_len(n_pcs)); X <- cbind(X,pcs)
    }
    if(qr(X)$rank<ncol(X)) stop("Collinear covariates/PCs in multi-trait GWAS")
    .sgwas_write_tsv(data.frame(Plant_ID=ids,Y,check.names=FALSE),file.path(folder,"ANALYSED_TRAITS.tsv"))
    .sgwas_write_tsv(data.frame(Plant_ID=ids,X,check.names=FALSE),file.path(folder,"ANALYSED_COVARIATES.tsv"))
    groups <- if(args$LOCO) unique(qc$map$CHROM) else "All"
    if(args$LOCO && length(groups)<2L) stop("LOCO requires retained markers on at least two chromosomes")
    allfits[[pkey]] <- list(); imputation_manifest <- list()
    for(ch in groups) {
      selected <- if(args$LOCO) which(qc$map$CHROM==ch) else seq_len(ncol(qc$G))
      eg <- if(args$LOCO) .sgwas_kinship_eigen(cp-tcrossprod(qc$centered[,selected,drop=FALSE]),sum(denom[-selected])) else full
      message("[systemsGWAS] Multi-trait (",paste(traits,collapse=", "),"); ",pkey,"; ",length(ids)," samples; ",
              sum(!is.finite(Y))," missing phenotype cells; strategy=",policy$strategy,"; ",
              length(selected)," markers; ",if(args$LOCO) paste0("LOCO ",ch) else "whole-genome kinship")
      fit <- tryCatch(.sgwas_multitrait_scan(Y,X,qc$centered[,selected,drop=FALSE],qc$map[selected,,drop=FALSE],eg,
                                             chunk_size,maxit,args$model_effect,ploidy,qc$G[,selected,drop=FALSE],max_geno_freq),error=function(e) {
                                               .sgwas_write_tsv(data.frame(Ploidy=ploidy,Chromosome=ch,Error=conditionMessage(e)),file.path(outdir,"FAILED.tsv"))
                                               if(!is.null(e$attempts)) .sgwas_write_tsv(e$attempts,file.path(folder,"FAILED_REML_ATTEMPTS.tsv"))
                                               stop(e)
                                             })
      results[[length(results)+1L]] <- cbind(Ploidy=ploidy,Kinship_Chromosome=ch,fit$results)
      allfits[[pkey]][[ch]] <- fit$fits
      if(file.exists(file.path(folder,"FAILED_REML_ATTEMPTS.tsv")))
        unlink(file.path(folder,"FAILED_REML_ATTEMPTS.tsv"))
      f <- fit$fits$Multitrait
      # Keep observed values and masks distinct from completed prediction exports.
      imp_folder <- file.path(folder,"imputation",if(args$LOCO) paste0("LOCO_",match(ch,groups)) else "All")
      dir.create(imp_folder,recursive=TRUE,showWarnings=FALSE)
      .sgwas_write_tsv(data.frame(Plant_ID=ids,f$completed,check.names=FALSE),file.path(imp_folder,"COMPLETED_TRAITS.tsv"))
      .sgwas_write_tsv(data.frame(Plant_ID=ids,!is.finite(Y),check.names=FALSE),file.path(imp_folder,"IMPUTATION_MASK.tsv"))
      imp <- f$imputation
      .sgwas_write_tsv(cbind(data.frame(Ploidy=rep(ploidy,nrow(imp)),Kinship_Chromosome=rep(ch,nrow(imp)),
                                        Plant_ID=ids[imp$Sample_Index]),imp),file.path(imp_folder,"IMPUTED_VALUES.tsv"))
      imputation_manifest[[length(imputation_manifest)+1L]] <- data.frame(Kinship_Chromosome=ch,
                                                                          Directory=file.path("imputation",basename(imp_folder)),N_Imputed=nrow(imp),Inference=f$inference)
      .sgwas_write_tsv(do.call(rbind,imputation_manifest),file.path(folder,"IMPUTATION_MANIFEST.tsv"))
      attempts[[length(attempts)+1L]] <- cbind(Ploidy=ploidy,Chromosome=ch,f$attempts)
      fs <- c(list(f),fit$fits$Marginal); labels <- c("Multitrait",traits)
      for(i in seq_along(fs)) {
        a <- fs[[i]]
        nulls[[length(nulls)+1L]] <- data.frame(Ploidy=ploidy,Chromosome=ch,Model=if(i==1L) "Multitrait" else "Marginal",
                                                Trait=if(i==1L) "All_selected_traits" else labels[i],N=a$n_samples,N_Traits=if(i==1L) length(traits) else 1L,
                                                N_Observed_Phenotypes=a$n_observed,N_Missing_Phenotypes=if(i==1L) a$n_missing else 0L,
                                                Converged=a$converged,Boundary=a$boundary,Optimizer=a$optimizer,Objective=a$objective)
      }
      for(component in c("Genetic","Residual")) {
        V <- if(component=="Genetic") f$G_original else f$R_original
        correlation <- cov2cor(V)
        if(f$boundary && component=="Genetic") correlation[] <- NA_real_
        grid <- expand.grid(Trait1=traits,Trait2=traits,stringsAsFactors=FALSE)
        covariance[[length(covariance)+1L]] <- cbind(Ploidy=ploidy,Chromosome=ch,Component=component,grid,
                                                     Covariance=as.vector(V),Correlation=as.vector(correlation),Boundary=f$boundary)
      }
    }
  }
  result <- do.call(rbind,results); rownames(result) <- NULL
  result$FDR_BH <- result$Bonferroni <- NA_real_
  # Three declared families. The omnibus occurs once per marker/model/ploidy.
  for(ix in split(seq_len(nrow(result)),result$Test)) {
    good <- ix[is.finite(result$Pvalue[ix])]
    result$FDR_BH[good] <- p.adjust(result$Pvalue[good],"BH")
    result$Bonferroni[good] <- p.adjust(result$Pvalue[good],"bonferroni")
  }
  summary <- do.call(rbind,nulls); cov_table <- do.call(rbind,covariance); attempt_table <- do.call(rbind,attempts)
  .sgwas_write_tsv(result,file.path(outdir,"ASSOCIATIONS.tsv"))
  .sgwas_write_tsv(summary,file.path(outdir,"NULL_MODELS.tsv"))
  .sgwas_write_tsv(cov_table,file.path(outdir,"TRAIT_COVARIANCES.tsv"))
  .sgwas_write_tsv(attempt_table,file.path(outdir,"REML_ATTEMPTS.tsv"))
  .sgwas_write_tsv(data.frame(Trait_Index=seq_along(traits),Trait=traits),file.path(outdir,"TRAIT_SET.tsv"))
  adapted <- .sgwas_mt_report_data(result); manifest <- list()
  total_tests <- vapply(c("P_Trait","P_Trait_joint","P_Multitrait"),function(p) sum(is.finite(adapted[[p]])),integer(1))
  for(ploidy in unique(result$Ploidy)) for(t in unique(result$Trait_Index)) {
    data <- result[result$Ploidy==ploidy & result$Trait_Index==t,,drop=FALSE]
    subdir <- if(t==0L) "Omnibus" else paste0("trait_",sprintf("%02d",t),"_",gsub("[^A-Za-z0-9_.-]","_",traits[t]))
    folder <- file.path(outdir,paste0(ploidy,"x"),subdir); dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    .sgwas_write_tsv(data,file.path(folder,"ASSOCIATIONS.tsv")); r <- .sgwas_mt_report_data(data)
    if(output_style=="verbose") .sgwas_gwas_results_reports(r,folder,adapted,args$fdr,args$bonferroni,args$suggestive,plot_formats,args$fdr_level)
    if("pdf" %in% plot_formats) .sgwas_association_plots(r,file.path(folder,"Manhattan_QQ.pdf"),args$fdr,args$bonferroni,args$suggestive,total_tests,fdr_level=args$fdr_level)
    if(output_style=="compact" && "tiff" %in% plot_formats) .sgwas_association_plots(r,file.path(folder,"Manhattan_QQ_%03d.tiff"),args$fdr,args$bonferroni,args$suggestive,total_tests,"tiff",args$fdr_level)
    manifest[[length(manifest)+1L]] <- data.frame(Ploidy=ploidy,Trait_Index=t,Trait=if(t==0L) "All_selected_traits" else traits[t],
                                                  Directory=file.path(paste0(ploidy,"x"),subdir))
  }
  .sgwas_write_tsv(do.call(rbind,manifest),file.path(outdir,"REPORT_MANIFEST.tsv"))
  metadata <- list(mode="sgwas-multitrait",arguments=args,traits=traits,
                   missing_strategy=missing_strategy,phenotype_missing_threshold=phenotype_missing_threshold,
                   missingness=do.call(rbind,missingness),
                   missing="Conditional null-model imputation with observed-data REML/marker information; complete-case fallback if overall or any trait exceeds threshold; marginal tests use observed trait values only",
                   imputation_uncertainty="Missing cells integrated out in association tests; export SEs are plug-in prediction errors including fixed-effect estimation, excluding covariance-parameter uncertainty",
                   n_pcs=n_pcs,pc_source="full-genome centred dosage Gram matrix before kinship shifting; also in LOCO",
                   max_missing=max_missing,max_geno_freq=max_geno_freq,chunk_size=chunk_size,maxit=maxit,
                   output_style=output_style,plot_formats=plot_formats,covariance="unstructured genetic and residual; trait-specific fixed coefficients",
                   method="Spectral REML null covariance; P3D null residual scale; approximate F tests, not marker-refitted likelihood-ratio tests",
                   multiple_testing="BH FDR and Bonferroni across all finite marker/model/trait/ploidy tests in each of Trait_marginal, Trait_joint, Multitrait_joint separately",
                   models=allfits,session=sessionInfo())
  saveRDS(metadata,file.path(outdir,"MODEL_DETAILS.rds"))
  if(file.exists(file.path(outdir,"FAILED.tsv"))) unlink(file.path(outdir,"FAILED.tsv"))
  complete <- TRUE; status("complete")
  message("[systemsGWAS] Multi-trait results: ",file.path(outdir,"ASSOCIATIONS.tsv"))
  invisible(list(results=result,null_models=summary,trait_covariances=cov_table,reml_attempts=attempt_table,
                 sample_audit=do.call(rbind,audits),missingness=do.call(rbind,missingness),metadata=metadata,outdir=outdir))
}


# Public entry point. Legacy named arguments remain accepted; prefer
# microbial_features + microbial_method for new microbial-phenotype analyses.
systemsGWAS <- function(
    wdir="./", projname="GWAS", ploidy_levels=c(2,4,6,8),
    trait_names=c("trait1","trait2"), alternate_trait=NULL,
    trait_microbial_proxy=c("NULL"), model_effect=c("Add","Dom"),
    fdr=TRUE, fdr_level=0.05, bonferroni=TRUE, suggestive=5, perm=1, cores=1,
    genofile_2x=NULL, genofile_4x=NULL, genofile_6x=NULL, genofile_8x=NULL,
    phenofile=NULL, method=c("GLM","MLM"), covariate_pheno=NULL,
    covariate_metag=FALSE, maf=0.02, LOCO=FALSE,
    pheno_taxa_strain=NULL, pheno_taxa_species=NULL,
    metag_data_strains=NULL, metag_data_species=NULL,
    microbial_features=NULL, microbial_method="auto", microbial_response=NULL,
    microbial_weights=NULL, microbial_trait_name="Microbial_composite",
    microbial_transform="auto", microbial_prevalence=0.05,
    microbial_bootstraps=100L, microbial_stability=0.30,
    microbial_seed=123L, microbial_pseudocount=1e-6, run_gwas=TRUE,
    association_mode="sgwas-marginal", association_n_pcs=0L,
    association_max_missing=0.50, association_chunk_size=256L,
    association_maxit=300L, association_geno_freq=0.99,
    association_output_style="verbose", association_plot_formats=c("tiff"),
    multitrait_missing="impute", multitrait_missing_threshold=0.30) {
  new_mode <- length(.sgwas_names(microbial_features)) > 0L || length(.sgwas_names(microbial_weights)) > 0L
  if (missing(association_mode) && new_mode) {
    association_mode <- "sgwas-all"
  }
  association_mode <- match.arg(
    tolower(association_mode),
    c("sgwas-marginal", "gwaspoly-marginal", "sgwas-joint", "sgwas-interaction", "sgwas-all", "sgwas-multitrait")
  )
  integrated <- association_mode != "gwaspoly-marginal"
  multitrait <- identical(association_mode,"sgwas-multitrait")
  if(multitrait) {
    multitrait_missing <- match.arg(tolower(multitrait_missing),c("impute","complete_case"))
    if(length(multitrait_missing_threshold)!=1L || !is.numeric(multitrait_missing_threshold) ||
       !is.finite(multitrait_missing_threshold) || multitrait_missing_threshold<0 || multitrait_missing_threshold>1)
      stop("multitrait_missing_threshold must be between 0 and 1")
    if(length(.sgwas_names(trait_names))<2L)
      stop("sgwas-multitrait requires at least two distinct trait_names from phenofile")
    if(new_mode || length(.sgwas_names(trait_microbial_proxy)) || length(.sgwas_names(alternate_trait)) ||
       length(.sgwas_names(pheno_taxa_species)) || length(.sgwas_names(pheno_taxa_strain)) || !identical(covariate_metag,FALSE))
      stop("sgwas-multitrait reads outcomes from phenofile; omit microbial/proxy selectors and use covariate_metag = FALSE")
  }
  if (association_mode %in% c("sgwas-joint", "sgwas-interaction", "sgwas-all") && !new_mode)
    stop("sgwas-joint/sgwas-interaction/sgwas-all require microbial_features or microbial_weights")
  if (association_mode == "sgwas-marginal" && new_mode)
    stop("association_mode = 'sgwas-marginal' tests phenotypes without microbiome data; omit microbial_features/microbial_weights")
  if(integrated) {
    if(missing(method)) method <- "MLM"
    if(missing(model_effect)) model_effect <- "Add"
    .sgwas_validate_integrated_models(method,model_effect)
    association_n_pcs <- .sgwas_integer(association_n_pcs,"association_n_pcs",0L)
    association_chunk_size <- .sgwas_integer(association_chunk_size,"association_chunk_size")
    association_maxit <- .sgwas_integer(association_maxit,"association_maxit")
    association_output_style <- match.arg(association_output_style,c("verbose","compact"))
    association_plot_formats <- tolower(.sgwas_names(association_plot_formats))
    if(length(association_geno_freq)!=1L || !is.numeric(association_geno_freq) ||
       !is.finite(association_geno_freq) || association_geno_freq<=0.5 || association_geno_freq>1)
      stop("association_geno_freq must be between 0.5 and 1")
    if(any(!association_plot_formats %in% c("pdf","tiff")))
      stop("association_plot_formats must contain pdf, tiff, both, or NULL to skip plots")
    if(association_output_style=="verbose" && length(association_plot_formats) && !requireNamespace("ggplot2",quietly=TRUE))
      stop("Install ggplot2 for verbose integrated plots, or use association_output_style = 'compact'")
    if(length(association_max_missing)!=1L || !is.numeric(association_max_missing) ||
       !is.finite(association_max_missing) || association_max_missing<0 || association_max_missing>1)
      stop("association_max_missing must be between 0 and 1")
  }
  if(length(fdr_level)!=1L || !is.numeric(fdr_level) || !is.finite(fdr_level) || fdr_level<=0 || fdr_level>1)
    stop("fdr_level must be between 0 and 1")
  if(length(run_gwas)!=1L || !is.logical(run_gwas) || is.na(run_gwas))
    stop("run_gwas must be TRUE or FALSE")
  wdir <- normalizePath(wdir,mustWork=TRUE)
  if(length(projname)!=1L || is.na(projname) || !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$",projname))
    stop("projname must be a simple project prefix (letters, digits, _, . or -)")
  numeric_cores <- suppressWarnings(as.numeric(cores))
  if(length(numeric_cores)!=1L || !is.finite(numeric_cores) || numeric_cores<1 || numeric_cores!=floor(numeric_cores))
    stop("cores must be a positive integer")
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
    is_weighted <- !is.null(.sgwas_names(microbial_weights)) ||
      identical(tolower(trimws(microbial_method)), "weighted")
    if (integrated && is_weighted) {
      outdir <- file.path(wdir, paste0(projname, "_integrated_gwas"), paste0(projname, "_microbial_traits"))
    } else {
      outdir <- file.path(wdir, paste0(projname, "_microbial_traits"))
    }
    
    prepared <- systemsGWAS_prepare_microbiome(
      microbiome_file=.sgwas_path(metag_data_species,wdir), features=microbial_features,
      method=microbial_method, response=microbial_response,
      phenofile=if(needs_pheno) .sgwas_training_phenofile(phenofile,wdir,tolower(microbial_method)=="spls") else NULL,
      weights_file=.sgwas_path(microbial_weights,wdir), trait_name=microbial_trait_name,
      transform=microbial_transform, prevalence=microbial_prevalence,
      bootstraps=microbial_bootstraps, stability=microbial_stability,
      seed=microbial_seed, pseudocount=microbial_pseudocount,
      outdir=outdir, covariates=covariate_pheno)
    if(!run_gwas) return(invisible(prepared))
    if(!integrated) {
      args <- .sgwas_standard_proxy_args(args,prepared)
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
  gwas_result <- if(multitrait) .sgwas_multitrait_gwas(args,association_n_pcs,association_max_missing,
                                                       association_chunk_size,association_maxit,association_geno_freq,
                                                       association_output_style,association_plot_formats,
                                                       multitrait_missing,multitrait_missing_threshold) else
                                                         if(integrated) .sgwas_integrated_gwas(args,prepared,association_mode,
                                                                                               association_n_pcs,association_max_missing,association_chunk_size,association_maxit,association_geno_freq,
                                                                                               association_output_style,association_plot_formats) else
                                                                                                 do.call(.systemsGWAS_engine,args)
  if(!is.null(prepared)) {prepared$gwas_result <- gwas_result;return(invisible(prepared))}
  invisible(gwas_result)
}
