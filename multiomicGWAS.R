#!/usr/bin/env Rscript

multiomicGWAS <- function(
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
  setwd(wdir)
  dir_prefix_name <- projname
  pop <- projname
  ploidy_levels <- as.numeric(ploidy_levels)
  gwas_method <- strsplit(method, ",")
  maf <- as.numeric(maf)
  file_ploidy_2 <- genofile_2x
  file_ploidy_4 <- genofile_4x
  file_ploidy_6 <- genofile_6x
  file_ploidy_8 <- genofile_8x
  phenotype_data <- phenofile
  alternate_trait=alternate_trait
  if (is.null(covariate_pheno) || covariate_pheno == "NULL") {
    covariatename <- NULL
  } else {
    covariatename <- unlist(strsplit(covariate_pheno, ","))
  }
  if (is.null(trait_microbial_proxy) || trait_microbial_proxy == "NULL") {
    trait_microbial_proxy <- NULL
  } else {
    trait_microbial_proxy <- unlist(strsplit(trait_microbial_proxy, ","))
  }
  covariate <- covariate_metag
  corr_coeff <- if(covariate_metag) "full" else NULL
  # Microbiome as phenotypes
  taxa_strain <- pheno_taxa_strain
  taxa_species <- pheno_taxa_species
  metagenome_data_strains <- metag_data_strains
  metagenome_data_species <- metag_data_species
  if (!is.null(model_effect)) {
    model_effect <- strsplit(model_effect, ",")
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
  cores=1                                                 # If running on Linux machine, increase number of cores to increase computational speed
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
        if(!is.null(metagenome_data_species)) {metag <- read.table(metagenome_data_species, header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)}
        if(!is.null(metagenome_data_species)) {metag <- read.table(metagenome_data_species, header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)}

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
          geno <- read.table(paste("../",genotype_data,sep=""), header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
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
            traits[, 2:ncol(traits)] <- sapply(traits[, 2:ncol(traits)], as.factor)
            samples <- as.data.frame(colnames(geno[,-c(1:5)]))
            colnames(samples) <- "Plant_ID"
            traits <- merge(traits, samples, by="Plant_ID")
          }
          if (!is.null(metagenome_data)){traits_metag <- traits}

          if (!is.null(phenotype_data)){
            traits_pheno <- read.delim(paste("../",phenotype_data,sep=""), header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
            colnames(traits_pheno)[1] <- "Plant_ID"
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
          names(traits) <- sub(".x$", "", names(traits))
          names(traits) <- sub(".y$", "", names(traits))
          traits_hold <- traits[,c(1,2)]; traits_hold[,2] <- as.numeric(traits_hold[,2])
          traits_hold <- aggregate(traits_hold[,2],by=list(name=traits_hold$Plant_ID),data=traits_hold,FUN=mean)
          colnames(traits_hold) <- colnames(traits)[c(1,2)]
          if(ncol(traits) > 2){
            for(i in c(3:ncol(traits))){
              traits_holdtrans <- traits[,c(1,i)]; traits_holdtrans[,2] <- as.numeric(traits_holdtrans[,2])
              traits_holdtrans <- aggregate(traits_holdtrans[,2],by=list(name=traits_holdtrans$Plant_ID),data=traits_holdtrans,FUN=mean)
              colnames(traits_holdtrans) <- colnames(traits)[c(1,i)]
              traits_hold <- merge(traits_hold, traits_holdtrans, by="Plant_ID")
            }
          }
          traits <- traits_hold
          names(traits) <- sub(".x$", "", names(traits))
          names(traits) <- sub(".y$", "", names(traits))
          traits <- traits[, !duplicated(colnames(traits))]
          write.csv(traits, "traits.csv", row.names=F, quote = FALSE)


          #############################################################################################################################################################################
          # Perform Genome-Wide Association Analysis
          ##########################################
          traits <- read.csv("traits.csv", header=T, sep=",", check.names=FALSE,stringsAsFactors=FALSE)

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
              if(trait_microbial_proxy[1] == "auto"){
                set.seed(123)
                X <- metag_proxy[, colnames(metag_proxy) != train_traitname, drop = FALSE]
                Y <- metag_proxy[, train_traitname, drop = FALSE]
                test.keepX <- c(1,2,3,4,5,10,20,30,40,50,100,200,300,400,500)
                test.keepX <- test.keepX[test.keepX <= ncol(X)]
                max_ncomp <- min(5, ncol(X))
                folds <- sample(rep(1:5, length.out = nrow(X)))
                tuning_results <- list()
                for(ncomp in 1:max_ncomp){
                  for(k in test.keepX){
                    cv_cor <- c()
                    for(f in sort(unique(folds))){
                      train_idx <- folds != f
                      test_idx  <- folds == f
                      fit <- try(mixOmics::spls(
                          X = X[train_idx, , drop = FALSE],
                          Y = Y[train_idx, , drop = FALSE],
                          ncomp = ncomp,
                          keepX = rep(k, ncomp)),silent = TRUE)
                      if(inherits(fit, "try-error"))
                        next
                      pred <- try(predict(fit, X[test_idx, , drop = FALSE])$predict[, , ncomp],silent = TRUE)
                      if(inherits(pred, "try-error"))
                        next
                      pred <- as.numeric(pred)
                      obs <- as.numeric(Y[test_idx, 1])
                      r <- suppressWarnings(cor(pred, obs, method = "spearman", use = "complete.obs"))
                      cv_cor <- c(cv_cor, r)
                    }
                    tuning_results[[length(tuning_results)+1]] <- data.frame(ncomp = ncomp, keepX = k,CV_R = mean(cv_cor, na.rm = TRUE))
                  }
                }
                tuning_results <- do.call(rbind, tuning_results)
                best_row <- tuning_results[which.max(tuning_results$CV_R),]
                best_ncomp <- best_row$ncomp
                best_keepX <- rep(best_row$keepX, best_ncomp)
                cat("\nSelected:", "ncomp =", best_ncomp, "keepX =", best_row$keepX, "CV_R =", round(best_row$CV_R, 3), "\n")

                spls_model <- mixOmics::spls(X = X, Y = Y, ncomp = best_ncomp, keepX = best_keepX)
                Y_vec <- as.numeric(Y[,1])
                for(h in seq_len(best_ncomp)){
                  r <- cor(spls_model$variates$X[,h], Y_vec, method = "spearman", use = "complete.obs")
                  if(!is.na(r) && r < 0){
                    spls_model$variates$X[,h] <- -spls_model$variates$X[,h]
                    spls_model$loadings$X[,h] <- -spls_model$loadings$X[,h]
                    if(!is.null(spls_model$loadings$Y))
                      spls_model$loadings$Y[,h] <- -spls_model$loadings$Y[,h]
                  }
                }
                pheno <- as.data.frame(spls_model$variates$X)
                selected_taxa_list <- lapply(seq_len(best_ncomp),
                  function(h){tmp <- selectVar(spls_model,comp = h)$X$value
                    data.frame(proxy_trait = rownames(tmp), weight = tmp[,1], comp = h,row.names = NULL
                    )
                  }
                )
                selected_taxa_weight <- do.call(rbind, selected_taxa_list)
                select_proxy_taxa <- unique(selected_taxa_weight$proxy_trait)

                if(length(select_proxy_taxa) == 0){stop(paste0(train_traitname, "\tFAILED: no associated taxa selected\n"))}
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
                if(length(overlapping_taxa) > 1){
                  X <- metag_proxy[,overlapping_taxa, drop = FALSE]
                  Y <- metag_proxy[, train_traitname , drop = FALSE]
                  spls_model <- spls(X = X, Y = Y, ncomp = 1, keepX = ncol(X))
                  proxy_trait <- spls_model$variates$X[, 1]
                  Y_vec <- as.numeric(Y[, 1])
                  r <- cor(proxy_trait, Y_vec, method = "spearman", use = "complete.obs")
                  if (!is.na(r) && r < 0) {
                    proxy_trait <- -proxy_trait
                    spls_model$variates$X[, 1] <- -spls_model$variates$X[, 1]
                    spls_model$loadings$X[, 1] <- -spls_model$loadings$X[, 1]
                  }
                  pheno <- data.frame(proxy_trait = proxy_trait)
                  selected_taxa_weight <- as.data.frame(selectVar(spls_model, comp = 1)$X$value)
                  select_proxy_taxa <- rownames(selected_taxa_weight)
                  selected_taxa_weight <- as.data.frame(cbind(proxy_trait = rownames(selected_taxa_weight), weight = selected_taxa_weight[,1]))
                }
              }
              shared_samples <- intersect(traits[,1], row.names(pheno))
              pheno <- pheno[row.names(pheno) %in% shared_samples, , drop = FALSE]
              pheno <- as.data.frame(pheno)
              pheno <- pheno[, colnames(pheno) != train_traitname , drop = FALSE]
              head(pheno)

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
                tcorr <- read.table(paste("../",corr_coeff,sep=""), header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
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

            geno <- read.table(paste("../",genotype_data,sep=""), header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
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

              # Calculate PVE
              geno_mat_imputed <- apply(geno_mat, 2, function(col) {
                col[is.na(col)] <- mean(col, na.rm = TRUE)
                return(col)
              })
              geno_mat_imputed <- as.matrix(geno_mat_imputed)
              max_dom <- ploidy / 2  # maximum dominance level
              # compute SNP variance (dosage variance)
              varX <- apply(geno_mat_imputed, 2, var, na.rm = TRUE)
              # compute additive PVE
              snps <- GWAS_scores_effects$SNP
              vX <- varX[match(snps, colnames(geno_mat_imputed))]
              GWAS_scores_effects$additive_PVE <- pmin((GWAS_scores_effects$additive_effects^2 * vX) / var_y, 1)
              # compute dominance PVE
              for(d in 1:max_dom){
                alt_col <- paste0(d, "-dom-alt_effects")
                ref_col <- paste0(d, "-dom-ref_effects")
                # alt dominance PVE
                if(alt_col %in% colnames(GWAS_scores_effects)) {
                  snps <- GWAS_scores_effects$SNP
                  # use corresponding genotype variance
                  vX <- varX[match(snps, colnames(geno_mat_imputed))]
                  GWAS_scores_effects[[paste0(d,"-dom-alt_PVE")]] <- pmin((GWAS_scores_effects[[alt_col]]^2 * vX) / var_y, 1)
                }
                # ref dominance PVE
                if(ref_col %in% colnames(GWAS_scores_effects)) {
                  snps <- GWAS_scores_effects$SNP
                  vX <- varX[match(snps, colnames(geno_mat_imputed))]
                  GWAS_scores_effects[[paste0(d,"-dom-ref_PVE")]] <- pmin((GWAS_scores_effects[[ref_col]]^2 * vX) / var_y, 1)
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
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape=22) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(fill = guide_legend(title = "Thresholds: ", nrow=1, keywidth = 3, keyheight = 3)) +
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
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape=22) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(fill = guide_legend(title = "Thresholds: ", nrow=1, keywidth = 3, keyheight = 3)) +
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
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape=22) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(fill = guide_legend(title = "Thresholds: ", nrow=1, keywidth = 3, keyheight = 3)) +
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
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape=22) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(fill = guide_legend(title = "Thresholds: ", nrow=1, keywidth = 3, keyheight = 3)) +
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
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape=22) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(fill = guide_legend(title = "Thresholds: ", nrow=1, keywidth = 3, keyheight = 3)) +
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
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape=22) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(fill = guide_legend(title = "Thresholds: ", nrow=1, keywidth = 3, keyheight = 3)) +
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
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape=22) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(fill = guide_legend(title = "Thresholds: ", nrow=1, keywidth = 3, keyheight = 3)) +
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
                    geom_point(data = legend_df, aes(x = x, y = y, fill = threshold), inherit.aes = FALSE, size = 10, shape=22) +
                    scale_fill_manual(values = setNames(legend_df$color, legend_df$threshold)) +
                    guides(fill = guide_legend(title = "Thresholds: ", nrow=1, keywidth = 3, keyheight = 3)) +
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
