#' Intensity Imputation for Proteomics Data
#'
#' This function performs imputation of missing intensity values in a proteomics data matrix. 
#' It uses PhosR's imputation methods or falls back to median imputation if PhosR's method fails. 
#' The function works on intensity matrices for both proteome and phosphoproteome data.
#' If select `pcaMethods` need to provide the normalized data from `normalization_ProTN`.
#'
#' @param proteome_data List; a list containing proteomics data. It must include either `c_anno` (sample annotations for proteome) or `c_anno_proteome` and `c_anno_phospho` (sample annotations for proteome and phosphoproteome, respectively). It also requires `psm_log_prot_df` and `psm_log_pet_df` (log-transformed intensity matrices for proteome and phosphoproteome).
#' @param type character; Type of imputation selected (`phosr`, `gaussian`, `missForest`, `pcaMethods`). Default: `phosr`. NOTE: pcaMethods require data after normalization.
#'
#' @return A list with the imputed proteomics data, including the updated `psm_log_prot_df` and `psm_log_pet_df` matrices.
#'
#' @examples
#' \dontrun{
#' proteome_data_imputed <- impute_intensity(proteome_data)
#' }
#'
#' @import data.table
#' @import PhosR
#' @import missForest
#' @import doParallel
#' @import parallel
#' @import pcaMethods
#' @export
impute_intensity <- function(proteome_data, type="phosr") {
  message("Imputation of intensity matrix in progress...")
  
  if(("c_anno" %in% names(proteome_data))){
    ## PROTEINS
    if(type=="pcaMethods"){
      # Add missing values
      if("dat_gene" %in% names(proteome_data)){
        samples_pep <- data.table("sample" = proteome_data$c_anno$sample, 
                                  "names_sample" =make.names(proteome_data$c_anno$sample))
        numeric_df <- proteome_data$dat_gene[, .(GeneName, is.na(.SD)), .SDcols = samples_pep$sample]
        numeric_df <- numeric_df[, .(GeneName,
                                     protein_impute_level = rowSums(.SD) / length(samples_pep$names_sample) * 100), .SDcols = samples_pep$names_sample]
        psm_prot <- copy(proteome_data$psm_anno_df)
        proteome_data$psm_anno_df <- psm_prot[numeric_df, on= .(symbol = GeneName)]
        
        proteome_data$dat_gene <- impute_matrix(proteome_data$dat_gene, proteome_data$c_anno, type="pcaMethods_protein")
      } else{
        stop("'pcaMethods' selected, please provide the normalized data after 'normalization_ProTN'")
      }
    } else{
      # Add missing values
      # Take correct data.table if pre-normalized or raw data
      if("dat_gene" %in% names(proteome_data)){
        message("Imputing normalized data...")
        dt_to_impute <- copy(proteome_data$dat_gene)
        setnames(dt_to_impute, "GeneName", "ID_peptide")
      } else{
        message("Imputing raw data...")
        dt_to_impute <- copy(proteome_data$psm_log_prot_df)
      }
      
      samples_pep <- data.table("sample" = proteome_data$c_anno$sample, 
                                "names_sample" =make.names(proteome_data$c_anno$sample))
      numeric_df <- dt_to_impute[, .(ID_peptide, is.na(.SD)), .SDcols = samples_pep$sample]
      numeric_df <- numeric_df[, .(ID_peptide,
                                   impute_level = rowSums(.SD) / length(samples_pep$names_sample) * 100), .SDcols = samples_pep$names_sample]
      numeric_df <- numeric_df[proteome_data$psm_anno_df[, c("ID_peptide","symbol")], on="ID_peptide"]
      numeric_df <- numeric_df[, .(protein_impute_level = mean(impute_level)), by="symbol"]
      psm_prot <- copy(proteome_data$psm_anno_df)
      proteome_data$psm_anno_df <- psm_prot[numeric_df, on="symbol"]
      
      dt_imputed <- impute_matrix(mat = dt_to_impute, c_anno = proteome_data$c_anno, type=type)
      
      if("dat_gene" %in% names(proteome_data)){
        setnames(dt_imputed, "ID_peptide", "GeneName")
        proteome_data$dat_gene <- dt_imputed
      } else{
        proteome_data$psm_log_prot_df <- dt_imputed
      }
    }
    
    ## PEPTIDES
    if("dat_pep" %in% names(proteome_data) | "psm_log_pet_df" %in% names(proteome_data)){
      if(type=="pcaMethods"){
        if("dat_pep" %in% names(proteome_data)){
          # Add missing values
          samples_pep <- data.table("sample" = proteome_data$c_anno$sample, 
                                    "names_sample" =make.names(proteome_data$c_anno$sample))
          numeric_df <- proteome_data$dat_pep[, .(ID_peptide, is.na(.SD)), .SDcols = samples_pep$sample]
          numeric_df <- numeric_df[, .(ID_peptide,
                                       impute_level = rowSums(.SD) / length(samples_pep$names_sample) * 100), .SDcols = samples_pep$names_sample]
          psm_peptide <- copy(proteome_data$psm_peptide_table)
          proteome_data$psm_peptide_table <- psm_peptide[numeric_df, on="ID_peptide"]
          
          proteome_data$dat_pep <- impute_matrix(proteome_data$dat_pep, proteome_data$c_anno, type=type)
        } else{
          stop("'pcaMethods' selected, please provide the normalized data after 'normalization_ProTN'")
        }
      } else{
        # Add missing values
        # Take correct data.table if pre-normalized or raw data
        if("dat_pep" %in% names(proteome_data)){
          message("Imputing normalized data...")
          dt_to_impute <- copy(proteome_data$dat_pep)
        } else{
          message("Imputing raw data...")
          dt_to_impute <- copy(proteome_data$psm_log_pet_df)
        }
        
        samples_pep <- data.table("sample" = proteome_data$c_anno$sample, 
                                  "names_sample" =make.names(proteome_data$c_anno$sample))
        numeric_df <- dt_to_impute[, .(ID_peptide, is.na(.SD)), .SDcols = samples_pep$sample]
        numeric_df <- numeric_df[, .(ID_peptide,
                                     impute_level = rowSums(.SD) / length(samples_pep$names_sample) * 100), .SDcols = samples_pep$names_sample]
        psm_peptide <- copy(proteome_data$psm_peptide_table)
        proteome_data$psm_peptide_table <- psm_peptide[numeric_df, on="ID_peptide"]
        
        
        if("dat_pep" %in% names(proteome_data)){
          dt_imputed <- impute_matrix(mat = dt_to_impute, c_anno = proteome_data$c_anno, type=type)
          
          proteome_data$dat_pep <- dt_imputed
        } else{
          proteome_data$psm_log_pet_df <- copy(proteome_data$psm_log_prot_df)
        }
      }
    }
  } else if(("c_anno_proteome" %in% names(proteome_data)) & ("c_anno_phospho" %in% names(proteome_data))){
    ## PROTEINS
    if(type=="pcaMethods"){
      # Add missing values
      if("dat_gene" %in% names(proteome_data)){
        samples_pep <- data.table("sample" = proteome_data$c_anno_proteome$sample, 
                                  "names_sample" =make.names(proteome_data$c_anno_proteome$sample))
        numeric_df <- proteome_data$dat_gene[, .(GeneName, is.na(.SD)), .SDcols = samples_pep$sample]
        numeric_df <- numeric_df[, .(GeneName,
                                     protein_impute_level = rowSums(.SD) / length(samples_pep$names_sample) * 100), .SDcols = samples_pep$names_sample]
        psm_prot <- copy(proteome_data$psm_anno_df)
        proteome_data$psm_anno_df <- psm_prot[numeric_df, on = .(symbol = GeneName)]
        
        
        proteome_data$dat_gene <- impute_matrix(proteome_data$dat_gene, proteome_data$c_anno_proteome, type="pcaMethods_protein")
      } else{
        stop("'pcaMethods' selected, please provide the normalized data after 'normalization_ProTN'")
      }
    } else{
      # Add missing values
      # Take correct data.table if pre-normalized or raw data
      if("dat_gene" %in% names(proteome_data)){
        message("Imputing normalized data...")
        dt_to_impute <- copy(proteome_data$dat_gene)
        setnames(dt_to_impute, "GeneName", "ID_peptide")
      } else{
        message("Imputing raw data...")
        dt_to_impute <- copy(proteome_data$psm_log_prot_df)
      }
      
      samples_pep <- data.table("sample" = proteome_data$c_anno_proteome$sample, 
                                "names_sample" =make.names(proteome_data$c_anno_proteome$sample))
      numeric_df <- dt_to_impute[, .(ID_peptide, is.na(.SD)), .SDcols = samples_pep$sample]
      numeric_df <- numeric_df[, .(ID_peptide,
                                   impute_level = rowSums(.SD) / length(samples_pep$names_sample) * 100), .SDcols = samples_pep$names_sample]
      numeric_df <- numeric_df[proteome_data$psm_anno_df[, c("ID_peptide","symbol")], on="ID_peptide"]
      numeric_df <- numeric_df[, .(protein_impute_level = mean(impute_level)), by="symbol"]
      psm_prot <- copy(proteome_data$psm_anno_df)
      proteome_data$psm_anno_df <- psm_prot[numeric_df, on="symbol"]
      
      dt_imputed <- impute_matrix(mat = dt_to_impute, c_anno = proteome_data$c_anno_proteome, type=type)
      
      if("dat_gene" %in% names(proteome_data)){
        setnames(dt_imputed, "ID_peptide", "GeneName")
        proteome_data$dat_gene <- dt_imputed
      } else{
        proteome_data$psm_log_prot_df <- dt_imputed
      }
    }
    
    ## PEPTIDES
    if("dat_pep" %in% names(proteome_data) | "psm_log_pet_df" %in% names(proteome_data)){
      if(type=="pcaMethods"){
        # Add missing values
        if("dat_pep" %in% names(proteome_data)){
          samples_pep <- data.table("sample" = proteome_data$c_anno_phospho$sample, 
                                    "names_sample" =make.names(proteome_data$c_anno_phospho$sample))
          numeric_df <- proteome_data$dat_pep[, .(ID_peptide, is.na(.SD)), .SDcols = samples_pep$sample]
          numeric_df <- numeric_df[, .(ID_peptide,
                                       impute_level = rowSums(.SD) / length(samples_pep$names_sample) * 100), .SDcols = samples_pep$names_sample]
          psm_peptide <- copy(proteome_data$psm_peptide_table)
          proteome_data$psm_peptide_table <- psm_peptide[numeric_df, on="ID_peptide"]
          
          proteome_data$dat_pep <- impute_matrix(proteome_data$dat_pep, proteome_data$c_anno_phospho, type=type)
        } else{
          stop("'pcaMethods' selected, please provide the normalized data after 'normalization_ProTN'")
        }
      } else{
        # Add missing values
        # Take correct data.table if pre-normalized or raw data
        if("dat_pep" %in% names(proteome_data)){
          message("Imputing normalized data...")
          dt_to_impute <- copy(proteome_data$dat_pep)
        } else{
          message("Imputing raw data...")
          dt_to_impute <- copy(proteome_data$psm_log_pet_df)
        }
        
        samples_pep <- data.table("sample" = proteome_data$c_anno_phospho$sample, 
                                  "names_sample" =make.names(proteome_data$c_anno_phospho$sample))
        numeric_df <- dt_to_impute[, .(ID_peptide, is.na(.SD)), .SDcols = samples_pep$sample]
        numeric_df <- numeric_df[, .(ID_peptide,
                                     impute_level = rowSums(.SD) / length(samples_pep$names_sample) * 100), .SDcols = samples_pep$names_sample]
        psm_peptide <- copy(proteome_data$psm_peptide_table)
        proteome_data$psm_peptide_table <- psm_peptide[numeric_df, on="ID_peptide"]
        
        dt_imputed <- impute_matrix(mat = dt_to_impute, c_anno = proteome_data$c_anno_phospho, type=type)
        
        if("dat_pep" %in% names(proteome_data)){
          proteome_data$dat_pep <- dt_imputed
        } else{
          proteome_data$psm_log_pet_df <- dt_imputed
        }
      }
    }
    
  } else{
    stop("Missing sample annotation!")
  }
  return(proteome_data)
}

#' Matrix Imputation Function
#'
#' This function imputes missing values in a matrix using either PhosR's imputation method (scImpute + tImpute) for data or a median imputation approach if PhosR fails.
#'
#' @param mat Data frame or matrix; the intensity matrix to be imputed. The first column must be a peptide ID column, and the rest should contain intensity values.
#' @param c_anno Data frame; a data frame containing sample annotations, with columns for the sample ID and condition.
#' @param type character; Type of imputation selected (`phosr`, `gaussian`, `missForest`, `pcaMethods`).
#'
#' @return A data.table with the imputed intensity matrix, retaining the original peptide IDs.
#'
#' @import data.table
#' @import PhosR
#' @import missForest
#' @import doParallel
#' @import parallel
#' @import pcaMethods
#' @export
impute_matrix <- function(mat, c_anno, type) {
  mat<-tryCatch({
    if(type == "phosr"){
      mat_df <- as.data.frame(mat[,-1], row.names = as.character(mat$ID_peptide))
      ppe <- suppressWarnings(PhosphoExperiment(assays = list(Quantification = as.matrix(mat_df))))
      grps <- as.factor(c_anno[sample %in% colnames(mat_df), condition])
      set.seed(42)
      ppe <- scImpute(ppe, 0.6, grps)
      set.seed(42)
      ppe <- tImpute(ppe, assay = "imputed")
      mat_df[is.na(mat_df)] <- ppe@assays@data@listData$imputed[is.na(mat_df)]
      message("Performed PhosR imputation")
      as.data.table(mat_df, keep.rownames = "ID_peptide")
    } else if(type == "gaussian"){
      mat_df <- as.data.frame(mat[,-1], row.names = as.character(mat$ID_peptide))
      vec_all <- unlist(mat_df)
      vec_numeric <- vec_all[!is.na(vec_all)]
      orig_stats <- c(mean(vec_numeric), sd(vec_numeric))
      imp_stats <- c(orig_stats[1] -(orig_stats[2]*2),orig_stats[2]*0.3)
      
      set.seed(42)
      imp_matrix <- matrix(rnorm(length(vec_all), mean = imp_stats[1], sd = imp_stats[2]), ncol = ncol(mat_df))
      mat_df[is.na(mat_df)] <- imp_matrix[is.na(mat_df)]
      message("Performed Gaussian imputation")
      as.data.table(mat_df, keep.rownames = "ID_peptide")
    } else if(type == "missForest"){
      mat_df <- as.data.frame(mat[,-1], row.names = as.character(mat$ID_peptide))
      set.seed(42)
      n_cores <- detectCores() - 2
      cl <- makeCluster(min(n_cores, ncol(mat_df)))
      registerDoParallel(cl)
      imp_mat_dt <- missForest(mat_df, verbose = T, parallelize = "forest", 
                               maxiter = 4, ntree = 100)
      stopCluster(cl)
      mat_df <- imp_mat_dt$ximp
      message("Performed missForest imputation")
      as.data.table(mat_df, keep.rownames = "ID_peptide")
    } else if(type == "pcaMethods"){
      mat <- mat[, lapply(.SD, function(x) replace(x, is.nan(x), NA))]
      mat_df <- as.data.frame(mat[,-1], row.names = as.character(mat$ID_peptide))
      y_pcaMethods_impute <- pca(mat_df, method="svdImpute", nPcs=4, center = TRUE, verbose=TRUE)
      mat_df <- completeObs(y_pcaMethods_impute)
      message("Performed pcaMethods imputation")
      as.data.table(mat_df, keep.rownames = "ID_peptide")
    } else if(type == "pcaMethods_protein"){
      mat <- mat[, lapply(.SD, function(x) replace(x, is.nan(x), NA))]
      mat_df <- as.data.frame(mat[,-1], row.names = as.character(mat$GeneName))
      y_pcaMethods_impute <- pca(mat_df, method="svdImpute", nPcs=4, center = TRUE, verbose=TRUE)
      mat_df <- completeObs(y_pcaMethods_impute)
      as.data.table(mat_df, keep.rownames = "GeneName")
    }
  }, error = function(cond) {
    stop(paste0("Error Imputation method: Try another one.", as.character(cond)))
  })
  return(mat)
}
