#' Intensity Imputation for Proteomics Data
#'
#' This function performs imputation of missing intensity values in a proteomics data matrix. 
#' It uses PhosR's imputation methods or falls back to median imputation if PhosR's method fails. 
#' The function works on intensity matrices for both proteome and phosphoproteome data.
#'
#' @param proteome_data List; a list containing proteomics data. It must include either `c_anno` (sample annotations for proteome) or `c_anno_proteome` and `c_anno_phospho` (sample annotations for proteome and phosphoproteome, respectively). It also requires `psm_log_prot_df` and `psm_log_pet_df` (log-transformed intensity matrices for proteome and phosphoproteome).
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
#' @export
impute_intensity <- function(proteome_data) {
  message("Imputation of intensity matrix in progress...")
  # Load PhosR functions
  suppressWarnings(lapply(list.files("R/PhosR", pattern = ".R", full.names = T), function(x){if(!file.info(x)$isdir){source(x)}}))
  
  if(("c_anno" %in% names(proteome_data))){
    proteome_data$psm_log_prot_df <- impute_matrix(mat = proteome_data$psm_log_prot_df, c_anno = proteome_data$c_anno)
    
    proteome_data$psm_log_pet_df <- impute_matrix(proteome_data$psm_log_pet_df, proteome_data$c_anno)
  } else if(("c_anno_proteome" %in% names(proteome_data)) & ("c_anno_phospho" %in% names(proteome_data))){
    proteome_data$psm_log_prot_df <- impute_matrix(mat = proteome_data$psm_log_prot_df, c_anno = proteome_data$c_anno_proteome)
    
    proteome_data$psm_log_pet_df <- impute_matrix(proteome_data$psm_log_pet_df, proteome_data$c_anno_phospho)
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
#'
#' @return A data.table with the imputed intensity matrix, retaining the original peptide IDs.
#'
#' @import data.table
#' @import PhosR
#' @export
impute_matrix <- function(mat, c_anno) {
  mat_df <- as.data.frame(mat[,-1], row.names = as.character(mat$ID_peptide))
  mat<-tryCatch({
    ppe <- suppressWarnings(PhosphoExperiment(assays = list(Quantification = as.matrix(mat_df))))
    grps <- as.factor(c_anno[sample %in% colnames(mat_df), condition])
    set.seed(42)
    ppe <- scImpute(ppe, 0.6, grps)
    set.seed(42)
    ppe <- tImpute(ppe, assay = "imputed")
    mat_df[is.na(mat_df)] <- ppe@assays@data@listData$imputed[is.na(mat_df)]
    message("Performed PhosR imputation")
    as.data.table(mat_df, keep.rownames = "ID_peptide")
  }, error = function(cond) {
    message(paste0("Error PhosR: ",as.character(cond)))
    vec_all <- unlist(mat_df)
    vec_numeric <- vec_all[!is.na(vec_all)]
    orig_stats <- c(mean(vec_numeric), sd(vec_numeric))
    imp_stats <- c(orig_stats[1] - (orig_stats[2] * 1.8), orig_stats[2] * 0.3)
    
    set.seed(42)
    imp_matrix <- matrix(rnorm(length(vec_all), mean = imp_stats[1], sd = imp_stats[2]), ncol = ncol(mat_df))
    mat_df[is.na(mat_df)] <- imp_matrix[is.na(mat_df)]
    message("Performed Median imputation")
    as.data.table(mat_df, keep.rownames = "ID_peptide")
  })
  return(mat)
}
