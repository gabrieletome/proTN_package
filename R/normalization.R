#' Normalize Data
#'
#' This function normalizes data by summarizing peptide intensities into protein-level data, applying median-based normalization. 
#' It also performs median centering for peptide-level log2 data.
#'
#' @param proteome_data A list containing proteome data, including peptide annotations, protein-level data, and peptide-level data.
#'
#' @return A list containing the normalized proteome data:
#'   \item{dat_gene}{Data table containing the protein-level normalized data.}
#'   \item{dat_pep}{Data table containing the peptide-level normalized data.}
#'
#' @examples
#' \dontrun{
#' proteome_data <- normalization_ProTN(proteome_data)
#' }
#'
#' @import data.table
#' @import DEqMS
#' @export
normalization_ProTN <- function(proteome_data) {

  # Summarize into proteins and normalize by median
  dat_psm <- merge.data.table(proteome_data$psm_anno_df[, .(ID_peptide, symbol)], proteome_data$psm_log_prot_df, by = "ID_peptide")
  proteome_data$dat_gene <- as.data.table(medianSweeping(dat_psm, group_col = 2), keep.rownames = "GeneName")
  
  # Median centering log2 table for peptides
  if("psm_log_pet_df" %in% names(proteome_data)){
    dat_pep <- copy(proteome_data$psm_log_pet_df)
    dat_pep <- as.data.frame(dat_pep[,-1], row.names = as.character(dat_pep$ID_peptide))
    proteome_data$dat_pep <- as.data.table(equalMedianNormalization(dat_pep), keep.rownames = "ID_peptide")
  }
  
  message("Normalization DONE.")
  return(proteome_data)
}
