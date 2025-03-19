#' Normalization
#'
#' @param proteome_data Object proTN
#' @return A list of data tables.
#' @details \strong{ProTN}
#' @examples 
#' ## ## Example:
#' ## example
#' ## example2
#' @import data.table
#' @export
normalization_ProTN <- function(proteome_data) {

  # Summarize into proteins and normalize by median
  dat_psm <- merge.data.table(proteome_data$psm_anno_df[, .(ID_peptide, symbol)], proteome_data$psm_log_prot_df, by = "ID_peptide")
  proteome_data$dat_gene <- as.data.table(medianSweeping(dat_psm, group_col = 2), keep.rownames = "GeneName")
  
  # Median centering log2 table for peptides
  dat_pep <- copy(proteome_data$psm_log_pet_df)
  dat_pep <- as.data.frame(dat_pep[,-1], row.names = as.character(dat_pep$ID_peptide))
  proteome_data$dat_pep <- as.data.table(equalMedianNormalization(dat_pep), keep.rownames = "ID_peptide")
  
  message("Normalization DONE.")
  return(proteome_data)
}
