#' Function to perform batch correction
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
batch_correction <- function(proteome_data, batch_col="batch") {
  message("A batch effect correction is applied to the dataset using the proBatch package.")
  
  batch_annotation <- copy(proteome_data$c_anno)
  batch_annotation[, order := .I]
  
  # Batch correction of peptides
  if(("c_anno" %in% names(proteome_data))){
    batch_annotation <- copy(proteome_data$c_anno)
    batch_annotation[, order := .I]
  } else if(("c_anno_proteome" %in% names(proteome_data)) & ("c_anno_phospho" %in% names(proteome_data))){
    batch_annotation <- copy(proteome_data$c_anno_phospho)
    batch_annotation[, order := .I]
  } else{
    stop("Missing sample annotation!")
  }
  dat_pep <- copy(proteome_data$dat_pep)
  dat_pep <- as.data.frame(dat_pep[,-1], row.names = as.character(dat_pep$ID_peptide))
  dat_pep_long <- matrix_to_long(dat_pep)
  setnames(dat_pep_long, "FullRunName", "sample")
  comBat_df_pep <- correct_with_ComBat_df(dat_pep_long, batch_annotation, batch_col = batch_col, sample_id_col = "sample")
  proteome_data$dat_pep <- as.data.table(as.data.frame(long_to_matrix(comBat_df_pep, sample_id_col = "sample")), keep.rownames = "ID_peptide")
  
  # Batch correction of genes
  if(("c_anno" %in% names(proteome_data))){
    batch_annotation <- copy(proteome_data$c_anno)
    batch_annotation[, order := .I]
  } else if(("c_anno_proteome" %in% names(proteome_data)) & ("c_anno_phospho" %in% names(proteome_data))){
    batch_annotation <- copy(proteome_data$c_anno_proteome)
    batch_annotation[, order := .I]
  } else{
    stop("Missing sample annotation!")
  }
  dat_gene <- copy(proteome_data$dat_gene)
  dat_gene <- as.data.frame(dat_gene[,-1], row.names = as.character(dat_gene$ID_peptide))
  dat_gene_long <- matrix_to_long(dat_gene)
  setnames(dat_gene_long, "FullRunName", "sample")
  comBat_df_gene <- correct_with_ComBat_df(dat_gene_long, batch_annotation, batch_col = batch_col, sample_id_col = "sample")
  proteome_data$dat_gene <- as.data.table(as.data.frame(long_to_matrix(comBat_df_gene, sample_id_col = "sample")), keep.rownames = "GeneName")
  
  return(proteome_data)
}
