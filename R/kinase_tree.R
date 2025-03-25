kinase_activity_calculation <- function(dirOutput_kinase, formule_CORAL, comp, dat_pep, deps_pep_l_df, psm_peptide_table, c_anno) {
  data("KinaseMotifs")
  data("PhosphoSitePlus")
  
  # lapply(list.files("R/PhosR", full.names = TRUE), source)
  
  deps_filtered <- deps_pep_l_df[comp == deps_pep_l_df$comp & deps_pep_l_df$class %in% c("+", "-"), .(id)]
  selected_peptides <- psm_peptide_table[ID_peptide %in% deps_filtered$id,]
  tmp_dat_pep <- dat_pep[ID_peptide %in% selected_peptides$ID_peptide, ]
  tmp_dat_pep <- unique(tmp_dat_pep[selected_peptides[, c("ID_peptide","GeneName","Annotated_Sequence")], on="ID_peptide"][, ID_peptide := NULL])
  dat_pep_matrix <- data.frame(tmp_dat_pep[,-c("GeneName","Annotated_Sequence")], row.names = make.names(tmp_dat_pep$GeneName, unique=T))
  
  ppe <- PhosphoExperiment(assays = list(Quantification = as.matrix(dat_pep_matrix)))
  
  GeneSymbol(ppe) <- tmp_dat_pep$GeneName
  Sequence(ppe) <- tmp_dat_pep$Annotated_Sequence
  
  mat <- SummarizedExperiment::assay(ppe, "Quantification")
  substrate.list <- sapply(PhosphoSite.human, function(y) unlist(lapply(y, function(x) gsub(";.+", "", x))))
  
  mat.std <- standardise(mat)
  rownames(mat.std) <- toupper(sub("\\..*", "", rownames(mat.std)))
  
  kssMat <- kinaseSubstrateScore(substrate.list, mat.std, na.omit(ppe@Sequence), 5, 1, "human", TRUE)
  set.seed(42)
  predMat <- kinaseSubstratePred(kssMat, inclusion = 5)
  
  design <- model.matrix(~0 + c_anno$condition)
  colnames(design) <- unique(c_anno$condition)
  rownames(design) <- c_anno$sample
  contrast <- makeContrasts(contrasts = formule_CORAL[comp], levels = design)
  
  # Differential analysis for the svg
  mean_kinase_activity <- lapply(rownames(contrast)[which((contrast != 0)[, 1])], 
                                 function(x) {
                                   rowMeans(kssMat$ksActivityMatrix[, grepl(x, colnames(kssMat$ksActivityMatrix))])[colnames(predMat)] * contrast[x, ]
                                 }
  )
  kinase_Act <- Reduce("+", mean_kinase_activity)
  
  dt <- as.data.table(kssMat$ksActivityMatrix, keep.rownames = "GeneName")
  fwrite(data.table(kinase_Act, keep.rownames = T), file = file.path(dirOutput_kinase, paste0(comp, "_kinase_activity_differential.txt")), col.names = FALSE, quote = FALSE)
  writexl::write_xlsx(dt, file.path(dirOutput_kinase, paste0(comp, "_kinase_activity_matrix.xlsx")), col_names = T)
  
  renderSvgPanZoom(comp, kinase_Act, dirOutput_kinase)
  
  return(list(dt))
}


#' Kinase Tree
#'
#' @param software Must be \strong{PD} for Proteome Discoverer or \strong{MQ} for MaxQuant
#' @param folder Path with files
#' @return A list of data tables.
#' @details \strong{ProTN}
#' @examples 
#' ## ## Example:
#' ## example
#' ## example2
#' @import data.table
#' @export
kinase_tree <- function(proteome_data, differential_results, formule_CORAL, 
                        dirOutput = "results_ProTN", subfold="pics", phospho_ctrl = FALSE) {
  # oldwd <- getwd()
  source("R/functions.R")
  lapply(list.files("R/PhosR", full.names = TRUE), source)
  source("R/CORAL/global.R")
  source("R/CORAL/server.R")
  
  dir.create(file.path( dirOutput, subfold, "kinaseTree"), showWarnings = FALSE)
  dirOutput_kinase <- file.path(dirOutput, subfold, "kinaseTree")
  
  dat_pep = copy(proteome_data$dat_pep)
  deps_pep_l_df <- copy(differential_results$peptide_results_long)
  psm_peptide_table <- copy(proteome_data$psm_peptide_table)
  
  if(("c_anno" %in% names(proteome_data))){
    c_anno <- copy(proteome_data$c_anno)
  } else if(("c_anno_phospho" %in% names(proteome_data))){
    c_anno <- copy(proteome_data$c_anno_phospho)
  } else{
    stop("Missing sample annotation!")
  }
  
  if(!phospho_ctrl){
    compToKeep <- unique(grep("_Phospho_CTRL", names(formule_CORAL), value = T, invert = T))
    formule_CORAL <- formule_CORAL[compToKeep]
  }
  
  list_svg <- list()
  for(comparison in names(formule_CORAL)){
    message(paste0("\nDoing: ",comparison,"\n"))
    list_svg[[comparison]] <- kinase_activity_calculation(dirOutput_kinase = dirOutput_kinase, 
                                                          formule_CORAL = formule_CORAL, 
                                                          comp = comparison, 
                                                          dat_pep = dat_pep, 
                                                          deps_pep_l_df = deps_pep_l_df, 
                                                          psm_peptide_table = psm_peptide_table, 
                                                          c_anno = c_anno)
  }
  
  # setwd(oldwd)
  return(list_svg)
}


