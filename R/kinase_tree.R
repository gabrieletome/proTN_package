kinase_activity_calculation <- function(dirOutput_kinase, 
                                        formule_CORAL, 
                                        comparison, 
                                        dat_pep, 
                                        deps_pep_l_df, 
                                        psm_peptide_table, 
                                        c_anno,
                                        phosR_thr = 0.7,
                                        species = "Homo sapiens") {
  data("KinaseMotifs", envir = environment())
  data("PhosphoSitePlus", envir = environment())
  
  deps_filtered <- deps_pep_l_df[comparison == deps_pep_l_df$comp & deps_pep_l_df$class %in% c("+", "-"), .(id)]
  selected_peptides <- psm_peptide_table[ID_peptide %in% deps_filtered$id,]
  tmp_dat_pep <- dat_pep[ID_peptide %in% selected_peptides$ID_peptide, ]
  tmp_dat_pep <- unique(tmp_dat_pep[selected_peptides[, c("ID_peptide","GeneName","Annotated_Sequence")], on="ID_peptide"][, ID_peptide := NULL])
  dat_pep_matrix <- data.frame(tmp_dat_pep[,-c("GeneName","Annotated_Sequence")], row.names = make.names(tmp_dat_pep$GeneName, unique=T))
  
  ppe <- PhosphoExperiment(assays = list(Quantification = as.matrix(dat_pep_matrix)))
  
  GeneSymbol(ppe) <- tmp_dat_pep$GeneName
  Sequence(ppe) <- tmp_dat_pep$Annotated_Sequence
  
  mat <- assay(ppe, "Quantification")
  
  if(species == "Mus musculus"){
    substrate.list <- sapply(PhosphoSite.mouse, function(y) unlist(lapply(y, function(x) gsub(";.+", "", x))))
  } else if(species == "Homo sapiens"){
    substrate.list <- sapply(PhosphoSite.human, function(y) unlist(lapply(y, function(x) gsub(";.+", "", x))))
  } else{
    warning("Species not recognised. Select between: 'Homo sapiens' and 'Mus musculus'. Select default 'Homo sapiens'")
    substrate.list <- sapply(PhosphoSite.human, function(y) unlist(lapply(y, function(x) gsub(";.+", "", x))))
  }
  
  mat.std <- standardise(mat)
  rownames(mat.std) <- toupper(sub("\\..*", "", rownames(mat.std)))
  kssMat <- kinaseSubstrateScore_local(substrate.list, comp = comparison, mat = mat.std, seqs=na.omit(ppe@Sequence), 5, 1, "human", TRUE, file.path(dirOutput_kinase, paste0(comparison, "_kinase_heatmap.pdf")))
  set.seed(42)
  
  tryCatch({
    
    predMat <- kinaseSubstratePred(kssMat, inclusion = 5, cs = phosR_thr)
    
    design <- model.matrix(~0 + c_anno$condition)
    colnames(design) <- unique(c_anno$condition)
    rownames(design) <- c_anno$sample
    contrast <- makeContrasts(contrasts = formule_CORAL[comparison], levels = design)
    
    # Differential analysis for the svg
    mean_kinase_activity <- lapply(rownames(contrast)[which((contrast != 0)[, 1])], 
                                   function(x) {
                                     rowMeans(kssMat$ksActivityMatrix[, make.names(c_anno[x == condition, sample])])[colnames(predMat)] * contrast[x, ]
                                   }
    )
    kinase_Act <- Reduce("+", mean_kinase_activity)
    
    dt <- as.data.table(kssMat$ksActivityMatrix, keep.rownames = "GeneName")
    message("Saving differential kinase activity.")
    fwrite(as.data.table(kinase_Act, keep.rownames = "GeneName"), file = file.path(dirOutput_kinase, paste0(comparison, "_kinase_activity_differential.txt")), col.names = FALSE, quote = FALSE)
    message("Saving kinase activity matrix.")
    write_xlsx(dt, file.path(dirOutput_kinase, paste0(comparison, "_kinase_activity_matrix.xlsx")), col_names = T)
    
    if(species == "Homo sapiens"){
      message("Preparing svg tree...")
      renderSvg(comparison, kinase_Act, dirOutput_kinase)
    }
    
    return(list(dt))
  }, error = function(err){
    print("Error: probably kinase are not passing the filter. We suggest to decrease the threshold.")
    return(list("Error: probably kinase are not passing the filter. We suggest to decrease the threshold."))
  })
}


#' Kinase Activity Tree Calculation
#'
#' This function calculates kinase activity for differential analysis results, creating a tree-based visualization of kinase activity. It uses CORAL-based formulae and generates kinase activity plots. The function allows the user to specify whether to include phospho-control data or not.
#'
#' @param proteome_data List; a list containing proteomics data.
#' @param differential_results List; a list containing differential analysis results, specifically `peptide_results_long` for peptide-wise differential analysis.
#' @param formule_CORAL List; a list of formulae for CORAL-based kinase activity calculations. It must include formulae for the comparisons of interest.
#' @param species Character. The species name "Homo sapiens" or "Mus musculus" (default: `"Homo sapiens"`). If selected "Mus musculus" the visial CORAL kinome tree will not be done.
#' @param dirOutput Character; the directory where the results will be stored. Default is "results_ProTN".
#' @param subfold Character; the subdirectory within `dirOutput` to store the results, default is "pics".
#' @param phospho_ctrl Logical; whether to include phospho-control data in the calculations. Default is FALSE (excludes phospho-control data).
#' @param phosR_thr Integer; Threshold for kinase activity (0-1). Default is 0.7
#'
#' @return A list of SVG plots, one for each comparison in the `formule_CORAL` list, showing the kinase activity tree for each comparison.
#'
#' @examples
#' \dontrun{
#' kinase_tree_results <- kinase_tree(proteome_data = proteome_data, 
#'                                    differential_results = differential_results, 
#'                                    formule_CORAL = formule_CORAL)
#' }
#'
#' @import data.table
#' @import ggplot2
#' @import writexl
#' @importFrom SummarizedExperiment assay
#' @importFrom readr read_tsv
#' @import rsvg
#' @import RColorBrewer
#' @import colourpicker
#' @import rjson
#' @import data.tree
#' @import pheatmap
#' 
#' @export
kinase_tree <- function(proteome_data, differential_results, formule_CORAL, species="Homo sapiens", 
                        dirOutput = "results_ProTN", subfold="pics", phospho_ctrl = FALSE, phosR_thr = 0.7) {

  dir.create(file.path( dirOutput, subfold, "kinaseTree"), showWarnings = FALSE, recursive = T)
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
                                                          comparison = comparison, 
                                                          dat_pep = dat_pep, 
                                                          deps_pep_l_df = deps_pep_l_df, 
                                                          psm_peptide_table = psm_peptide_table, 
                                                          c_anno = c_anno, 
                                                          phosR_thr = phosR_thr,
                                                          species = species)
  }
  
  return(list_svg)
}


