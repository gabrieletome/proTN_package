#' Function to generate barplots for differential proteins and peptides
#'
#' @param differential_results Object proTN
#' @details \strong{ProTN}
#' @examples 
#' ## ## Example:
#' ## example
#' ## example2
#' @import data.table
#' @export
generate_differential_barplots <- function(differential_results, data_type="protein", 
                                           color_contrast=NULL, phospho_ctrl = FALSE) {
  if(("protein_results_long" %in% names(differential_results))){
    phospho_with_proteome = FALSE
  } else if(("peptide_results_long" %in% names(differential_results))){
    phospho_with_proteome = TRUE
    data_type="peptide"
  } else{
    stop("Error in differential results paramenter! Verify the presence of protein_results_long or peptide_results_long")
  }
  
  deps_l_df <- if (data_type == "protein"){ copy(differential_results$protein_results_long)
  } else if(data_type == "peptide") { copy(differential_results$peptide_results_long)
      } else{stop("data_type MUST be \"protein\" or \"peptide\"")}
  
  if(!phospho_ctrl){
    compToKeep <- unique(grep("_Phospho_CTRL", deps_l_df$comp, value = T, invert = T))
    deps_l_df <- deps_l_df[comp %in% compToKeep]
  }
  
  if (nrow(deps_l_df[class != "="]) > 0) {
    lolli_df <- data.table(
      "comp" = deps_l_df$comp[deps_l_df$class %in% c("+", "-")],
      "class" = factor(deps_l_df$class[deps_l_df$class %in% c("+", "-")])
    )
    lolli_df <- lolli_df[, .(N = .N), by = c("comp", "class")]
    lolli_df$id <- paste(lolli_df$comp, lolli_df$class, sep = "_")
    
    if(is.null(color_contrast)){
      color_contrast=c("#664069","#8A628D")
      col_vec <- rep(as.vector(t(color_contrast)), uniqueN(lolli_df$comp))
      message("Set default colors.")
    } else{
      col_vec <- as.vector(t(color_contrast))
    }
    
    tryCatch({
      names(col_vec)[seq(1, length(lolli_df$comp), by = 2)] <- paste0(lolli_df$comp[seq(1, length(lolli_df$comp), by = 2)], "_+")
      names(col_vec)[seq(2, length(lolli_df$comp), by = 2)] <- paste0(lolli_df$comp[seq(1, length(lolli_df$comp), by = 2)], "_-")
    }, error = function(cond){
      stop("Color must be a vector of the same length of the number of comparison.")
    })
    
    
    pDEPs <- deps_b2b_lollipop(
      input_df = lolli_df,
      break_vec = seq(0, plyr::round_any(max(lolli_df$N) * 1.2, if (max(lolli_df$N) > 90) 100 else 10, f = ceiling),
                      round(plyr::round_any(max(lolli_df$N) * 1.2, if (max(lolli_df$N) > 90) 100 else 10, f = ceiling) / 5, 0)),
      fill_vec = col_vec,
      color_vec = col_vec,
      char_max = 30,
      shape_vec = c(16, 21),
      position_dodge = 0
    )
    
    message(sprintf("%s \n", toupper(data_type)))
    return(list("plot"=pDEPs))
  } else {
    stop("No gene differentiated discovered.")
  }
}


#' Function to generate volcano plots for differential proteins and peptides
#'
#' @param differential_results Object proTN
#' @details \strong{ProTN}
#' @examples 
#' ## ## Example:
#' ## example
#' ## example2
#' @import data.table
#' @export
generate_volcano_plots <- function(differential_results, data_type=NULL, 
                                   comparison=NULL, fc_thr=0.75, pval_fdr = "p_val", 
                                   pval_thr=0.05, color_contrast=NULL) {
  if(is.null(comparison)){
    stop("Provide a valid comparison")
  }
  if(("protein_results_long" %in% names(differential_results))){
    phospho_with_proteome = FALSE
  } else if(("peptide_results_long" %in% names(differential_results))){
    phospho_with_proteome = TRUE
    data_type="peptide"
  } else{
    stop("Error in differential results paramenter! Verify the presence of protein_results_long or peptide_results_long")
  }
  
  deps_l_df <- if (data_type == "protein"){ copy(differential_results$protein_results_long)
  } else if(data_type == "peptide") { copy(differential_results$peptide_results_long)
  } else{stop("data_type MUST be \"protein\" or \"peptide\"")}
  
  if(!(comparison %in% unique(deps_l_df$comp))){
    stop("Provide a valid comparison. Comparison not present in the differential analysis")
  }
    
  if (nrow(deps_l_df[class != "="]) > 0) {
    plotlist <- list()
    
    col <- if (pval_fdr == "p_adj") "p_adj" else if(pval_fdr == "p_val") "p_val" else stop("pval_fdr MUST be p_adj or p_val")
    input_df <- deps_l_df[comp == comparison]
    input_df[, log2_FC := round(log2_FC, 2)]
    input_df[, (col) := round(-log10(get(col)), 2)]
    
    if(is.null(color_contrast)){
      color_contrast=c("#664069","#8A628D")
      col_vec <- rep(as.vector(t(color_contrast)), uniqueN(input_df$comp))
      message("Set default colors.")
    } else{
      col_vec <- as.vector(t(color_contrast))
    }
    
    if(!(length(col_vec < 2))){
      stop("Color must be a vector of the same length == 2")
    }
    
    bs=20
    cmd <- ggplot(input_df, aes(x = log2_FC, y = get(col), col = class, text = paste('</br>Gene: ', id, '</br>Class: ', class, '</br>Log2_FC: ', log2_FC, '</br>', col, ': ', get(col)))) +
      geom_point(pch = 20, cex = 2) +
      geom_hline(yintercept = -log10(pval_thr), col = "black") +
      geom_vline(xintercept = c(-fc_thr, fc_thr), col = "black") +
      ggtitle(paste0("Volcano Plot of ", comparison)) +
      ylab(if (pval_fdr=="p_adj") "-log10(fdr)" else "-log10(p_val)") +
      scale_color_manual(values = c("+" = col_vec[1], "-" = col_vec[2], "=" = "grey70")) +
      scale_x_continuous(limits = c(min(-max(abs(input_df$log2_FC)), -3), max(max(abs(input_df$log2_FC)), 3))) +
      scale_y_continuous(limits = c(0, max(max(abs(input_df[,get(col)])), 4))) +
      theme_bw(base_size = bs)
    
    plotlist[[comparison]] <- ggplotly(cmd, tooltip = c("text"))
    return(plotlist)
  } else{
    stop("No gene differentiated discovered.")
  }
}



#' Function to perform MDS analysis
#'
#' @param proteome_data Object proTN
#' @param type \strong{protein} or \strong{peptide}
#' @return A list of data tables.
#' @details \strong{ProTN}
#' @examples 
#' ## ## Example:
#' ## example
#' ## example2
#' @import data.table
#' @export
mds_differential_analysis_plot <- function(differential_analysis, proteome_data, type){
  if(("protein_results_long" %in% names(differential_analysis))){
    phospho_with_proteome = FALSE
  } else if(!("protein_results_long" %in% names(differential_analysis)) & 
            ("peptide_results_long" %in% names(differential_analysis))){
    phospho_with_proteome = TRUE
    type="peptide"
  } else{
    stop("Error in differential results paramenter! Verify the presence of protein_results_long or peptide_results_long")
  }
  
  deps_l_df <- if (type == "protein"){ 
    copy(differential_analysis$protein_results_long)
  } else if(type == "peptide") { 
    copy(differential_analysis$peptide_results_long)
  } else{stop("type MUST be \"protein\" or \"peptide\"")}
  
  if(phospho_with_proteome){
    compToKeep <- unique(grep("_Phospho_CTRL", deps_l_df$comp, value = T, invert = T))
    deps_l_df <- deps_l_df[comp %in% compToKeep]
  }
  
  proteome_data_copy <- copy(proteome_data)
  deps_vec_all<- unique(deps_l_df[class != "=", id])
  
  if (type == "protein"){ 
    proteome_data_copy$dat_gene<-proteome_data_copy$dat_gene[GeneName %in% deps_vec_all,]
  } else if(type == "peptide") { 
    proteome_data_copy$dat_pep<-proteome_data_copy$dat_pep[ID_peptide %in% deps_vec_all,]
  } else{stop("type MUST be \"protein\" or \"peptide\"")}
  
  return(mds_plot(proteome_data = proteome_data_copy, type = type))

}


#' Function to perform pca analysis
#'
#' @param proteome_data Object proTN
#' @param type \strong{protein} or \strong{peptide}
#' @return A list of data tables.
#' @details \strong{ProTN}
#' @examples 
#' ## ## Example:
#' ## example
#' ## example2
#' @import data.table
#' @export
pca_differential_analysis_plot <- function(differential_analysis, proteome_data, type){
  if(("protein_results_long" %in% names(differential_analysis))){
    phospho_with_proteome = FALSE
  } else if(!("protein_results_long" %in% names(differential_analysis)) & 
            ("peptide_results_long" %in% names(differential_analysis))){
    phospho_with_proteome = TRUE
    type="peptide"
  } else{
    stop("Error in differential results paramenter! Verify the presence of protein_results_long or peptide_results_long")
  }
  
  deps_l_df <- if (type == "protein"){ 
    copy(differential_results$protein_results_long)
  } else if(type == "peptide") { 
    copy(differential_results$peptide_results_long)
  } else{stop("type MUST be \"protein\" or \"peptide\"")}
  
  if(phospho_with_proteome){
    compToKeep <- unique(grep("_Phospho_CTRL", deps_l_df$comp, value = T, invert = T))
    deps_l_df <- deps_l_df[comp %in% compToKeep]
  }
  
  proteome_data_copy <- copy(proteome_data)
  deps_vec_all<- unique(deps_l_df[class != "=", id])
  
  if (type == "protein"){ 
    proteome_data_copy$dat_gene<-proteome_data_copy$dat_gene[GeneName %in% deps_vec_all,]
  } else if(type == "peptide") { 
    proteome_data_copy$dat_pep<-proteome_data_copy$dat_pep[ID_peptide %in% deps_vec_all,]
  } else{stop("type MUST be \"protein\" or \"peptide\"")}
  
  return(pca_plot(proteome_data = proteome_data_copy, type = type))
  
}
