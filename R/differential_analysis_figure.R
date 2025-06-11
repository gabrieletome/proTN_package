#' Generate Differential Barplots
#'
#' This function generates differential barplots based on the differential analysis results. 
#' It creates a lollipop plot to visualize the comparison of protein or peptide differential expression results.
#'
#' @param differential_results A list containing the differential analysis results.
#' @param data_type A character string indicating the type of data to plot. Must be either "protein" or "peptide".
#' @param color_contrast A vector of colors to use for the plot (optional).
#' @param phospho_ctrl A logical value indicating whether to include phospho-control comparisons in the plot (default is FALSE).
#'
#' @return A list containing the lollipop plot for differential expression results.
#'
#' @details
#' This function filters the differential results based on the specified parameters, 
#' such as whether to include phospho-control comparisons or not, and generates a 
#' lollipop plot to visualize the differential expression of proteins or peptides.
#'
#' @examples
#' \dontrun{
#' plot_results <- generate_differential_barplots(differential_results = differential_results,
#'                                                data_type = "protein")
#' }
#' 
#' @import ggplot2
#' @importFrom plyr round_any
#' @importFrom dplyr pull select arrange
#' @import stringr
#' @import scales
#' @import RColorBrewer
#' @import data.table
#' @import ggthemes
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
      color_contrast=data.table("class"=c("+","-"), "color"=c("#664069","#8A628D"))
      col_dt <- merge.data.table(lolli_df[,c("comp","class","id")], color_contrast, by = "class")
      col_vec <- col_dt$color
      names(col_vec) <- col_dt$id
      message("Set default colors.")
    } else{
      col_vec <- as.vector(t(color_contrast))
    }
    
    pDEPs <- deps_b2b_lollipop(
      input_df = lolli_df,
      break_vec = seq(0, round_any(max(lolli_df$N) * 1.2, if (max(lolli_df$N) > 90) 100 else 10, f = ceiling),
                      round(round_any(max(lolli_df$N) * 1.2, if (max(lolli_df$N) > 90) 100 else 10, f = ceiling) / 5, 0)),
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


#' Generate Volcano Plots for Differential Expression Analysis
#'
#' This function generates volcano plots for visualizing the differential expression of proteins or peptides. 
#' It plots the log2 fold change (log2_FC) against the -log10 p-value (or FDR adjusted p-value),
#' highlighting significant changes based on fold-change and p-value thresholds.
#'
#' @param differential_results A list containing the differential analysis results.
#' @param data_type A character string indicating the type of data to plot. Must be either "protein" or "peptide" (default is NULL).
#' @param comparison A character string specifying the comparison to plot (e.g., "comparison1").
#' @param fc_thr A numeric value for the fold-change threshold (default is 0.75).
#' @param pval_fdr A character string indicating whether to use "p_val" or "p_adj" for the p-value threshold (default is "p_val").
#' @param pval_thr A numeric value for the p-value threshold (default is 0.05).
#' @param color_contrast A vector of two colors to use for plotting (optional).
#' @param interactomics Logical. Set to TRUE if is an interactomics analysis (default is FALSE).
#'
#' @return A list containing the volcano plot for the specified comparison.
#'
#' @details
#' This function generates volcano plots to visualize differential expression analysis, 
#' with fold-change and p-value (or FDR adjusted p-value) on the axes. 
#' It highlights genes/peptides based on the fold-change and p-value thresholds.
#'
#' @examples
#' \dontrun{
#' volcano_plot <- generate_volcano_plots(differential_results = differential_results,
#'                                        comparison = "comparison1",
#'                                        fc_thr = 1, pval_thr = 0.05)
#' }
#' 
#' @import tidyverse
#' @import data.table
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @import stringr
#' @import RColorBrewer
#' @export
generate_volcano_plots <- function(differential_results, data_type=NULL, 
                                   comparison=NULL, fc_thr=0.75, pval_fdr = "p_val", 
                                   pval_thr=0.05, color_contrast=NULL, interactomics = FALSE) {
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
    
    bs=16
    cmd <- ggplot(input_df, aes(x = log2_FC, y = get(col), col = class, text = paste('</br>Gene: ', id, '</br>Class: ', class, '</br>Log2_FC: ', log2_FC, '</br>', col, ': ', get(col)))) +
      geom_point(pch = 20, cex = 2) +
      geom_hline(yintercept = -log10(pval_thr), col = "black") +
      ggtitle(paste0("Volcano Plot of ", comparison)) +
      ylab(if (pval_fdr=="p_adj") "-log10(fdr)" else "-log10(p_val)") +
      scale_color_manual(values = c("+" = col_vec[1], "-" = col_vec[2], "=" = "grey70")) +
      scale_x_continuous(limits = c(min(-max(abs(input_df$log2_FC)), -3), max(max(abs(input_df$log2_FC)), 3))) +
      scale_y_continuous(limits = c(0, max(max(abs(input_df[,get(col)])), 4))) +
      theme_bw(base_size = bs)
    
    if(!interactomics){
      cmd <- cmd + geom_vline(xintercept = c(-fc_thr, fc_thr), col = "black")
    } else{
      cmd <- cmd + geom_vline(xintercept = c(fc_thr), col = "black")
    }
    
    plotlist[[comparison]] <- ggplotly(cmd, tooltip = c("text"))
    return(plotlist)
  } else{
    stop("No gene differentiated discovered.")
  }
}



#' MDS Plot
#'
#' This function generates a Multi-Dimensional Scaling (MDS) plot for either protein or peptide data, based on a Euclidean distance matrix of sample similarities.
#' It supports the inclusion of phosphoproteomics data if available.
#'
#' @param proteome_data A list containing proteome data, including sample annotation and either protein or peptide data matrices.
#' @param differential_analysis A list containing the differential analysis results.
#' @param type Character; specifies the type of data to plot. Options are "protein" or "peptide".
#'
#' @return A list containing:
#'   \item{mds_dt}{Data table containing the MDS coordinates for each sample.}
#'   \item{plot}{ggplot2 object representing the MDS plot with sample labels and color-coded conditions.}
#'
#' @examples
#' \dontrun{
#' result <- mds_plot(proteome_data, type = "protein")
#' print(result$plot)
#' }
#'
#' @import data.table
#' @import ggplot2
#' @import ggrepel
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


#' PCA Plot
#'
#' This function generates a Principal Component Analysis (PCA) plot for either protein or peptide data based on the first two principal components.
#' It supports the inclusion of phosphoproteomics data if available and provides a variance explanation for the axes.
#'
#' @param differential_analysis A list containing the differential analysis results.
#' @param proteome_data A list containing proteome data, including sample annotation and either protein or peptide data matrices.
#' @param type Character; specifies the type of data to plot. Options are "protein" or "peptide".
#'
#' @return A list containing:
#'   \item{pca_dt}{Data table containing the PCA coordinates for each sample.}
#'   \item{plot}{ggplot2 object representing the PCA plot.}
#'
#' @examples
#' \dontrun{
#' result <- pca_plot(proteome_data, type = "protein")
#' print(result$plot)
#' }
#'
#' @import data.table
#' @import ggplot2
#' @import ggrepel
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
  
  return(pca_plot(proteome_data = proteome_data_copy, type = type))
  
}
