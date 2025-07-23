#' Create Phosphosite Distribution Plot
#'
#' This function generates a lollipop plot showing the distribution of phosphorylated residues (S, T, Y)
#' in phosphoproteomics data, based on the output of either Proteome Discoverer (PD) or MaxQuant (MQ).
#'
#' @param proteome_data A list containing `psm_peptide_table`, which includes peptide modification information.
#' @param software Character. The software used for processing (`"PD"` for Proteome Discoverer, `"MQ"` for MaxQuant).
#' @param size_text Numeric. Size of the label in the plot. Default: `3`
#'
#' @return A list containing:
#'   - `dt`: A data table with phosphosite percentages.
#'   - `plot`: A ggplot2 object visualizing the phosphosite distribution.
#' @export
#'
#' @import data.table
#' @import stringi
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#'   phospho_plot <- create_phosphosite_plot(proteome_data, software="MQ")
#' }
create_phosphosite_plot <- function(proteome_data, software, size_text=3) {
  psm_peptide_table <- copy(proteome_data$psm_peptide_table)
  
  if (software == "PD") {
    active_sites <- table(unlist(lapply(lapply(
      stri_extract_all_regex(psm_peptide_table$Modifications, "\\w\\d+\\(1|\\w\\d+\\.|\\ \\w\\]|\\w\\(1"), 
      function(x) { stri_replace_all(x, regex = " ", replacement = "") }), 
      function(y) { unlist(stri_extract_all_regex(unlist(y), "S|T|Y")) }
    )))
  } else if( software=="MQ"){
    active_sites <- table(unlist(lapply(
      stri_extract_all_regex(psm_peptide_table$Annotated_Sequence, "\\w\\(1"), 
      function(x) { stri_sub_all(x, from = 1, to = 1) }
    )))
  } else{
    stop("Software must be MQ or PD")
  }
  
  numeric_dt <- data.table(
    Site = c("S", "T", "Y"),
    `Percentage of phosphosite (%)` = c(
      (active_sites["S"] / sum(active_sites) * 100),
      (active_sites["T"] / sum(active_sites) * 100),
      (active_sites["Y"] / sum(active_sites) * 100)
    ),
    label = c(
      paste0(round(active_sites["S"] / sum(active_sites) * 100, 2), "%"),
      paste0(round(active_sites["T"] / sum(active_sites) * 100, 2), "%"),
      paste0(round(active_sites["Y"] / sum(active_sites) * 100, 2), "%")
    )
  )
  
  # Create the plot
  hs <- enrichment_lollipop(numeric_dt, 
                            x_col = "Percentage of phosphosite (%)", 
                            y_col = "Site", 
                            size_col = "Percentage of phosphosite (%)", size_vec = c(6,6),
                            shape_col = "Site", shape_vec = c(16,16,16),
                            color_col = "Site", color_vec = c("#F79256", "#FBD1A2", "#7DCFB6"), 
                            fill_col = "Site", fill_vec = c("#F79256", "#FBD1A2", "#7DCFB6"), 
                            text_col = "label", size_text=size_text) + theme(legend.position = "none", axis.text.y = element_text(size = 14))
  
  
  return(list("dt"=numeric_dt, "plot"=hs))
}
