#' From Proteome Discoverer or MaxQuant files to lists of data tables.
#'
#' This function reads phosphoproteomic data.
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
create_phosphosite_plot <- function(proteome_data, software) {
  psm_peptide_table <- copy(proteome_data$psm_peptide_table)
  
  if (software == "PD") {
    active_sites <- lapply(lapply(
      stri_extract_all_regex(psm_peptide_table$Modifications, "\\w\\d+\\(1|\\w\\d+\\.|\\ \\w\\]|\\w\\(1"), 
      function(x) { stri_replace_all(x, regex = " ", replacement = "") }), 
      function(y) { unlist(stri_extract_all_regex(unlist(y), "S|T|Y")) }
    ) |> unlist() |> table()
  } else if( software=="MQ"){
    active_sites <- lapply(
      stri_extract_all_regex(psm_peptide_table$Annotated_Sequence, "\\w\\(1"), 
      function(x) { stri_sub_all(x, from = 1, to = 1) }
    ) |> unlist() |> table()
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
                            text_col = "label") + theme(legend.position = "none", axis.text.y = element_text(size = 14))
  
  
  return(list("dt"=numeric_dt, "plot"=hs))
}
