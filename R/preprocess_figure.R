#' Plot missing data
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
generate_abundance_plot <- function(proteome_data) {
  c_anno <- copy(proteome_data$c_anno)
  psm_sig_prot_df <- copy(proteome_data$psm_log_prot_df)
  psm_sig_prot_df[, ID_peptide := NULL]
  
  numeric_df <- c_anno[order(c_anno$condition)]
  numeric_df <- numeric_df[.(numeric_df$sample), on = "sample"]
  
  numeric_df[, numeric_values := colSums(is.na(psm_sig_prot_df))]
  numeric_df[, `% of available abundances` := 100 - (numeric_values / nrow(psm_sig_prot_df) * 100)]
  
  plot <- ggplot(data = numeric_df, aes(x = factor(sample, levels = unique(sample)), 
                                      y = `% of available abundances`, 
                                      fill = sample, colour = sample)) +
    coord_flip() +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.8) +
    theme_bw() +
    theme(legend.position = "none", axis.title.y = element_blank()) +
    scale_fill_manual(values = setNames(as.list(c_anno$color), c_anno$sample)) +
    scale_colour_manual(values = setNames(as.list(c_anno$color), c_anno$sample)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
  
  return(list("dt" = numeric_df, "plot" = plot))
}

#' generate_peptide_distribution_plot
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
generate_peptide_distribution_plot <- function(proteome_data) {
  psm_anno_df <- copy(proteome_data$psm_anno_df)
  numeric_df <- as.data.table(table(psm_anno_df$symbol))
  numeric_df[N > 20, N := 20]
  numeric_df <- numeric_df[order(N), .(NPeptides = .N), by = .(N)]
  numeric_df[, N := factor(str_sort(str_replace(as.character(N), pattern = "20", replacement = "20+"), numeric = TRUE), 
                           levels = c(1:19, "20+"))]
  
  plot <- ggplot(data = numeric_df, aes(x = N, y = NPeptides, color = NPeptides, fill = NPeptides)) +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.8) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "# peptides per protein", y = "# proteins") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank())
  
  return(list("dt" = numeric_df, "plot" = plot))
}


#' Function to generate abundance distribution plots
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
plot_abundance_distribution <- function(proteome_data, type) {
  if(type == "protein"){
    dat <- copy(proteome_data$dat_gene)
    c_anno <- copy(proteome_data$c_anno)
  } else if(type == "peptide"){
    dat <- copy(proteome_data$dat_pep)
    c_anno <- copy(proteome_data$c_anno)
  } else{
    stop("Invalid type parameter. Must be \"protein\" or \"peptide\"")
  }
  
  ordered_samples <- c_anno[order(condition), sample]
  pg_long_df <- melt(dat[, ..ordered_samples], variable.name = "sample", value.name = "log2 normalized abundance")
  pg_long_df[, id := rep(rownames(dat), length(ordered_samples))]
  setcolorder(pg_long_df, c("id", "sample", "log2 normalized abundance"))
  
  pg_long_df <- merge(pg_long_df, c_anno, by = "sample", all.x = TRUE)
  
  hs <- ggplot(data = pg_long_df, aes(x = factor(sample, levels = unique(pg_long_df$sample)), 
                                      y = `log2 normalized abundance`, fill = sample, colour = sample)) +
    coord_flip() +
    geom_violin(alpha = 0.5, scale = "width", trim = FALSE) +
    geom_boxplot(alpha = 1, fill = "white", width = 0.2, outlier.shape = NA, notch = FALSE) +
    theme_bw() +
    theme(legend.position = "none", axis.title.y = element_blank()) +
    scale_fill_manual(values = setNames(as.list(c_anno$color), c_anno$sample)) +
    scale_colour_manual(values = setNames(as.list(c_anno$color), c_anno$sample)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
  
  return(list("dt" = pg_long_df, "plot" = hs))
}

