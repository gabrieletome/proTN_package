#' Generate Abundance Plot
#'
#' This function generates an abundance plot.
#'
#' @param proteome_data A list containing `c_anno` or `c_anno_phospho`, `psm_log_prot_df`, and `psm_log_pet_df`.
#' @param phospho_with_proteome Logical. Whether to include phosphoproteome data alongside proteome data (default: `FALSE`).
#'
#' @return A list containing:
#'   - `dt`: Data table for abundance.
#'   - `plot`: A ggplot2 object visualizing abundance.
#' @export
#'
#' @import data.table
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#'   abundance_plot <- generate_abundance_plot(proteome_data)
#' }
generate_abundance_plot <- function(proteome_data, phospho_with_proteome=FALSE) {
  if(("c_anno" %in% names(proteome_data))){
    phospho_with_proteome = FALSE
  } else if(("c_anno_phospho" %in% names(proteome_data))){
    phospho_with_proteome = TRUE
  } else{
    stop("Missing sample annotation!")
  }
  
  result <- list()
  if(phospho_with_proteome){
    prot_list <- list("c_anno" = proteome_data$c_anno_proteome, "psm_log_prot_df" = proteome_data$psm_log_prot_df)
    res_prot <- generate_abundance_subplot(proteome_data = prot_list)
    phospho_list <- list("c_anno" = proteome_data$c_anno_phospho, "psm_log_prot_df" = proteome_data$psm_log_pet_df)
    res_phospho <- generate_abundance_subplot(proteome_data = phospho_list)
    result <- list("proteome_dt" = res_prot$dt, "proteome_plot" = res_prot$plot,
                   "phospho_dt" = res_phospho$dt, "phospho_plot" = res_phospho$plot)
  } else{
    result <- generate_abundance_subplot(proteome_data)
  }

  return(result)
}


generate_abundance_subplot <- function(proteome_data) {
  c_anno <- copy(proteome_data$c_anno)
  psm_sig_prot_df <- copy(proteome_data$psm_log_prot_df)
  psm_sig_prot_df[, ID_peptide := NULL]
  
  numeric_df <- c_anno[order(c_anno$condition)]
  numeric_df <- numeric_df[.(numeric_df$sample), on = "sample"]
  
  numeric_df[, numeric_values := colSums(is.na(psm_sig_prot_df))]
  numeric_df[, `Available` := 100 - (numeric_values / nrow(psm_sig_prot_df) * 100)]
  numeric_df[, `Missing` := 100 - `Available`]
  
  numeric_df <- melt(numeric_df, id.vars = c("sample"), measure.vars = c("Available","Missing"), variable.name = "Coverage", value.name = "% covered abundance")
  
  plot <- ggplot(data = numeric_df, aes(x = factor(sample, levels = unique(sample)), 
                                        y = `% covered abundance`, 
                                        fill = Coverage, colour = Coverage)) +
    coord_flip() +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.8) +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    scale_fill_manual(values = list("Available"="darkgreen", "Missing"="darkred")) +
    scale_colour_manual(values = list("Available"="darkgreen", "Missing"="darkred")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    ylim(0, 100)
  
  return(list("dt" = numeric_df, "plot" = plot))
}


#' Generate Peptide Distribution Plot
#'
#' This function generates a peptide distribution plot for a given proteome dataset.
#'

#' @param proteome_data A list containing proteome data, including peptide annotation tables.
#' @param phospho_with_proteome Logical; whether the dataset contains phosphoproteomics data alongside proteomics data. Default is FALSE.
#' 
#' @return A list containing:
#'   - `dt`: Data table for abundance distribution.
#'   - `plot`: A ggplot2 object visualizing abundance distribution.
#' @export
#'
#' @import data.table
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' result <- generate_peptide_distribution_plot(proteome_data)
#' }
generate_peptide_distribution_plot <- function(proteome_data, phospho_with_proteome=FALSE) {
  if(("c_anno" %in% names(proteome_data))){
    phospho_with_proteome = FALSE
  } else if(("c_anno_phospho" %in% names(proteome_data))){
    phospho_with_proteome = TRUE
  } else{
    stop("Missing sample annotation!")
  }
  
  result <- list()
  if(phospho_with_proteome){
    prot_list <- list("psm_anno_df" = proteome_data$psm_anno_df)
    res_prot <- generate_peptide_distribution_subplot(proteome_data = prot_list)
    phospho_list <- list("psm_anno_df" = proteome_data$psm_peptide_table)
    setnames(phospho_list$psm_anno_df, "GeneName", "symbol")
    res_phospho <- generate_peptide_distribution_subplot(proteome_data = phospho_list)
    setnames(proteome_data$psm_peptide_table, "symbol", "GeneName", skip_absent = T)
    result <- list("proteome_dt" = res_prot$dt, "proteome_plot" = res_prot$plot,
                   "phospho_dt" = res_phospho$dt, "phospho_plot" = res_phospho$plot)
  } else{
    result <- generate_peptide_distribution_subplot(proteome_data)
  }
  
  return(result)
}


generate_peptide_distribution_subplot <- function(proteome_data) {
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


#' Plot Abundance Distribution
#'
#' This function generates a violin and boxplot representation of log2-normalized abundance values for proteins or peptides
#' 
#' @param proteome_data A list containing proteome data, including peptide annotation tables.
#' @param type Character; specifies the type of data to plot. Options are "protein" or "peptide".
#'
#' @return A list containing:
#'   - `dt`: A data.table of log2-normalized abundance values with sample annotations
#'   - `plot`: A ggplot2 object representing the abundance distribution as a violin and boxplot.
#' @export
#'
#' @import data.table
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' result <- plot_abundance_distribution(proteome_data, type = "protein")
#' }
plot_abundance_distribution <- function(proteome_data, type) {
  if(("c_anno" %in% names(proteome_data))){
    phospho_with_proteome = FALSE
  } else if(("c_anno_proteome" %in% names(proteome_data)) & ("c_anno_phospho" %in% names(proteome_data))){
    phospho_with_proteome = TRUE
  } else{
    stop("Missing sample annotation!")
  }
  
  if(type == "protein"){
    dat <- copy(proteome_data$dat_gene)
    if(phospho_with_proteome){
      c_anno <- copy(proteome_data$c_anno_proteome)
    } else{
      c_anno <- copy(proteome_data$c_anno)
    }
    ordered_samples <- c_anno[order(condition), sample]
    pg_long_df <- melt(dat, id.vars = "GeneName",
                       variable.name = "sample", value.name = "log2 normalized abundance")
    setcolorder(pg_long_df, c("GeneName", "sample", "log2 normalized abundance"))
    
  } else if(type == "peptide"){
    dat <- copy(proteome_data$dat_pep)
    if(phospho_with_proteome){
      c_anno <- copy(proteome_data$c_anno_phospho)
    } else{
      c_anno <- copy(proteome_data$c_anno)
    }
    ordered_samples <- c_anno[order(condition), sample]
    pg_long_df <- melt(dat, id.vars = "ID_peptide",
                       variable.name = "sample", value.name = "log2 normalized abundance")
    setcolorder(pg_long_df, c("ID_peptide", "sample", "log2 normalized abundance"))
    
  } else{
    stop("Invalid type parameter. Must be \"protein\" or \"peptide\"")
  }
  
  pg_long_df <- merge(pg_long_df, c_anno, by = "sample", all.x = TRUE)
  
  hs <- ggplot(data = pg_long_df, aes(x = factor(sample, levels = unique(pg_long_df$sample)), 
                                      y = `log2 normalized abundance`, fill = sample, colour = sample)) +
    coord_flip() +
    geom_violin(alpha = 0.5, scale = "width", trim = FALSE) +
    geom_boxplot(alpha = 1, fill = "white", width = 0.2, outlier.shape = NA, notch = FALSE) +
    theme_bw() +
    theme(legend.position = "none", axis.title.y = element_blank()) +
    scale_fill_manual(values = setNames(c_anno$color, c_anno$sample)) +
    scale_colour_manual(values = setNames(c_anno$color, c_anno$sample)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
  
  return(list("dt" = pg_long_df, "plot" = hs))
}


#' MDS Plot
#'
#' This function generates a Multi-Dimensional Scaling (MDS) plot for either protein or peptide data, based on a Euclidean distance matrix of sample similarities.
#' It supports the inclusion of phosphoproteomics data if available.
#'
#' @param proteome_data A list containing proteome data, including sample annotation and either protein or peptide data matrices.
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
mds_plot <- function(proteome_data, type) {
  if(("c_anno" %in% names(proteome_data))){
    phospho_with_proteome = FALSE
  } else if(("c_anno_proteome" %in% names(proteome_data)) & ("c_anno_phospho" %in% names(proteome_data))){
    phospho_with_proteome = TRUE
  } else{
    stop("Missing sample annotation!")
  }
  
  if(type == "protein"){
    data_matrix <- copy(proteome_data$dat_gene)
    data_matrix <- as.matrix(data_matrix[,-1])
    if(phospho_with_proteome){
      c_anno <- copy(proteome_data$c_anno_proteome)
    } else{
      c_anno <- copy(proteome_data$c_anno)
    }
  } else if(type == "peptide"){
    data_matrix <- copy(proteome_data$dat_pep)
    data_matrix <- as.matrix(data_matrix[,-1])
    if(phospho_with_proteome){
      c_anno <- copy(proteome_data$c_anno_phospho)
    } else{
      c_anno <- copy(proteome_data$c_anno)
    }
  } else{
    stop("Invalid type parameter. Must be \"protein\" or \"peptide\"")
  }
  
  sample_distances <- dist(t(data_matrix), method = "euclidean")
  mds_cmdscale <- as.data.table(cmdscale(as.matrix(sample_distances)), keep.rownames = "sample")
  setnames(mds_cmdscale, c("V1", "V2"), c("MDS_1", "MDS_2"))
  mds_cmdscale <- merge(mds_cmdscale, c_anno, by = "sample", all.x = TRUE)
  
  cc <- unique(mds_cmdscale$color)
  names(cc) <- unique(mds_cmdscale$condition)
  
  cmd <- ggplot(mds_cmdscale, aes(MDS_1, MDS_2, colour = condition)) +
    geom_point(size = 2, alpha = .9) +
    geom_text_repel(aes(label = sample), size = 4, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = cc) +
    theme_bw() +
    theme(legend.position = "right", panel.grid.minor = element_blank())
  
  return(list("mds_dt" = mds_cmdscale, "plot" = cmd))
}


#' PCA Plot
#'
#' This function generates a Principal Component Analysis (PCA) plot for either protein or peptide data based on the first two principal components.
#' It supports the inclusion of phosphoproteomics data if available and provides a variance explanation for the axes.
#'
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
pca_plot <- function(proteome_data, type) {
  if(("c_anno" %in% names(proteome_data))){
    phospho_with_proteome = FALSE
  } else if(("c_anno_proteome" %in% names(proteome_data)) & ("c_anno_phospho" %in% names(proteome_data))){
    phospho_with_proteome = TRUE
  } else{
    stop("Missing sample annotation!")
  }
  
  if(type == "protein"){
    data_matrix <- copy(proteome_data$dat_gene)
    data_matrix <- as.matrix(data_matrix[,-1])
    if(phospho_with_proteome){
      c_anno <- copy(proteome_data$c_anno_proteome)
    } else{
      c_anno <- copy(proteome_data$c_anno)
    }
  } else if(type == "peptide"){
    data_matrix <- copy(proteome_data$dat_pep)
    data_matrix <- as.matrix(data_matrix[,-1])
    if(phospho_with_proteome){
      c_anno <- copy(proteome_data$c_anno_phospho)
    } else{
      c_anno <- copy(proteome_data$c_anno)
    }
  } else{
    stop("Invalid type parameter. Must be \"protein\" or \"peptide\"")
  }
  
  apca_prot <- prcomp(t(data_matrix), scale. = TRUE, center = TRUE)
  pc <- as.data.table(apca_prot$x[, 1:2], keep.rownames = "sample")
  pc <- merge(pc, c_anno, by = "sample", all.x = TRUE)
  
  ve <- (apca_prot$sdev^2) / sum(apca_prot$sdev^2)
  labs <- paste0(names(pc)[2:3], " (", round(ve[1:2] * 100, 2), "%)")
  
  cc <- unique(pc$color)
  names(cc) <- unique(pc$condition)
  
  cmd <- ggplot(pc, aes(PC1, PC2, colour = condition)) +
    geom_point(size = 2, alpha = .9) +
    geom_text_repel(aes(label = sample), size = 4, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = cc) +
    theme_bw() +
    theme(legend.position = "right", panel.grid.minor = element_blank()) +
    xlab(labs[1]) +
    ylab(labs[2])
  
  return(list("pca_dt" = pc, "plot" = cmd))
}


#' Plot Selected Proteins
#'
#' This function generates a plot for the selected proteins' abundance across different conditions. It creates a crossbar plot 
#' with error bars representing the average intensity and standard error of selected proteins. If phosphoproteomics data is available, 
#' it uses the corresponding data, otherwise, it uses the proteomics data.
#'
#' @param proteome_data A list containing proteome data, including gene intensity data and sample annotations.
#' @param list_protein A character vector of gene names to be plotted.
#'
#' @return A list containing:
#'   \item{dt}{Data table of protein intensity values and the corresponding sample and condition information.}
#'   \item{plot}{ggplot2 object representing the protein abundance plot for the selected proteins.}
#'
#' @examples
#' \dontrun{
#' result <- plot_selected_proteins(proteome_data, list_protein = c("ProteinA", "ProteinB"))
#' print(result$plot)
#' }
#'
#' @import data.table
#' @import ggplot2
#' @import ggbeeswarm
#' @import stringr
#' @export
plot_selected_proteins <- function(proteome_data, list_protein) {
  if(("c_anno" %in% names(proteome_data))){
    phospho_with_proteome = FALSE
  } else if(("c_anno_phospho" %in% names(proteome_data))){
    phospho_with_proteome = TRUE
  } else{
    stop("Missing sample annotation!")
  }
  
  if(phospho_with_proteome){
    dat_gene <- copy(proteome_data$dat_pep)
    dat_gene[, GeneName := tstrsplit(ID_peptide, "_", keep = 1)[[1]]][, ID_peptide := NULL]
    c_anno = copy(proteome_data$c_anno_phospho)
  } else{
    dat_gene <- copy(proteome_data$dat_gene)
    c_anno = copy(proteome_data$c_anno)
  }
  
  prot_find <- unique(c(
    intersect(list_protein, dat_gene$GeneName),
    intersect(str_to_title(list_protein), dat_gene$GeneName),
    intersect(str_to_upper(list_protein), dat_gene$GeneName),
    intersect(str_to_lower(list_protein), dat_gene$GeneName)
  ))
  
  if (length(prot_find) > 0) {
    message("Selected proteins with available abundances: \n")
    message(paste(prot_find, collapse = ", "))

    prot_intensity_long <- melt(dat_gene[GeneName %in% prot_find,], id.vars = "GeneName", variable.name = "sample", value.name = "Intensity")
    prot_intensity_long <- merge.data.table(prot_intensity_long, 
                                            c_anno,
                                            by = "sample")
    
    prot_avg_se_long <- prot_intensity_long[, .(avg = mean(Intensity),
                                                se = sd(Intensity)/sqrt(.N)), by = c("GeneName", "condition")]
    
    cc <- unique(prot_intensity_long$color)
    names(cc) <- unique(prot_intensity_long$condition)
    
    bs = 23
    g<-ggplot(prot_avg_se_long,aes(condition,avg,fill=condition,colour=condition))+
      geom_crossbar(aes(ymin=avg,ymax=avg),position = "dodge",width=.8,alpha=.9,fatten=1.5)+
      geom_errorbar(aes(ymin=(avg-se), ymax=(avg+se)), width=.4,position=position_dodge(),show.legend=F,alpha=.8)+
      geom_quasirandom(data=prot_intensity_long, aes(condition,Intensity), alpha=.7,width=.1,shape=16,size=0.2*bs)+
      scale_fill_manual(name="condition",values=cc[sort(unique(prot_intensity_long$condition))]) +
      scale_colour_manual(name="condition",values=cc[sort(unique(prot_intensity_long$condition))])+
      theme_bw(base_size = bs) +
      theme(axis.title.x=element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1),
            panel.grid.major.x=element_blank(),
            panel.grid.minor.y=element_blank(),
            legend.text = element_text(size = 0.7*bs),
            legend.key.size = unit((0.015*bs),"in"),
            legend.position="none",
            legend.title=element_blank(),
            legend.background = element_rect(fill = NA),
            strip.text=element_text(colour="white",face="bold"),
            panel.border=element_rect(colour=c("grey40"),size=0.03*bs),
            strip.background=element_rect(fill="grey40",colour="grey40",size=0.03*bs),
            plot.title = element_text(hjust = 0.5))+
      facet_wrap(~GeneName, scales = "free",ncol = if(length(prot_find)>4){round(length(prot_find)/1.9)}else{4})+
      labs(y="Abundance")

    return(list("dt" = prot_intensity_long, "plot" = g))
  } else {
    warning("Selected proteins not found in the normalized matrix. Check the spell of the proteins.")
    return(NULL)
  }
}
