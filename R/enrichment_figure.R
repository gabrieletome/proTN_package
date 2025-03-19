#' Enrichment Analysis figure
#'
#' @param proteome_data Object proTN
#' @param pval_fdr "p_adj" else "p_val"
#' @return A list of data tables.
#' @details \strong{ProTN}
#' @examples 
#' ## ## Example:
#' ## example
#' ## example2
#' @import data.table
#' @export
enrichment_figure <- function(enr_df, category = c("all"), enrich_filter_term=NULL, overlap_size_enrich_thr = 5, save=F,
                              color_contrast=NULL, dirOutput="results_ProTN", subfolder="pics", namefile="enrichement_plot"){
  # Esecuzione delle funzioni
  enr_sele_df <- enrichment_filter(enr_df, category, overlap_size_enrich_thr, enrich_filter_term)
  category_db <- generate_category_db(enr_sele_df)
  plotlist <- generate_plotlist(enr_sele_df, category_db, color_contrast)
  if(save){save_plotlist(plotlist, enr_sele_df, category_db, dirOutput, subfolder, namefile)}
  
  return(plotlist)
}


# Funzione per filtrare i dati di arricchimento
enrichment_filter <- function(enr_df, category = c("all"), overlap_size_enrich_thr = 5, enrich_filter_term=NULL) {
  if(!all(category %in% c("all","up","down"))){stop("category not valid. MUST be a vector containing all, down, up")}
  
  if(!is.null(enrich_filter_term)){
    message("Filtering by terms...")
    lookup_words<-enrich_filter_term
    
    enr_sele_names<-NULL
    #For each pattern search in DEPs name
    for(word in lookup_words){
      matches <- enr_df[,c("anno_class","anno_name")]
      matches <- matches[grepl(paste0("\\b",word,"\\b"),matches$anno_name,ignore.case=TRUE),]
      enr_sele_names<-unique(rbind(enr_sele_names,matches))
    }
    
    enr_df <- enr_df[anno_name %in% enr_sele_names$anno_name]
  }
  
  if(nrow(enr_df)==0 & !is.null(enrich_filter_term)){
    stop("Zero terms detectect with the words in enrich_filter_term")
  }
  
  if(length(category)>1){
    enr_sele_df <- enr_df[overlap_size >= overlap_size_enrich_thr & 
                            anno_size <= 10000 & 
                            Significant == "TRUE"]
    enr_sele_df_tmp <- data.table()
    for(cc in category){
      enr_sele_df_tmp <- rbind(enr_sele_df_tmp, enr_sele_df[grepl(paste0("_",cc,"$"), input_name)])
    }
    enr_sele_df <- enr_sele_df_tmp                   
  }else{
    enr_sele_df <- enr_df[overlap_size >= overlap_size_enrich_thr & 
                          anno_size <= 10000 & 
                          Significant == "TRUE" & 
                          grepl(paste0("_",category,"$"), input_name)]
  }
  enr_sele_df <- enr_sele_df[, .SD[overlap_size %in% tail(sort(unique(overlap_size)), 10)], by=c("input_name","anno_class")]
  setorder(enr_sele_df, -overlap_size)
  enr_sele_names <- unique(enr_sele_df[, .(anno_class, anno_name)])
  enr_sele_df <- enr_sele_df[enr_sele_names, on = .(anno_class, anno_name)]
  return(enr_sele_df)
}

# Funzione per generare le liste di categorie
generate_category_db <- function(enr_sele_df) {
  if (!is.null(enr_sele_df)) {
    # LOAD category EnrichR
    dbs_default <- fread("R/enrichR/dbs_enrichR.txt", header = F)
    dbs_category <- split(dbs_default, f = as.factor(dbs_default$V2))
    
    category_db <- lapply(dbs_category, function(x) {
      x[x[[1]] %in% unique(enr_sele_df$anno_class)]
    })
    return(category_db)
  } else{
    stop("No significant terms detected.")
  }
}

# Funzione per generare la lista di plot
generate_plotlist <- function(enr_sele_df, category_db, color_contrast=NULL) {
  plotlist <- list()
  
  if(is.null(color_contrast)){
    color_contrast=c("#664069","#8A628D")
    col_vec <- rep(as.vector(t(color_contrast)), uniqueN(enr_sele_df$input_name))
    message("Set default colors.")
  } else{
    col_vec <- as.vector(t(color_contrast))
  }
  
  tryCatch({
    names(col_vec) <- unique(enr_sele_df$input_name)
  }, error = function(cond){
    stop("Color must be a vector of the same length of the number of comparison.")
  })
  
  for (db in names(category_db)) {
    if (nrow(category_db[[db]]) > 0) {
      enr_sele_df_db <- enr_sele_df[anno_class %in% category_db[[db]][[1]]]
      plotlist[[db]] <- enrichment_dotmatrix(
        enr_sele_df_db,
        size_col = "log2_OR",
        color_col = "input_name",
        color_vec = col_vec,
        shape_vec = c(16, 21),
        shape_col = "Significant",
        fill_col = "Significant",
        char_max = 60
      ) + guides(color = FALSE, shape = guide_legend(override.aes = list(size = 4)))
    }
  }
  return(plotlist)
}

# Funzione per generare e salvare i plot
save_plotlist <- function(plotlist, enr_sele_df, category_db, dirOutput="results_ProTN", subfolder="pics", namefile="enrichment_plot") {
  tryCatch({
    name_list <- vector()
    for (db in names(plotlist)) {
      name_list[db] <- paste0(dirOutput,"/",subfolder,"/", "enr_DE_keysources_", db, ".pdf")
      ggsave(name_list[db], plotlist[[db]], 
             device = cairo_pdf, width = 7, 
             height = max(min(20, length(unique(enr_sele_df[anno_class %in% category_db[[db]][[1]], anno_name])) * 0.3), 3),
             units = "in", create.dir = T)
    }
    pdf_combine(input = name_list, output = paste0(dirOutput,"/",subfolder,"/", namefile, ".pdf"))
    unlink(name_list)
  }, error = function(cond) {
    stop("Error: No enriched term found with current parameters\n")
  })
}
