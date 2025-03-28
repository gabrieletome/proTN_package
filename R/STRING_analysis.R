#' STRINGdb Network Analysis
#'
#' This function performs a STRINGdb network analysis on differentially expressed proteins or peptides.  
#' It processes input differential expression results and evaluates STRING network interactions.  
#'
#' @param differential_results A list containing differential expression results, either at the protein  
#'   or peptide level. It should include `protein_results_long` or `peptide_results_long`.
#' @param species Character. The species name (default: `"Homo sapiens"`).
#' @param dirOutput Character. The main output directory (default: `"results_ProTN"`).
#' @param subfold_Fig Character. Subfolder for figures (default: `"pics"`).
#' @param subfold_net Character. Subfolder for STRINGdb network results (default: `"STRINGdb"`).
#' @param phospho_ctrl Logical. Whether to include phospho-control comparisons (default: `FALSE`).
#'
#' @return A data frame with STRINGdb network results.
#' @export
#'
#' @import data.table
#' @import stringr
#' @import STRINGdb
#'
#' @examples
#' \dontrun{
#'   results <- STRINGdb_network(differential_results, species="Mus musculus")
#' }
STRINGdb_network <- function(differential_results, species="Homo sapiens", 
                             dirOutput="results_ProTN", subfold_Fig="pics", subfold_net="STRINGdb",
                             phospho_ctrl = FALSE){
  #Read Taxonomy
  codtax <- fread("R/NCBI_taxID/subset_tax.csv")
  tryCatch({
    taxonomy_NCBI <- species
    taxonomy_NCBI <- codtax[name == taxonomy_NCBI, taxid]
  }, error = function(cond){stop(paste0("Species not present in the db. Select one of the following: ",paste(codtax$name,collapse = ",")))})
  message(paste0("Taxonomy ID: ",taxonomy_NCBI))
  if(("protein_results_long" %in% names(differential_results))){
    phospho_with_proteome = FALSE
    deps_l_df <- copy(differential_results$protein_results_long)
  } else if(!("protein_results_long" %in% names(differential_results)) & 
            ("peptide_results_long" %in% names(differential_results))){
    phospho_with_proteome = TRUE
    deps_l_df <- copy(differential_results$peptide_results_long)
    deps_l_df[, id := tstrsplit(id, "_", keep = 1)[[1]]]
    
    if(!phospho_ctrl){
      compToKeep <- unique(grep("_Phospho_CTRL", deps_l_df$comp, value = T, invert = T))
      deps_l_df <- deps_l_df[comp %in% compToKeep]
    }
    
  } else{
    stop("Error in differential results paramenter! Verify the presence of protein_results_long or peptide_results_long")
  }
  
  res <- select_regulated_genes(deps_l_df)
  
  dirOutput_net=NULL
  if (length(res$doComp) == 0 | all(res$doComp == FALSE)) {
    stop("Not possible to continue with the STRING network evaluation. No protein up- or down-regulated.")
  } else {
    dir.create(file.path(dirOutput, subfold_Fig), showWarnings = FALSE)
    dirOutput_net <- paste(dirOutput, subfold_Fig, sep = "/")
    message(paste0("Saving in: ",dirOutput_net))
  }
  
  if (!is.null(dirOutput_net)) {
    stringdb_results <- process_string_network(g_sel_comp = res$g_sel_comp, dirOutput_net, taxonomy_NCBI)
    return(stringdb_results)
  }else{
    stop("Error in creating output directory")
  }
  
}

# Function to select up/down-regulated genes
select_regulated_genes <- function(deps_l_df) {
  g_sel_comp <- list()
  doComp <- list()
  for (comp in unique(deps_l_df$comp)) {
    tmp_deps <- deps_l_df[class != "=" & comp == comp]
    genes <- tmp_deps[order(p_val)][1:min(500, .N), .(id)]
    if (nrow(genes) > 0) {
      g_sel_comp[[comp]] <- genes
      doComp[[comp]] <- TRUE
    } else {
      g_sel_comp[[comp]] <- list()
    }
  }
  return(list(g_sel_comp = g_sel_comp, doComp = doComp))
}

# Function to process STRING network for each comparison
process_string_network <- function(g_sel_comp, dirOutput_net, taxonomy_NCBI) {
  string_db <- STRINGdb$new(version="12", species=taxonomy_NCBI, score_threshold=500, input_directory="R/STRINGdb/")
  
  stringdb_results <- list()
  for (comp in names(g_sel_comp)) {
    gene_name <- unique(g_sel_comp[[comp]])
    setnames(gene_name, "id", "gene_id")
    string_mapped <- string_db$map(gene_name, "gene_id", removeUnmappedRows = TRUE)
    links_string <- as.data.table(string_db$get_interactions(string_mapped$STRING_id))
    links_string[, `:=`(
      from = string_mapped$gene_id[match(from, string_mapped$STRING_id)],
      to = string_mapped$gene_id[match(to, string_mapped$STRING_id)]
    )]
    string_gene_df <- unique(data.table(gene1 = links_string$from, gene2 = links_string$to, weight = links_string$combined_score))
    if (nrow(string_gene_df) > 0) {
      fwrite(string_gene_df, file = paste0(dirOutput_net, "/", gsub(comp, pattern = "\\/", replacement="vs"), "_connection.txt"), sep = "\t")
      link <- paste0("https://string-db.org/cgi/network?identifiers=", paste(unique(string_mapped$STRING_id), collapse = "%0d"), "&species=", taxonomy_NCBI, "&required_score=500")
      string_db$plot_network(string_mapped, required_score = 500)
      pdf(file = paste0(dirOutput_net,"/", comp, "_network.pdf"), width = 150, height = 150)
      string_db$plot_network(string_mapped, required_score = 500)
      dev.off()
      
      stringdb_results[[comp]] <- string_gene_df
    } else {
      stop("No strong interaction detected between the proteins. Usually too few proteins.")
    }
  }
  return(stringdb_results)
}
