#' Perform Enrichment Analysis on Differential Results
#'
#' This function performs enrichment analysis using the EnrichR tool based on the results
#' of differential expression analysis (proteomics or phospho-proteomics). The function
#' generates enrichment results and saves them as both an Excel file and an RData object
#' for further analysis.
#'
#' It supports custom EnrichR databases and allows the user to filter enrichment results 
#' based on p-value, FDR threshold, and overlap size.
#'
#' The function generates two types of output files:
#' - `enrichment.xlsx`: An Excel file containing the enrichment results.
#' - `enrichment.RData`: A saved RData file with the filtered enrichment results for further analysis.
#'
#' @param differential_results A list containing differential results, either `protein_results_long`
#'   for proteomics data or `peptide_results_long` for phospho-proteomics data.
#' @param dirOutput The directory where the results file will be saved. Default is "results_ProTN".
#' @param subfold_Tab The subfolder within `dirOutput` where the Excel file will be saved. Default is "table".
#' @param subfold_Dat The subfolder within `dirOutput` where the RData file will be saved. Default is "rdata".
#' @param pval_fdr_enrich The column name in the differential results containing the p-value after FDR correction.
#'   Default is `"p_adj"`.
#' @param pval_enrich_thr The threshold for the p-value (adjusted) to consider significant enrichment.
#'   Default is 0.05.
#' @param overlap_size_enrich_thr The threshold for the size of the overlap (number of DEPs in the term) 
#'   to consider significant enrichment. Default is 5.
#' @param enrichR_custom_DB A logical flag to indicate if custom EnrichR databases should be used. Default is `FALSE`.
#' @param enrich_filter_DBs A vector of custom EnrichR databases to filter enrichment results.
#'   This argument is used if `enrichR_custom_DB` is set to `TRUE`.
#' @param phospho_ctrl A logical flag to indicate if phospho control data should be excluded.
#'   Default is `FALSE`.
#'
#' @return A data.table containing the enrichment results, including terms, significance scores, 
#'   overlap sizes, and other relevant information.
#'
#' @import writexl
#' @import data.table
#' @import enrichR
#' @import parallel
#' @importFrom readr read_tsv
#' @importFrom plyr ldply
#' @importFrom tidyr separate_rows
#' @importFrom dplyr summarise rename select arrange group_by ungroup mutate left_join
#' @importFrom curl has_internet
#' @importFrom httr GET POST status_code http_status use_proxy
#'
#' @examples
#' \dontrun{
#' perform_enrichment_analysis(differential_results = differential_results_example,
#'                             dirOutput = "results_directory")
#' }
#' @export
perform_enrichment_analysis <- function(differential_results, dirOutput="results_ProTN", 
                                        subfold_Tab="table", subfold_Dat="rdata",
                                        pval_fdr_enrich="p_adj", pval_enrich_thr=0.05, 
                                        overlap_size_enrich_thr=5, enrichR_custom_DB=FALSE, enrich_filter_DBs=NULL,
                                        phospho_ctrl = FALSE) {
  
  source("R/enrichR/functions.R")
  set_enrichR()
  if(("protein_results_long" %in% names(differential_results))){
    phospho_with_proteome = FALSE
    diff_dt <- differential_results$protein_results_long
  } else if(("peptide_results_long" %in% names(differential_results))){
    phospho_with_proteome = TRUE
    diff_dt <- copy(differential_results$peptide_results_long)
    diff_dt[, id := tstrsplit(id, "_", keep = 1)[[1]]]
    
    if(!phospho_ctrl){
      compToKeep <- unique(grep("_Phospho_CTRL", diff_dt$comp, value = T, invert = T))
      diff_dt <- diff_dt[comp %in% compToKeep]
    }
    
  } else{
    stop("Error in differential results paramenter! Verify the presence of protein_results_long or peptide_results_long")
  }
  
  doNextChunk <- tryCatch({
    dbs <- NULL
    if (enrichR_custom_DB) {
      dbs <- enrich_filter_DBs
    }
    
    # Ensure directories exist
    dir.create(file.path(dirOutput, subfold_Dat), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(dirOutput, subfold_Tab), showWarnings = FALSE, recursive = TRUE)

    # Perform enriched analysis with EnrichR
    enr_df <- suppressMessages(enrichRfnc(in_df = diff_dt, pval_fdr_enrich, pval_enrich_thr, overlap_size_enrich_thr, dbs))
    enr_df <- as.data.table(enr_df)
    # Save in RData for possible further analysis
    enrich_df <- enr_df[overlap_size >= overlap_size_enrich_thr, .SD, .SDcols = 1:13]
    save(enrich_df, file = paste0(dirOutput, "/", subfold_Dat, "/", "enrichment.RData"))
    
    # Prepare the README
    readme_sheet <- data.table(INFO = c(
      NA,
      "Excel file containing a selection of enrichment results starting from differentially expressed proteins. Terms are selected according to significance thresholds specified in the input (Default: adj.P.Value < 0.05, Overlap Size >= 5)",
      NA,
      "Enrichment columns:",
      "1. *input_name*: comparison name,",
      "2. *anno_name*: enriched name term,",
      "3. *anno_class*: dataset,",
      "4. *overlap_size*: DEPs in term,",
      "5. *p_value*: p-value,",
      "6. *fdr*: adjusted p-value (FDR after BH correction),",
      "7. *odds_ratio*: ",
      "8. *combined_score*: combined score provided by EnrichR,",
      "9. *input_size*: DEPs of the comparison,",
      "10. *anno_size*: number protein of the term,",
      "11. *overlap_input_ratio*: overlap_size/input_size",
      "12. *overlap_anno_ratio*: overlap_size/anno_size",
      "13. *overlap_ids*: gene symbols identified in the term"
    ))
    
    write_xlsx(enrich_df, path = paste0(dirOutput, "/", subfold_Tab, "/", "enrichment.xlsx"))
    
    return(enr_df)
  }, error = function(cond) {
    stop(paste0("An error occurred when connecting to EnrichR. \n ",as.character(cond)))
  })
}
