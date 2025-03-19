#' Enrichment Analysis
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
perform_enrichment_analysis <- function(differential_results, dirOutput="results_ProTN", 
                                        subfold_Tab="table", subfold_Dat="rdata",
                                        pval_fdr_enrich="p_adj", pval_enrich_thr=0.05, 
                                        overlap_size_enrich_thr=5, enrichR_custom_DB=FALSE, enrich_filter_DBs=NULL) {
  doNextChunk <- tryCatch({
    dbs <- NULL
    if (enrichR_custom_DB) {
      dbs <- enrich_filter_DBs
    }
    
    # Ensure directories exist
    dir.create(file.path(dirOutput, subfold_Dat), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(dirOutput, subfold_Tab), showWarnings = FALSE, recursive = TRUE)

    # Perform enriched analysis with EnrichR
    enr_df <- enrichRfnc(in_df = differential_results$protein_results_long, pval_fdr_enrich, pval_enrich_thr, overlap_size_enrich_thr, dbs)
    enr_df <- as.data.table(enr_df)
    # Save in RData for possible further analysis
    enrich_df <- enr_df[overlap_size >= overlap_size_enrich_thr, .SD, .SDcols = 1:13]
    save(enrich_df, file = paste0(dirOutput, subfold_Dat, "enrichment.RData"))
    
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
    
    writexl::write_xlsx(enrich_df, path = paste0(dirOutput, "/", subfold_Tab, "/", "enrichment.xlsx"))
    
    return(enr_df)
  }, error = function(cond) {
    stop("An error occurred when connecting to EnrichR. \n ")
  })
}
