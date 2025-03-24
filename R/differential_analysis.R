#' Differential Analysis
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
differential_analysis <- function(proteome_data, formule_contrast, 
                                  fc_thr=0.75, pval_fdr = "p_val", pval_thr=0.05, 
                                  signal_thr=0) {
  if(("c_anno" %in% names(proteome_data))){
    phospho_with_proteome = FALSE
  } else if(("c_anno_phospho" %in% names(proteome_data))){
    phospho_with_proteome = TRUE
  } else{
    stop("Missing sample annotation!")
  }
  
  message("**Thresholds used:**\n")
  message(sprintf("- Fold Change threshold: log2_FC > %s (+), log2_FC < -%s (-)\n", fc_thr, fc_thr))
  message(sprintf("- Statistical significance threshold: p.val < %s\n", pval_thr))
  if (signal_thr != -Inf) message(sprintf("- Signal log2 intensity threshold: signal > %s\n", signal_thr))
  
  if(phospho_with_proteome){
    # Merging annotation
    c_anno_prot <- copy(proteome_data$c_anno_proteome)
    c_anno_phos <- copy(proteome_data$c_anno_phospho)
    c_anno_prot$sample <- str_c(c_anno_prot$sample, "_proteome")
    c_anno_phos$sample <- str_c(c_anno_phos$sample, "_phospho")
    c_anno_prot$condition <- str_c(c_anno_prot$condition, "_proteome")
    c_anno_phos$condition <- str_c(c_anno_phos$condition, "_phospho")
    c_anno <- rbind(c_anno_phos,c_anno_prot)
    
    # Change name samples
    dat_gene <- copy(proteome_data$dat_gene)
    setnames(dat_gene, old = proteome_data$c_anno_proteome$sample, 
             new = str_c(proteome_data$c_anno_proteome$sample, "_proteome"))
    dat_pep <- copy(proteome_data$dat_pep)
    setnames(dat_pep, old = proteome_data$c_anno_proteome$sample, 
             new = str_c(proteome_data$c_anno_proteome$sample, "_phospho"))
    dat_pep[, GeneName := tstrsplit(ID_peptide, "_", keep = 1)[[1]]]

    # Merge proteomic and phosphoproteomic
    dat_complete <- dat_pep[dat_gene, on="GeneName", nomatch = NULL][, GeneName := NULL]
    
    psm_count_table <- data.table(table(dat_complete$ID_peptide))
    
    # Process proteins
    expr_mat <- copy(dat_gene)
    setnames(expr_mat, "GeneName", "gene")
    expr_l_df <- melt(expr_mat, id.vars = "gene", variable.name = "sample", value.name = "expr")
    expr_l_df <- expr_l_df[c_anno_prot, on = "sample"]
    expr_cond_df <- expr_l_df[, .(N = .N, avg = mean(expr), sd = sd(expr), CV = sd(expr) / mean(expr), se = sd(expr) / sqrt(.N)), by = .(condition, gene)]
    
    expr_avg_df <- dcast(expr_cond_df[, .(gene, condition, avg)], gene ~ condition, value.var = "avg")
    colnames(expr_avg_df)[-1] <- paste0(colnames(expr_avg_df)[-1], "_avg")
    expr_se_df <- dcast(expr_cond_df[, .(gene, condition, se)], gene ~ condition, value.var = "se")
    colnames(expr_se_df)[-1] <- paste0(colnames(expr_se_df)[-1], "_se")
    expr_cv_df <- dcast(expr_cond_df[, .(gene, condition, CV)], gene ~ condition, value.var = "CV")
    colnames(expr_cv_df)[-1] <- paste0(colnames(expr_cv_df)[-1], "_Coef_Variant_(%)")
    
    expr_avgse_df <- Reduce(merge, list(expr_avg_df, expr_se_df, expr_cv_df))
    colnames(expr_avgse_df)[1] <- "GeneName"
    # Process peptides
    expr_mat <- copy(dat_pep)
    expr_mat[, GeneName := NULL]
    expr_l_df <- melt(expr_mat, id.vars = "ID_peptide", variable.name = "sample", value.name = "expr")
    expr_l_df <- expr_l_df[c_anno_phos, on = "sample"]
    expr_cond_df <- expr_l_df[, .(N = .N, avg = mean(expr), sd = sd(expr), CV = sd(expr) / mean(expr), se = sd(expr) / sqrt(.N)), by = .(condition, ID_peptide)]
    
    expr_avg_pep_df <- dcast(expr_cond_df[, .(ID_peptide, condition, avg)], ID_peptide ~ condition, value.var = "avg")
    colnames(expr_avg_pep_df)[-1] <- paste0(colnames(expr_avg_pep_df)[-1], "_avg")
    expr_se_pep_df <- dcast(expr_cond_df[, .(ID_peptide, condition, se)], ID_peptide ~ condition, value.var = "se")
    colnames(expr_se_pep_df)[-1] <- paste0(colnames(expr_se_pep_df)[-1], "_se")
    expr_cv_pep_df <- dcast(expr_cond_df[, .(ID_peptide, condition, CV)], ID_peptide ~ condition, value.var = "CV")
    colnames(expr_cv_pep_df)[-1] <- paste0(colnames(expr_cv_pep_df)[-1], "_Coef_Variant_(%)")
    
    expr_avgse_pep_df <- Reduce(merge, list(expr_avg_pep_df, expr_se_pep_df, expr_cv_pep_df))
    expr_avgse_pep_df[, GeneName := tstrsplit(ID_peptide, "_", keep = 1)[[1]]]
    
    # Merge proteomic and phosphoproteomic
    expr_avgse_complete_df <- expr_avgse_pep_df[expr_avgse_df, on="GeneName", nomatch = NULL][, GeneName := NULL]
    
    #Make the design
    df <- data.frame("phospho" = unique(sort(str_remove(c_anno_phos$condition, "_phospho\\b"))), 
                     "proteome" = unique(sort(str_remove(c_anno_prot$condition, "_proteome\\b"))))
    df <- unique(df[match(df$phospho, df$proteome),])
    df$rule <- str_c(str_c(df$phospho, "_phospho"), "-", str_c(df$proteome, "_proteome"))
    contro_list <- df$rule
    
    contro_list<-c(contro_list,
                   str_c(str_c("(", stringi::stri_replace_all_regex(formule_contrast, 
                                                                    str_c("\\b",df$phospho,"\\b"), 
                                                                    str_c(df$phospho, "_phospho"), vectorize = FALSE), ")"),
                         "-",
                         str_c("(", stringi::stri_replace_all_regex(formule_contrast, 
                                                                    str_c("\\b",df$phospho,"\\b"), 
                                                                    str_c(df$phospho, "_proteome"), vectorize = FALSE), ")"))
    )
    names(contro_list) <- make.names(c(str_c(df$phospho,"_Phospho_CTRL"), names(formule_contrast)), unique = T)
    
    if (length(contro_list) == 0) {
      stop("No valid contrast design given. Check the match between the spell of Condition and contrast design.")
    }
    
    message("The following table contains the names and formulas of the contrasts considered for differential expression analysis:\n")
    message(as.character(contro_list))
    
    deps_pep_df <- limmafnc_dt("PEP", c_anno = c_anno, dat_gene = dat_complete, psm_count_table = psm_count_table, formule_contrast = contro_list,
                               expr_avgse_df = expr_avgse_complete_df, signal_thr = signal_thr, fc_thr = fc_thr, pval_thr = pval_thr, pval_fdr = pval_fdr)
    
    deps_pep_l_df <- deps_pep_df$degs_l_df
    deps_pep_w_df <- deps_pep_df$degs_w_df
    
    return(list(peptide_results_long = deps_pep_l_df, peptide_results_wide = deps_pep_w_df))
  } else{
    # Differential analysis proteomic or phospho without proteome
    c_anno <- copy(proteome_data$c_anno)
    
    design <- model.matrix(~0 + c_anno$condition)
    colnames(design) <- levels(as.factor(c_anno$condition))
    rownames(design) <- c_anno$sample
    
    filt_contro_list <- list()
    for (i in seq_along(formule_contrast)) {
      if (all(stri_remove_empty(str_remove_all(str_extract_all(formule_contrast[i], "\\w+")[[1]], "^\\d+$")) %in% colnames(design))) {
        filt_contro_list <- c(filt_contro_list, i)
      }
    }
    
    contro_list <- formule_contrast[unlist(filt_contro_list)]
    if (length(contro_list) == 0) {
      stop("No valid contrast design given. Check the match between the spell of Condition and contrast design.")
    }
    
    contrast <- limma::makeContrasts(contrasts = contro_list, levels = design)
    colnames(contrast) <- names(contro_list)
    
    unique_cond <- colSums(design) == 1
    unique_cond <- sort(names(unique_cond[unique_cond]))
    res <- sapply(colnames(contrast), function(i) {
      all(sort(names((contrast[, i] != 0)[contrast[, i] != 0])) %in% unique_cond)
    })
    formule_contrast <- contro_list[!res]
    if (length(formule_contrast) == 0) {
      stop("No valid contrast design given. Condition with only 1 replica")
    }
    
    message("The following table contains the names and formulas of the contrasts considered for differential expression analysis:\n")
    message(contro_list)
    
    message("Differentiation analysis in progress...")
    dat_gene <- copy(proteome_data$dat_gene)
    psm_count_table <- data.table(table(proteome_data$psm_anno_df$symbol))
    
    # Process proteins
    expr_mat <- copy(proteome_data$dat_gene)
    setnames(expr_mat, "GeneName", "gene")
    expr_l_df <- melt(expr_mat, id.vars = "gene", variable.name = "sample", value.name = "expr")
    expr_l_df <- expr_l_df[c_anno, on = "sample"]
    expr_cond_df <- expr_l_df[, .(N = .N, avg = mean(expr), sd = sd(expr), CV = sd(expr) / mean(expr), se = sd(expr) / sqrt(.N)), by = .(condition, gene)]
    
    expr_avg_df <- dcast(expr_cond_df[, .(gene, condition, avg)], gene ~ condition, value.var = "avg")
    colnames(expr_avg_df)[-1] <- paste0(colnames(expr_avg_df)[-1], "_avg")
    expr_se_df <- dcast(expr_cond_df[, .(gene, condition, se)], gene ~ condition, value.var = "se")
    colnames(expr_se_df)[-1] <- paste0(colnames(expr_se_df)[-1], "_se")
    expr_cv_df <- dcast(expr_cond_df[, .(gene, condition, CV)], gene ~ condition, value.var = "CV")
    colnames(expr_cv_df)[-1] <- paste0(colnames(expr_cv_df)[-1], "_Coef_Variant_(%)")
    
    expr_avgse_df <- Reduce(merge, list(expr_avg_df, expr_se_df, expr_cv_df))
    colnames(expr_avgse_df)[1] <- "GeneName"
    
    deps_df <- limmafnc_dt(type = "PROT", c_anno = c_anno, dat_gene = dat_gene, psm_count_table = psm_count_table, formule_contrast = formule_contrast,
                           expr_avgse_df = expr_avgse_df, signal_thr = signal_thr, fc_thr = fc_thr, pval_thr = pval_thr, pval_fdr = pval_fdr)
    
    deps_l_df <- deps_df$degs_l_df
    deps_w_df <- deps_df$degs_w_df
    
    dat_pep <- copy(proteome_data$dat_pep)
    psm_count_table <- data.table(table(proteome_data$dat_pep$ID_peptide))
    
    # Process peptides
    expr_mat <- copy(proteome_data$dat_pep)
    expr_l_df <- melt(expr_mat, id.vars = "ID_peptide", variable.name = "sample", value.name = "expr")
    expr_l_df <- expr_l_df[c_anno, on = "sample"]
    expr_cond_df <- expr_l_df[, .(N = .N, avg = mean(expr), sd = sd(expr), CV = sd(expr) / mean(expr), se = sd(expr) / sqrt(.N)), by = .(condition, ID_peptide)]
    
    expr_avg_pep_df <- dcast(expr_cond_df[, .(ID_peptide, condition, avg)], ID_peptide ~ condition, value.var = "avg")
    colnames(expr_avg_pep_df)[-1] <- paste0(colnames(expr_avg_pep_df)[-1], "_avg")
    expr_se_pep_df <- dcast(expr_cond_df[, .(ID_peptide, condition, se)], ID_peptide ~ condition, value.var = "se")
    colnames(expr_se_pep_df)[-1] <- paste0(colnames(expr_se_pep_df)[-1], "_se")
    expr_cv_pep_df <- dcast(expr_cond_df[, .(ID_peptide, condition, CV)], ID_peptide ~ condition, value.var = "CV")
    colnames(expr_cv_pep_df)[-1] <- paste0(colnames(expr_cv_pep_df)[-1], "_Coef_Variant_(%)")
    
    expr_avgse_pep_df <- Reduce(merge, list(expr_avg_pep_df, expr_se_pep_df, expr_cv_pep_df))
    
    deps_pep_df <- limmafnc_dt("PEP", c_anno = c_anno, dat_gene = dat_pep, psm_count_table = psm_count_table, formule_contrast = formule_contrast,
                               expr_avgse_df = expr_avgse_df, signal_thr = signal_thr, fc_thr = fc_thr, pval_thr = pval_thr, pval_fdr = pval_fdr)
    
    deps_pep_l_df <- deps_pep_df$degs_l_df
    deps_pep_w_df <- deps_pep_df$degs_w_df
    
    return(list(protein_results_long = deps_l_df, protein_results_wide = deps_w_df, 
                peptide_results_long = deps_pep_l_df, peptide_results_wide = deps_pep_w_df))
  }
  return(NULL)
}
