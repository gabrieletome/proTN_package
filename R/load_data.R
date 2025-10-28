#' Extract case study
#'
#' @param path_phospho Character; path where extract case study
#' @param path_proteome Character; path where extract case study
#' @param path_interactn Character; path where extract case study
#' @examples
#' \dontrun{
#' extract_example(path_proteome = tempdir())
#' }
#' @import utils
#' @import stringr
#' @importFrom grDevices cairo_pdf dev.off pdf
#' @importFrom stats cmdscale dist fisher.test model.matrix na.omit p.adjust prcomp rnorm sd setNames
#' @importFrom stats cutree median pchisq predict
#' @export
extract_example = function(path_phospho = NULL, 
                           path_proteome = NULL,
                           path_interactn = NULL) {
  if(!is.null(path_phospho)){
    message("Extracting phosphoproteome example...")
    dir<-path_phospho
    phospho_zip <- system.file("extdata", "case_study_phospho.zip", package = "proTN")
    unzip(phospho_zip, exdir = path_phospho)  
    message(paste0("Phospho-proteome example saved in: ",normalizePath(dir),"/case_study_phospho"))
  }
  if(!is.null(path_proteome)){
    message("Extracting proteome example...")
    dir<-path_proteome
    proteome_zip <- system.file("extdata", "case_study_proteomics.zip", package = "proTN")
    unzip(proteome_zip, exdir = path_proteome)
    message(paste0("Proteome example saved in: ",normalizePath(dir),"/case_study_proteomics"))
  }
  if(!is.null(path_interactn)){
    message("Extracting interactomics example...")
    dir<-path_interactn
    interactn_zip <- system.file("extdata", "case_study_interactn.zip", package = "proTN")
    unzip(interactn_zip, exdir = path_interactn)
    message(paste0("Interactomics example saved in: ",normalizePath(dir),"/case_study_interactn"))
  }
  if(is.null(path_phospho) & is.null(path_proteome) & is.null(path_interactn)){
    stop("Missing path. Provide at least one path")
  }
}

#' Read Proteomics Data
#'
#' This function reads proteomics data from files based on the specified software (Proteome Discoverer or MaxQuant). 
#' It supports reading annotation, peptide, and protein group data, with options for batch correction and filtering absent values.
#'
#' @param software Character; the proteomics software used, either "PD" for Proteome Discoverer or "MQ" for MaxQuant.
#' @param folder Character; the folder containing the data files.
#' @param use_proteinGroups_MQ Logical; Only for MaxQuant. If FALSE it expect the evidence.txt file, if TRUE require peptide.txt and proteinGroups.txt files. Default is FALSE-
#' @param annotation_filename Character; the name or pattern of the annotation file. Default is "annotation".
#' @param peptide_filename Character; the name or pattern of the peptide file. Default is "pep".
#' @param proteinGroup_filename Character; the name or pattern of the protein group file. Default is "prot".
#' @param condition_col Character; the column name representing condition information in the annotation file. Default is "Condition".
#' @param sample_col Character; the column name representing sample information in the annotation file. Default is "Sample".
#' @param color_col Character; the column name representing color information for the plot in the annotation file. Default is "Color".
#' @param batch_corr_exe Logical; whether batch correction should be applied. Default is FALSE.
#' @param batch_col Character; the column name representing batch information in the annotation file. Default is "batch".
#' @param filt_absent_value Numeric; the value used to filter out absent data. Default is 0.
#' @param min_peptide_protein Numeric; the value used to filter out protein by N peptide. Default is 0.
#'
#' @return A list containing proteomics data, including annotation and peptide data, 
#'         with batch correction applied if specified.
#'
#' @examples
#' \dontrun{
#' proteome_data <- read_proteomics(software = "MQ", folder = "data_folder")
#' }
#'
#' @importFrom dplyr ungroup mutate filter group_by n
#' @import data.table
#' @import readxl
#' @export
read_proteomics <- function(software, folder, use_proteinGroups_MQ = FALSE, annotation_filename = "annotation",
                            peptide_filename = "pep", proteinGroup_filename = "prot", 
                            condition_col="Condition", sample_col="Sample", 
                            color_col="Color", batch_corr_exe = FALSE, batch_col="batch",
                            filt_absent_value = 0,min_peptide_protein = 0){
  
  if(!(software %in% c("PD","MQ","SP","FP"))){
    stop("Valid software is required. Write PD or MQ or SP or FP.",
         "\tPD: Proteome Discoverer",
         "\tMQ: MaxQuant",
         "\tSP: Spectronaut",
         "\tFP: FragPipe")
  }
  
  anno_filename = list.files(folder, pattern = annotation_filename, full.names = T, ignore.case = T)
  if(is.null(anno_filename) | length(anno_filename) > 1){
    if(software != "SP"){
      stop("Missing file annotation or wrong annotation filename parameter or multiple files with same pattern")
    } else{
      # Update Column name for Spectronaut Annotation
      condition_col="R.Condition"
      sample_col="R.FileName"
      anno_filename = NULL
    }
  }
  pep_filename = list.files(folder, pattern = peptide_filename, full.names = T, ignore.case = T)
  if(is.null(pep_filename) | length(pep_filename) > 1){
    stop("Missing file peptide or wrong peptide filename parameter or multiple files with same pattern")
  }
  if(software == "PD"){
    prot_filename = list.files(folder, pattern = proteinGroup_filename, full.names = T, ignore.case = T)
    if(is.null(prot_filename) | length(prot_filename) > 1){
      stop("Missing file proteinGroup or wrong proteinGroup filename parameter or multiple files with same pattern")
    }
  } else if(use_proteinGroups_MQ & software != "SP" & software != "FP"){
    prot_filename = list.files(folder, pattern = proteinGroup_filename, full.names = T, ignore.case = T)
    if(is.null(prot_filename) | length(prot_filename) > 1){
      stop("Missing file proteinGroup or wrong proteinGroup filename parameter or multiple files with same pattern")
    }
  }else{
    prot_filename = NULL
  }
  
  proteome_data = NULL
  if(software == "PD"){
    message("Reading Proteome Discover files...")
    proteome_data = read_PD_files(anno_filename, pep_filename, prot_filename, 
                                  condition_col=condition_col, 
                                  sample_col=sample_col, 
                                  color_col=color_col,
                                  batch_corr_exe= batch_corr_exe, 
                                  filt_absent_value = filt_absent_value, 
                                  batch_col = batch_col,
                                  min_peptide_protein = min_peptide_protein)
  } else if(software == "MQ" & !use_proteinGroups_MQ){
    message("Reading MaxQuant files (evidence.txt)...")
    proteome_data = read_MQ_files(anno_filename, pep_filename, 
                                  condition_col=condition_col, 
                                  sample_col=sample_col, 
                                  color_col=color_col,
                                  batch_corr_exe = batch_corr_exe, 
                                  filt_absent_value = filt_absent_value, 
                                  batch_col = batch_col,
                                  min_peptide_protein = min_peptide_protein)
  } else if(software == "MQ" & use_proteinGroups_MQ){
    message("Reading MaxQuant files (peptides.txt & proteinGroups.txt)...")
    proteome_data = read_MQ_prot_peptide_files(anno_filename, pep_filename, prot_filename, 
                                               condition_col=condition_col, 
                                               sample_col=sample_col, 
                                               color_col=color_col,
                                               batch_corr_exe = batch_corr_exe, 
                                               filt_absent_value = filt_absent_value, 
                                               batch_col = batch_col,
                                               min_peptide_protein = min_peptide_protein)
  } else if(software == "SP"){
    message("Reading Spectronaut files...")
    proteome_data = read_Spectronaut_files(anno_filename,
                                           pep_filename, 
                                           condition_col=condition_col, 
                                           sample_col=sample_col, 
                                           color_col=color_col,
                                           batch_corr_exe = batch_corr_exe, 
                                           filt_absent_value = filt_absent_value, 
                                           batch_col = batch_col,
                                           min_peptide_protein = min_peptide_protein)
    
  } else if(software == "FP"){
    message("Reading FragPipe files...")
    proteome_data = read_FragPipe_files(anno_filename,
                                        pep_filename, 
                                        condition_col=condition_col, 
                                        sample_col=sample_col, 
                                        color_col=color_col,
                                        batch_corr_exe = batch_corr_exe, 
                                        filt_absent_value = filt_absent_value, 
                                        batch_col = batch_col,
                                        min_peptide_protein = min_peptide_protein)
    
  }
  
  return(proteome_data)
}


# Read MaxQuant files
read_MQ_files <- function(anno_filename, pep_filename, 
                          condition_col="Condition", sample_col="Sample", 
                          color_col="Color", batch_corr_exe = FALSE, batch_col="batch",
                          filt_absent_value = 0, min_peptide_protein = 0){
  input_files <- list()
  
  message("Reading files...")
  # Read annotation file
  input_files[["annotation"]] <- tryCatch({
    as.data.table(read_xlsx(anno_filename))
  }, error=function(cond){
    stop(paste0("Missing file. The file \'ANNOTATION\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  
  # Read peptide file
  input_files[["PEP"]] <- tryCatch({
    fread(pep_filename)
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PEPTIDE\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  message("File read.")
  
  message("Starting preprocessing...")
  #Clean files and merge
  colToKeep<-intersect(colnames(input_files[["annotation"]]), c(condition_col, sample_col, color_col, batch_col))
  if(!(condition_col %in% colToKeep)){
    stop("\'Condition\' column  in \'ANNOTATION\' file.")
  }
  if(!(sample_col %in% colToKeep)){
    stop("\'Sample\' column  in \'ANNOTATION\' file.")
  }
  if(batch_corr_exe & !(batch_col %in% colToKeep)){ 
    stop(paste0("\'",batch_col,"\' column  in \'ANNOTATION\' file with batch correction activated"))
  }
  input_files[["annotation"]] <- input_files[["annotation"]][, ..colToKeep]
  setnames(input_files[["annotation"]], old = condition_col, new = "Condition")
  setnames(input_files[["annotation"]], old = sample_col, new = "Sample")
  
  # Cast to string if all condition are numeric
  input_files[["annotation"]][, Condition := as.character(Condition)]
  # Check if conditions start with number. It generate problem later in DEqMS
  if(any(str_starts(input_files[["annotation"]]$Condition, "[0-9]"))){
    input_files[["annotation"]][(str_starts(Condition, "[0-9]")), Condition := str_c("X.", Condition)]
  }
  
  # Check if "Gene Names" and "Protein names" columnare present. Otherwise add them
  # TODO: check if true for all cases
  if(!("Gene names" %in% colnames(input_files[["PEP"]]))){
    input_files[["PEP"]][, `Gene names` := tstrsplit(tstrsplit(`Leading razor protein`, "\\|", keep = 3)[[1]], "_", keep = 1)[[1]]]
    input_files[["PEP"]][, `Leading razor protein` := tstrsplit(`Leading razor protein`, "\\|", keep = 2)[[1]]]
  }
  if(!("Protein names" %in% colnames(input_files[["PEP"]]))){
    data("anno_uniprot", envir = environment())
    anno_uniprot[, `Gene Names`:=NULL]
    if(length(intersect(anno_uniprot$`Leading razor protein`, input_files[["PEP"]]$`Leading razor protein`)) == 0){
      input_files[["PEP"]][, `Protein names` := `Leading razor protein`]
    }else{
      input_files[["PEP"]] <- merge.data.table(input_files[["PEP"]], unique(anno_uniprot, by="Leading razor protein", fromLast=FALSE), 
                                               by="Leading razor protein", 
                                               all.x = TRUE)
    }
  } 
  if(nrow(input_files[["PEP"]][is.na(`Protein names`)])>0){
    data("anno_uniprot", envir = environment())
    anno_uniprot[, `Gene Names`:=NULL]
    unique_anno_uniprot <- unique(anno_uniprot, by="Leading razor protein", fromLast=FALSE)
    setnames(unique_anno_uniprot, old=c("Protein names"), new = c("desc"))
    input_files[["PEP"]][, `Protein names` := as.character(`Protein names`)]
    input_files[["PEP"]][is.na(`Protein names`), `Protein names` := unique_anno_uniprot[.SD, on=.(`Leading razor protein`), x.desc]]
  }
  
  
  # Manage Peptide file
  message("Cleaning data...")
  initial_peptide = nrow(input_files[["PEP"]])
  message(paste0("\tRaw number of peptide: ",initial_peptide))
  
  to_remove <- nrow(input_files[["PEP"]][grepl("Keratin|keratin", `Protein names`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!grepl("Keratin|keratin", `Protein names`)]
  message(paste0("\tKeratin peptide removed: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][grepl("CON_",`Leading proteins`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!grepl("CON_",`Leading proteins`)]
  message(paste0("\tCONTAMINANT peptide removed: ",to_remove))
  
  # to_remove <- nrow(input_files[["PEP"]][is.na(`Protein names`)])
  # input_files[["PEP"]] <- input_files[["PEP"]][!is.na(`Protein names`)]
  # message(paste0("\tPeptide removed due to missing Protein Names: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][is.na(`Gene names`) | stri_isempty(`Gene names`) | `Gene names` == ""])
  input_files[["PEP"]] <- input_files[["PEP"]][!(is.na(`Gene names`) | stri_isempty(`Gene names`) | `Gene names` == "")]
  message(paste0("\tPeptide removed due to missing Gene Names: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][!(`Raw file` %in% input_files[["annotation"]]$Sample)])
  input_files[["PEP"]] <- input_files[["PEP"]][(`Raw file` %in% input_files[["annotation"]]$Sample)]
  input_files[["annotation"]] <- input_files[["annotation"]][Sample %in% input_files[["PEP"]]$`Raw file`]
  message(paste0("\tPeptide removed since sample not in sample annotation: ",to_remove))
  
  # to_remove <- nrow(input_files[["PEP"]][(Type == "MSMS")])
  # input_files[["PEP"]] <- input_files[["PEP"]][(Type != "MSMS")]
  # message(paste0("\tPeptide removed since sample not in sample annotation: ",to_remove))
  
  # Keep only first gene name
  input_files[["PEP"]][str_detect(`Gene names`,";"), `Gene names` := tstrsplit(`Gene names`, ";", keep = 1)]
  input_files[["PEP"]][, Modifications:=gsub(" ", "_", Modifications)]
  #Made the matrix
  psm_sig_raw <- data.table("ID_peptide" = as.factor(paste(input_files[["PEP"]]$`Gene names`, input_files[["PEP"]]$Sequence, input_files[["PEP"]]$Modifications, sep="_")), 
                            "Sample" = as.factor(input_files[["PEP"]]$`Raw file`), 
                            "Intensity" = input_files[["PEP"]]$Intensity)
  suppressWarnings({
    psm_sig_raw <- psm_sig_raw[, .(Intensity = sum(Intensity)), by = c("ID_peptide", "Sample")]
    psm_sig_raw <- dcast(psm_sig_raw, formula = ID_peptide~Sample, value.var = "Intensity")
  })
  colnames(psm_sig_raw)[-1] <- input_files[["annotation"]][match(colnames(psm_sig_raw)[-1], input_files[["annotation"]]$Sample)]$Sample
  
  #Made the description of the mpeptides
  psm_peptide_table <- as.data.table(unique(input_files[["PEP"]][, c("Leading razor protein", "Protein names", "Gene names", "Sequence", "Modifications")]))
  colnames(psm_peptide_table) <- c("Accession","Description","GeneName","Annotated_Sequence","Modifications")
  psm_peptide_table[, ID_peptide := paste(GeneName, Annotated_Sequence, Modifications, sep="_")]
  
  # Extract annotation
  c_anno<-input_files[["annotation"]]
  colnames(c_anno)<-tolower(colnames(c_anno))
  
  if(!("color" %in% colnames(c_anno))){
    message("Color column not found! Setting default color")
    c_anno<-merge.data.table(c_anno,
                             data.table("color"=colour_vec[1:length(unique(c_anno$condition))],
                                        "condition"=unique(c_anno$condition)),
                             by = "condition")
  }
  colour_vec<-na.omit(c_anno$color)
  names(colour_vec)<-na.omit(c_anno$sample)
  
  # Add chunk peptide annotation
  psm_anno_raw <- data.table("Accession" = psm_peptide_table$Accession,
                             "ID_peptide"=psm_peptide_table$ID_peptide,
                             "symbol"=psm_peptide_table$GeneName,
                             "sequence"=psm_peptide_table$Annotated_Sequence,
                             "modifications"=psm_peptide_table$Modifications)
  
  psm_anno_raw$old_id<-psm_anno_raw$ID_peptide
  
  # create peptide names linked to symbols
  psm_anno_raw[, row := seq(1,nrow(psm_anno_raw))]
  psm_anno_raw <- psm_anno_raw %>% group_by(symbol) %>% mutate("min"=min(row),"card"=n()) %>% ungroup()
  psm_anno_raw$rank<- psm_anno_raw$row - psm_anno_raw$min +1
  psm_anno_raw$id<- paste(psm_anno_raw$symbol,psm_anno_raw$rank,psm_anno_raw$card,sep="_")
  psm_anno_raw <- as.data.table(psm_anno_raw)
  
  # Convert to data.table
  psm_sig_prot_raw <- as.data.table(psm_sig_raw)
  psm_sig_pet_raw <- as.data.table(psm_sig_raw)
  
  # Preprocess Protein intensities
  psm_sig_prot_raw[psm_sig_prot_raw == 0] <- NA  # Transform 0s into NAs
  sig_thr <- filt_absent_value  # NA threshold
  
  psm_long_dt <- melt(psm_sig_prot_raw, id.vars = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_prot_df <- psm_sig_prot_raw[ID_peptide %in% filter_ID_peptide]
  psm_anno_df <- as.data.table(psm_anno_raw)[ID_peptide %in% filter_ID_peptide]
  
  # Log2 transformation
  psm_log_prot_df <- copy(psm_sig_prot_df)
  psm_log_prot_df[, (setdiff(names(psm_log_prot_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_prot_df), "ID_peptide")]
  
  # Filter proteins with only 1 peptID_peptidee
  psm_anno_df[, N_peptide := .N, by = symbol]
  filter_df_single_pep <- psm_anno_df[N_peptide >= min_peptide_protein, ID_peptide]
  psm_log_prot_df <- psm_log_prot_df[ID_peptide %in% filter_df_single_pep]
  psm_anno_df <- psm_anno_df[ID_peptide %in% filter_df_single_pep]
  
  # Preprocess peptID_peptidee intensities
  psm_sig_pet_raw[psm_sig_pet_raw == 0] <- NA  # Transform 0s into NAs
  
  psm_long_dt <- melt(psm_sig_pet_raw, id = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_pet_df <- psm_sig_pet_raw[ID_peptide %in% filter_ID_peptide]
  psm_peptide_table <- as.data.table(psm_peptide_table)[ID_peptide %in% filter_ID_peptide]
  
  # Determine tryptic condition
  peptides_df <- psm_peptide_table[, .(Accession, Annotated_Sequence)]
  peptides_df[, preAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,1], regex("\\[\\w+\\]")),
                                 2,
                                 str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                        n = 3)[,1], regex("\\[\\w+\\]")))-1)]
  peptides_df[, endAA := substr(Annotated_Sequence, nchar(Annotated_Sequence)-1, nchar(Annotated_Sequence)-1)]
  peptides_df[, postAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,3], regex("\\[\\w+\\]")),
                                  2,
                                  str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                         n = 3)[,3], regex("\\[\\w+\\]")))-1)]
  
  peptides_df[, fully_TRI := preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, NSEMI_TRI := preAA %in% c("K","R") & !endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, CSEMI_TRI := !preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, non_TRI := !fully_TRI & !NSEMI_TRI & !CSEMI_TRI]
  
  peptides_df[, tryptic_cond := fifelse(fully_TRI, "fully tryptic", fifelse(NSEMI_TRI, "N-semi tryptic", fifelse(CSEMI_TRI, "C-semi tryptic", "non tryptic")))]
  psm_peptide_table[, tryptic_cond := peptides_df$tryptic_cond]
  
  # Log2 transformation for peptID_peptidees
  psm_log_pet_df <- copy(psm_sig_pet_df)
  psm_log_pet_df[, (setdiff(names(psm_log_pet_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_pet_df), "ID_peptide")]
  
  psm_log_pet_df <- psm_log_pet_df[ID_peptide %in% filter_df_single_pep]
  psm_peptide_table <- psm_peptide_table[ID_peptide %in% filter_df_single_pep]
  
  message(paste0("N Proteins (after filter): ", uniqueN(psm_peptide_table$Accession)))
  message(paste0("N Peptides (after filter): ", uniqueN(psm_peptide_table$ID_peptide)))
  
  return(list("c_anno" = c_anno,
              "psm_anno_df" = psm_anno_df,
              "psm_log_prot_df" = psm_log_prot_df,
              "psm_peptide_table" = psm_peptide_table,
              "psm_log_pet_df" = psm_log_pet_df,
              "colour_vec" = colour_vec))
}

# Read Proteome Discoverer
read_PD_files <- function(anno_filename, pep_filename, prot_filename,
                          condition_col="Condition", sample_col="Sample", 
                          color_col="Color", batch_corr_exe = FALSE, batch_col="batch",
                          filt_absent_value = 0, min_peptide_protein=0){
  input_files <- list()
  
  message("Reading files...")
  # Read annotation file
  input_files[["annotation"]] <- tryCatch({
    as.data.table(read_xlsx(anno_filename))
  }, error=function(cond){
    stop(paste0("Missing file. The file \'ANNOTATION\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  
  # Read peptide file
  input_files[["PEP"]] <- tryCatch({
    as.data.table(read_xlsx(pep_filename))
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PEPTIDE\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  # Read protein file
  input_files[["PROT"]] <- tryCatch({
    as.data.table(read_xlsx(prot_filename))
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PROTEIN\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  message("File read.")
  
  message("Starting preprocessing...")
  #Clean files and merge
  colToKeep<-intersect(colnames(input_files[["annotation"]]), c("File ID",condition_col, sample_col, color_col, batch_col))
  if(!(condition_col %in% colToKeep)){
    stop("\'Condition\' column  in \'ANNOTATION\' file.")
  }
  if(!(sample_col %in% colToKeep)){
    stop("\'Sample\' column  in \'ANNOTATION\' file.")
  }
  if(batch_corr_exe & !(batch_col %in% colToKeep)){ 
    stop(paste0("\'",batch_col,"\' column  in \'ANNOTATION\' file with batch correction activated"))
  }
  input_files[["annotation"]] <- input_files[["annotation"]][, ..colToKeep]
  setnames(input_files[["annotation"]], old = "File ID", new = "File_ID")
  setnames(input_files[["annotation"]], old = condition_col, new = "Condition")
  setnames(input_files[["annotation"]], old = sample_col, new = "Sample")
  
  # Cast to string if all condition are numeric
  input_files[["annotation"]][, Condition := as.character(Condition)]
  # Check if conditions start with number. It generate problem later in DEqMS
  if(any(str_starts(input_files[["annotation"]]$Condition, "[0-9]"))){
    input_files[["annotation"]][(str_starts(Condition, "[0-9]")), Condition := str_c("X.", Condition)]
  }
  
  # Prepare PROT table
  message("Cleaning data...")
  initial_peptide = nrow(input_files[["PROT"]])
  message(paste0("\tRaw number of protein: ",initial_peptide))
  
  to_remove <- nrow(input_files[["PROT"]][grepl("Keratin|keratin", Description)])
  input_files[["PROT"]] <- input_files[["PROT"]][!grepl("Keratin|keratin", Description)]
  message(paste0("\tKeratin protein removed: ",to_remove))
  
  to_remove <- nrow(input_files[["PROT"]][grepl("CON_",Accession) | grepl("CON_",Description)])
  input_files[["PROT"]] <- input_files[["PROT"]][!(grepl("CON_",Accession) | grepl("CON_",Description))]
  message(paste0("\tCONTAMINANT protein removed: ",to_remove))
  
  input_files$PROT <- na.omit(input_files$PROT[, .(Accession, Description)])
  
  input_files$PROT[, GeneName := sub(".*GN=(.+) PE=.*", "\\1", Description)]
  
  to_remove <- nrow(input_files[["PROT"]][is.na(GeneName) | stri_isempty(GeneName) | GeneName == ""])
  input_files[["PROT"]] <- input_files[["PROT"]][!(is.na(GeneName) | stri_isempty(GeneName) | GeneName == "")]
  message(paste0("\tProtein removed due to missing Gene Names: ",to_remove))
  
  if(!("Modifications" %in% names(input_files$PEP))){
    setnames(input_files$PEP, grep("Modification",names(input_files$PEP), value = T), "Modifications")
  }
  
  # Maintain only the first UNIPROT code
  input_files$PEP[, `Master Protein Accessions` := sapply(`Master Protein Accessions`, function(x) strsplit(x, ";")[[1]][1])]
  
  # Verify and clean PEP table
  abundance_cols <- grep("Abundance", names(input_files$PEP), value = TRUE)
  input_files$PEP <- input_files$PEP[, c("Annotated Sequence", "Modifications", 
                                         "Master Protein Accessions", "Positions in Master Proteins", abundance_cols), with = FALSE]
  # Merge PROT and PEP tables
  input_files$PD_PEP_matrix <- merge.data.table(input_files$PROT, input_files$PEP, by.x = "Accession", by.y = "Master Protein Accessions")
  
  # Filter sample in file annotation
  colToKeep <- c(c("Accession","Description","GeneName","Annotated Sequence","Modifications",
                   "Positions in Master Proteins"),
                 lapply(input_files$annotation$File_ID, function(x){
                   colnames(input_files$PD_PEP_matrix)[grepl(str_c(x,":"), colnames(input_files$PD_PEP_matrix))]
                 }) %>% unlist())
  input_files$PD_PEP_matrix<-input_files$PD_PEP_matrix[, ..colToKeep]
  
  # TODO: check this code
  if("TRUE" %in% grepl("\\.", colnames(input_files$PD_PEP_matrix))){
    input_files$PD_PEP_matrix<-input_files$PD_PEP_matrix[,-which(colnames(input_files$PD_PEP_matrix) %in%
                                                                   c(colnames(input_files$PD_PEP_matrix)[grepl("\\.",
                                                                                                               colnames(input_files$PD_PEP_matrix))]))]
  }
  
  colnames(input_files$PD_PEP_matrix)<-c(c("Accession","Description","GeneName","Annotated_Sequence","Modifications","Position_in_Master_Proteins"),
                                         unlist(lapply(colnames(input_files$PD_PEP_matrix), function(x){
                                           lapply(input_files$annotation$File_ID, function(y){if(str_detect(x,paste0(y,":"))){
                                             input_files$annotation$Sample[which(input_files$annotation$File_ID %in% y)]
                                           }})
                                         })))
  
  
  input_files[["PD_PEP_matrix"]][, Modifications:=gsub(" ", "_", Modifications)]
  
  #Maintain only the new complete matrix of abundance
  input_files$PD_PEP_matrix[, ID_peptide := as.factor(paste(input_files[["PD_PEP_matrix"]]$GeneName, input_files[["PD_PEP_matrix"]]$Annotated_Sequence, input_files[["PD_PEP_matrix"]]$Modifications, sep="_"))]
  
  
  colToKeep <- c("ID_peptide",input_files$annotation$Sample)
  psm_sig_raw<-input_files$PD_PEP_matrix[, ..colToKeep]
  #Sum duplicated peptides
  psm_sig_raw<-psm_sig_raw[, lapply(.SD, sum), by = ID_peptide, .SDcols = is.numeric]
  
  #Made the description of the mpeptides
  psm_peptide_table <- as.data.table(unique(input_files[["PD_PEP_matrix"]][, c("ID_peptide","Accession","Description","GeneName","Annotated_Sequence","Modifications","Position_in_Master_Proteins")]))
  
  # Extract annotation
  c_anno<-input_files[["annotation"]]
  colnames(c_anno)<-tolower(colnames(c_anno))
  
  if(!("color" %in% colnames(c_anno))){
    message("Color column not found! Setting default color")
    c_anno<-merge.data.table(c_anno,
                             data.table("color"=colour_vec[1:length(unique(c_anno$condition))],
                                        "condition"=unique(c_anno$condition)),
                             by = "condition")
  }
  colour_vec<-na.omit(c_anno$color)
  names(colour_vec)<-na.omit(c_anno$sample)
  
  # Add chunk peptide annotation
  psm_anno_raw <- data.table("Accession" = psm_peptide_table$Accession,
                             "ID_peptide"=psm_peptide_table$ID_peptide,
                             "symbol"=psm_peptide_table$GeneName,
                             "sequence"=psm_peptide_table$Annotated_Sequence,
                             "modifications"=psm_peptide_table$Modifications)
  
  psm_anno_raw$old_id<-psm_anno_raw$ID_peptide
  
  # create peptide names linked to symbols
  psm_anno_raw[, row := seq(1,nrow(psm_anno_raw))]
  psm_anno_raw <- psm_anno_raw %>% group_by(symbol) %>% mutate("min"=min(row),"card"=n()) %>% ungroup()
  psm_anno_raw$rank<- psm_anno_raw$row - psm_anno_raw$min +1
  psm_anno_raw$id<- paste(psm_anno_raw$symbol,psm_anno_raw$rank,psm_anno_raw$card,sep="_")
  psm_anno_raw <- as.data.table(psm_anno_raw)
  
  
  # Convert to data.table
  psm_sig_prot_raw <- as.data.table(psm_sig_raw)
  psm_sig_pet_raw <- as.data.table(psm_sig_raw)
  
  # Preprocess Protein intensities
  psm_sig_prot_raw[psm_sig_prot_raw == 0] <- NA  # Transform 0s into NAs
  sig_thr <- filt_absent_value  # NA threshold
  
  psm_long_dt <- melt(psm_sig_prot_raw, id.vars = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_prot_df <- psm_sig_prot_raw[ID_peptide %in% filter_ID_peptide]
  psm_anno_df <- as.data.table(psm_anno_raw)[ID_peptide %in% filter_ID_peptide]
  
  # Log2 transformation
  psm_log_prot_df <- copy(psm_sig_prot_df)
  psm_log_prot_df[, (setdiff(names(psm_log_prot_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_prot_df), "ID_peptide")]
  
  # Filter proteins with only 1 peptID_peptidee
  psm_anno_df[, N_peptide := .N, by = symbol]
  filter_df_single_pep <- psm_anno_df[N_peptide >= min_peptide_protein, ID_peptide]
  psm_log_prot_df <- psm_log_prot_df[ID_peptide %in% filter_df_single_pep]
  psm_anno_df <- psm_anno_df[ID_peptide %in% filter_df_single_pep]
  
  # Preprocess peptID_peptidee intensities
  psm_sig_pet_raw[psm_sig_pet_raw == 0] <- NA  # Transform 0s into NAs
  
  psm_long_dt <- melt(psm_sig_pet_raw, id = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_pet_df <- psm_sig_pet_raw[ID_peptide %in% filter_ID_peptide]
  psm_peptide_table <- as.data.table(psm_peptide_table)[ID_peptide %in% filter_ID_peptide]
  
  # Determine tryptic condition
  # TODO: check this triptic test
  peptides_df <- psm_peptide_table[, .(Accession, Annotated_Sequence)]
  peptides_df[, preAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,1], regex("\\[\\w+\\]")),
                                 2,
                                 str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                        n = 3)[,1], regex("\\[\\w+\\]")))-1)]
  peptides_df[, endAA := substr(Annotated_Sequence, nchar(Annotated_Sequence)-1, nchar(Annotated_Sequence)-1)]
  peptides_df[, postAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,3], regex("\\[\\w+\\]")),
                                  2,
                                  str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                         n = 3)[,3], regex("\\[\\w+\\]")))-1)]
  
  peptides_df[, fully_TRI := preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, NSEMI_TRI := preAA %in% c("K","R") & !endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, CSEMI_TRI := !preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, non_TRI := !fully_TRI & !NSEMI_TRI & !CSEMI_TRI]
  
  peptides_df[, tryptic_cond := fifelse(fully_TRI, "fully tryptic", fifelse(NSEMI_TRI, "N-semi tryptic", fifelse(CSEMI_TRI, "C-semi tryptic", "non tryptic")))]
  psm_peptide_table[, tryptic_cond := peptides_df$tryptic_cond]
  
  # Log2 transformation for peptID_peptidees
  psm_log_pet_df <- copy(psm_sig_pet_df)
  psm_log_pet_df[, (setdiff(names(psm_log_pet_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_pet_df), "ID_peptide")]
  
  psm_log_pet_df <- psm_log_pet_df[ID_peptide %in% filter_df_single_pep]
  psm_peptide_table <- psm_peptide_table[ID_peptide %in% filter_df_single_pep]
  
  message(paste0("N Proteins (after filter): ", uniqueN(psm_peptide_table$Accession)))
  message(paste0("N Peptides (after filter): ", uniqueN(psm_peptide_table$ID_peptide)))
  
  return(list("c_anno" = c_anno,
              "psm_anno_df" = psm_anno_df,
              "psm_log_prot_df" = psm_log_prot_df,
              "psm_peptide_table" = psm_peptide_table,
              "psm_log_pet_df" = psm_log_pet_df,
              "colour_vec" = colour_vec))
}

# Read MaxQuant files
read_MQ_prot_peptide_files <- function(anno_filename, pep_filename, prot_filename,
                                       condition_col="Condition", sample_col="Sample", 
                                       color_col="Color", batch_corr_exe = FALSE, batch_col="batch",
                                       filt_absent_value = 0, min_peptide_protein = 0){
  input_files <- list()
  
  message("Reading files...")
  # Read annotation file
  input_files[["annotation"]] <- tryCatch({
    as.data.table(read_xlsx(anno_filename))
  }, error=function(cond){
    stop(paste0("Missing file. The file \'ANNOTATION\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  
  # Read peptide file
  input_files[["PEP"]] <- tryCatch({
    fread(pep_filename)
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PEPTIDE\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  
  # Read peptide file
  input_files[["PROT"]] <- tryCatch({
    fread(prot_filename)
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PROTEIN_GROUPS\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  message("File read.")
  message("Starting preprocessing...")
  #Clean files and merge
  colToKeep<-intersect(colnames(input_files[["annotation"]]), c(condition_col, sample_col, color_col, batch_col))
  if(!(condition_col %in% colToKeep)){
    stop("\'Condition\' column  in \'ANNOTATION\' file.")
  }
  if(!(sample_col %in% colToKeep)){
    stop("\'Sample\' column  in \'ANNOTATION\' file.")
  }
  if(batch_corr_exe & !(batch_col %in% colToKeep)){ 
    stop(paste0("\'",batch_col,"\' column  in \'ANNOTATION\' file with batch correction activated"))
  }
  input_files[["annotation"]] <- input_files[["annotation"]][, ..colToKeep]
  setnames(input_files[["annotation"]], old = condition_col, new = "Condition")
  setnames(input_files[["annotation"]], old = sample_col, new = "Sample")
  
  # Cast to string if all condition are numeric
  input_files[["annotation"]][, Condition := as.character(Condition)]
  # Check if conditions start with number. It generate problem later in DEqMS
  if(any(str_starts(input_files[["annotation"]]$Condition, "[0-9]"))){
    input_files[["annotation"]][(str_starts(Condition, "[0-9]")), Condition := str_c("X.", Condition)]
  }
  
  # Check if "Gene Names" and "Protein names" columnare present. Otherwise add them
  # TODO: check if true for all cases
  if(!("Gene names" %in% colnames(input_files[["PEP"]]))){
    input_files[["PEP"]][, `Gene names` := tstrsplit(tstrsplit(`Leading razor protein`, "\\|", keep = 3)[[1]], "_", keep = 1)[[1]]]
    input_files[["PEP"]][, `Leading razor protein` := tstrsplit(`Leading razor protein`, "\\|", keep = 2)[[1]]]
  }
  if(!("Protein names" %in% colnames(input_files[["PEP"]]))){
    data("anno_uniprot", envir = environment())
    anno_uniprot[, `Gene Names`:=NULL]
    if(length(intersect(anno_uniprot$`Leading razor protein`, input_files[["PEP"]]$`Leading razor protein`)) == 0){
      input_files[["PEP"]][, `Protein names` := `Leading razor protein`]
    }else{
      input_files[["PEP"]] <- merge.data.table(input_files[["PEP"]], unique(anno_uniprot, by="Leading razor protein", fromLast=FALSE), 
                                               by="Leading razor protein", 
                                               all.x = TRUE)
    }
  } 
  if(nrow(input_files[["PEP"]][is.na(`Protein names`)])>0){
    data("anno_uniprot", envir = environment())
    anno_uniprot[, `Gene Names`:=NULL]
    unique_anno_uniprot <- unique(anno_uniprot, by="Leading razor protein", fromLast=FALSE)
    setnames(unique_anno_uniprot, old=c("Protein names"), new = c("desc"))
    input_files[["PEP"]][, `Protein names` := as.character(`Protein names`)]
    input_files[["PEP"]][is.na(`Protein names`), `Protein names` := unique_anno_uniprot[.SD, on=.(`Leading razor protein`), x.desc]]
  }
  
  
  # Manage Peptide file
  message("Cleaning data...")
  initial_peptide = nrow(input_files[["PEP"]])
  message(paste0("\tRaw number of peptide: ",initial_peptide))
  
  to_remove <- nrow(input_files[["PEP"]][grepl("Keratin|keratin", `Protein names`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!grepl("Keratin|keratin", `Protein names`)]
  message(paste0("\tKeratin peptide removed: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][grepl("CON_",`Leading razor protein`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!grepl("CON_",`Leading razor protein`)]
  message(paste0("\tCONTAMINANT peptide removed: ",to_remove))
  
  # to_remove <- nrow(input_files[["PEP"]][is.na(`Protein names`)])
  # input_files[["PEP"]] <- input_files[["PEP"]][!is.na(`Protein names`)]
  # message(paste0("\tPeptide removed due to missing Protein Names: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][is.na(`Gene names`) | stri_isempty(`Gene names`) | `Gene names` == ""])
  input_files[["PEP"]] <- input_files[["PEP"]][!(is.na(`Gene names`) | stri_isempty(`Gene names`) | `Gene names` == "")]
  message(paste0("\tPeptide removed due to missing Gene Names: ",to_remove))
  
  # Keep only first gene name
  input_files[["PEP"]][str_detect(`Gene names`,";"), `Gene names` := tstrsplit(`Gene names`, ";", keep = 1)]
  input_files[["PEP"]][, Modifications:=NA]
  
  # Filter protein files
  input_files[["PROT"]] <- input_files[["PROT"]][!grepl("Keratin|keratin", `Fasta headers`) & !grepl("CON_|Keratin|keratin", `Majority protein IDs`), .(`Majority protein IDs`, `Fasta headers`)]
  input_files[["PROT"]][, `Leading razor protein` := tstrsplit(`Majority protein IDs`, "\\;", keep = 1)[[1]]]
  input_files[["PROT"]][, `Majority protein IDs` := NULL]
  if(any(grepl("\\|", input_files[["PROT"]]$`Leading razor protein`))){
    input_files[["PROT"]][, `Leading razor protein` := tstrsplit(`Leading razor protein`, "\\|", keep = 2)[[1]]]
  }
  
  input_files[["PROT_PEP"]] <- merge.data.table(input_files[["PROT"]], input_files[["PEP"]], 
                                                by.x = "Leading razor protein", by.y = "Leading razor protein")
  suppressWarnings({
    input_files[["PROT_PEP"]] <- melt(input_files[["PROT_PEP"]], 
                                      id.vars = c("Leading razor protein", "Protein names", "Fasta headers", "Gene names", "Sequence", "Modifications"),
                                      measure.vars = grep("Intensity ", names(input_files$PEP), value = TRUE), 
                                      variable.name = "Raw file", 
                                      value.name = "Intensity")
  })
  input_files[["PROT_PEP"]][, `Raw file` := stri_replace(`Raw file`, regex = "Intensity ", replacement = "")]
  
  to_remove <- nrow(input_files[["PROT_PEP"]][!(`Raw file` %in% input_files[["annotation"]]$Sample)])
  input_files[["PROT_PEP"]] <- input_files[["PROT_PEP"]][(`Raw file` %in% input_files[["annotation"]]$Sample)]
  input_files[["annotation"]] <- input_files[["annotation"]][Sample %in% input_files[["PROT_PEP"]]$`Raw file`]
  message(paste0("\tPeptide removed since sample not in sample annotation: ",to_remove))
  
  #Made the matrix
  psm_sig_raw <- data.table("ID_peptide" = as.factor(paste(input_files[["PROT_PEP"]]$`Gene names`, input_files[["PROT_PEP"]]$Sequence, input_files[["PEP"]]$Modifications, sep="_")), 
                            "Sample" = as.factor(input_files[["PROT_PEP"]]$`Raw file`), 
                            "Intensity" = input_files[["PROT_PEP"]]$Intensity)
  suppressWarnings({
    psm_sig_raw <- psm_sig_raw[, .(Intensity = sum(Intensity)), by = c("ID_peptide", "Sample")]
    psm_sig_raw <- dcast(psm_sig_raw, formula = ID_peptide~Sample, value.var = "Intensity")
  })
  colnames(psm_sig_raw)[-1] <- input_files[["annotation"]][match(colnames(psm_sig_raw)[-1], input_files[["annotation"]]$Sample)]$Sample
  
  #Made the description of the mpeptides
  psm_peptide_table <- as.data.table(unique(input_files[["PROT_PEP"]][, c("Leading razor protein", "Protein names", "Fasta headers", "Gene names", "Sequence", "Modifications")]))
  colnames(psm_peptide_table) <- c("Accession","Description","FASTA_description","GeneName","Annotated_Sequence","Modifications")
  psm_peptide_table[, ID_peptide := paste(GeneName, Annotated_Sequence, Modifications, sep="_")]
  
  # Extract annotation
  c_anno<-input_files[["annotation"]]
  colnames(c_anno)<-tolower(colnames(c_anno))
  
  if(!("color" %in% colnames(c_anno))){
    message("Color column not found! Setting default color")
    c_anno<-merge.data.table(c_anno,
                             data.table("color"=colour_vec[1:length(unique(c_anno$condition))],
                                        "condition"=unique(c_anno$condition)),
                             by = "condition")
  }
  colour_vec<-na.omit(c_anno$color)
  names(colour_vec)<-na.omit(c_anno$sample)
  
  # Add chunk peptide annotation
  psm_anno_raw <- data.table("Accession" = psm_peptide_table$Accession,
                             "ID_peptide"=psm_peptide_table$ID_peptide,
                             "symbol"=psm_peptide_table$GeneName,
                             "sequence"=psm_peptide_table$Annotated_Sequence,
                             "modifications"=psm_peptide_table$Modifications)
  
  psm_anno_raw$old_id<-psm_anno_raw$ID_peptide
  
  # create peptide names linked to symbols
  psm_anno_raw[, row := seq(1,nrow(psm_anno_raw))]
  psm_anno_raw <- psm_anno_raw %>% group_by(symbol) %>% mutate("min"=min(row),"card"=n()) %>% ungroup()
  psm_anno_raw$rank<- psm_anno_raw$row - psm_anno_raw$min +1
  psm_anno_raw$id<- paste(psm_anno_raw$symbol,psm_anno_raw$rank,psm_anno_raw$card,sep="_")
  psm_anno_raw <- as.data.table(psm_anno_raw)
  
  # Convert to data.table
  psm_sig_prot_raw <- as.data.table(psm_sig_raw)
  psm_sig_pet_raw <- as.data.table(psm_sig_raw)
  
  # Preprocess Protein intensities
  psm_sig_prot_raw[psm_sig_prot_raw == 0] <- NA  # Transform 0s into NAs
  sig_thr <- filt_absent_value  # NA threshold
  
  psm_long_dt <- melt(psm_sig_prot_raw, id.vars = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)

  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_prot_df <- psm_sig_prot_raw[ID_peptide %in% filter_ID_peptide]
  psm_anno_df <- as.data.table(psm_anno_raw)[ID_peptide %in% filter_ID_peptide]
  
  # Log2 transformation
  psm_log_prot_df <- copy(psm_sig_prot_df)
  psm_log_prot_df[, (setdiff(names(psm_log_prot_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_prot_df), "ID_peptide")]
  
  # Filter proteins with only 1 peptID_peptidee
  psm_anno_df[, N_peptide := .N, by = symbol]
  filter_df_single_pep <- psm_anno_df[N_peptide >= min_peptide_protein, ID_peptide]
  psm_log_prot_df <- psm_log_prot_df[ID_peptide %in% filter_df_single_pep]
  psm_anno_df <- psm_anno_df[ID_peptide %in% filter_df_single_pep]
  
  # Preprocess peptID_peptidee intensities
  psm_sig_pet_raw[psm_sig_pet_raw == 0] <- NA  # Transform 0s into NAs
  
  psm_long_dt <- melt(psm_sig_pet_raw, id = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_pet_df <- psm_sig_pet_raw[ID_peptide %in% filter_ID_peptide]
  psm_peptide_table <- as.data.table(psm_peptide_table)[ID_peptide %in% filter_ID_peptide]
  
  # Determine tryptic condition
  peptides_df <- psm_peptide_table[, .(Accession, Annotated_Sequence)]
  peptides_df[, preAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,1], regex("\\[\\w+\\]")),
                                 2,
                                 str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                        n = 3)[,1], regex("\\[\\w+\\]")))-1)]
  peptides_df[, endAA := substr(Annotated_Sequence, nchar(Annotated_Sequence)-1, nchar(Annotated_Sequence)-1)]
  peptides_df[, postAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,3], regex("\\[\\w+\\]")),
                                  2,
                                  str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                         n = 3)[,3], regex("\\[\\w+\\]")))-1)]
  
  peptides_df[, fully_TRI := preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, NSEMI_TRI := preAA %in% c("K","R") & !endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, CSEMI_TRI := !preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, non_TRI := !fully_TRI & !NSEMI_TRI & !CSEMI_TRI]
  
  peptides_df[, tryptic_cond := fifelse(fully_TRI, "fully tryptic", fifelse(NSEMI_TRI, "N-semi tryptic", fifelse(CSEMI_TRI, "C-semi tryptic", "non tryptic")))]
  psm_peptide_table[, tryptic_cond := peptides_df$tryptic_cond]
  
  # Log2 transformation for peptID_peptidees
  psm_log_pet_df <- copy(psm_sig_pet_df)
  psm_log_pet_df[, (setdiff(names(psm_log_pet_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_pet_df), "ID_peptide")]
  
  
  psm_log_pet_df <- psm_log_pet_df[ID_peptide %in% filter_df_single_pep]
  psm_peptide_table <- psm_peptide_table[ID_peptide %in% filter_df_single_pep]
  
  message(paste0("N Proteins (after filter): ", uniqueN(psm_peptide_table$Accession)))
  message(paste0("N Peptides (after filter): ", uniqueN(psm_peptide_table$ID_peptide)))
  
  return(list("c_anno" = c_anno,
              "psm_anno_df" = psm_anno_df,
              "psm_log_prot_df" = psm_log_prot_df,
              "psm_peptide_table" = psm_peptide_table,
              "psm_log_pet_df" = psm_log_pet_df,
              "colour_vec" = colour_vec))
}

# Read Spectronaut files
read_Spectronaut_files <- function(anno_filename, pep_filename, 
                                   condition_col="R.Condition", sample_col="R.FileName", 
                                   color_col="Color", batch_corr_exe = FALSE, batch_col="batch",
                                   filt_absent_value = 0, min_peptide_protein = 0){
  input_files <- list()
  
  message("Reading files...")
  # Read peptide file
  input_files[["PEP"]] <- tryCatch({
    fread(pep_filename)
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PEPTIDE\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  
  # Read annotation file
  if(!is.null(anno_filename)){
    message("Provide annotation excel. Use that as annotator of the sample.")
    input_files[["annotation"]] <- tryCatch({
      as.data.table(read_xlsx(anno_filename))
    }, error=function(cond){
      stop(paste0("Missing file. The file \'ANNOTATION\' is missing or not have the pattern in the filename or there are duplicates files."))
    })
  } else{
    message(paste0("Not provide annotation excel. Use column '", sample_col, "' for Sample and column '", condition_col, "' for the condition"))
  }
  
  message("File read.")
  
  message("Starting preprocessing...")
  
  # Read annotation file
  if(is.null(anno_filename)){
    colToKeep<-intersect(colnames(input_files[["PEP"]]), c(condition_col, sample_col, color_col, batch_col))
  } else{
    colToKeep<-intersect(colnames(input_files[["annotation"]]), c(condition_col, sample_col, color_col, batch_col))
  }
  #Clean files and merge
  if(!(condition_col %in% colToKeep)){
    stop("\'Condition\' column missing.")
  }
  if(!(sample_col %in% colToKeep)){
    stop("\'Sample\' column missing.")
  }
  if(batch_corr_exe & !(batch_col %in% colToKeep)){ 
    stop(paste0("\'",batch_col,"\' column missing with batch correction activated"))
  }
  
  if(is.null(anno_filename)){ input_files[["annotation"]] <- unique(input_files[["PEP"]][, ..colToKeep]) }
  setnames(input_files[["annotation"]], old = condition_col, new = "Condition")
  setnames(input_files[["annotation"]], old = sample_col, new = "Sample")
  
  # Cast to string if all condition are numeric
  input_files[["annotation"]][, Condition := as.character(Condition)]
  # Check if conditions start with number. It generate problem later in DEqMS
  if(any(str_starts(input_files[["annotation"]]$Condition, "[0-9]"))){
    input_files[["annotation"]][(str_starts(Condition, "[0-9]")), Condition := str_c("X.", Condition)]
  }
  
  
  # Remove Fragments columns (F.*)
  colToKeep_Spect <- grep("^F\\.", names(input_files[["PEP"]]), value = T, invert = T)
  input_files[["PEP"]] <- unique(input_files[["PEP"]][, ..colToKeep_Spect])
  # Add Accession and GeneName
  input_files[["PEP"]][, Accession := tstrsplit(PG.ProteinAccessions, "\\|", keep = 2)[[1]]]
  data("anno_uniprot", envir = environment())
  input_files[["PEP"]] <- merge.data.table(input_files[["PEP"]], anno_uniprot, by.x = "Accession", by.y = "Leading razor protein")
  
  
  # Manage Peptide file
  message("Cleaning data...")
  initial_peptide = nrow(input_files[["PEP"]])
  message(paste0("\tRaw number of peptide: ",initial_peptide))
  
  to_remove <- nrow(input_files[["PEP"]][grepl("Keratin|keratin", `Protein names`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!grepl("Keratin|keratin", `Protein names`)]
  message(paste0("\tKeratin peptide removed: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][grepl("CON_",`PG.ProteinAccessions`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!grepl("CON_",`PG.ProteinAccessions`)]
  message(paste0("\tCONTAMINANT peptide removed: ",to_remove))
  
  # to_remove <- nrow(input_files[["PEP"]][is.na(`Protein names`)])
  # input_files[["PEP"]] <- input_files[["PEP"]][!is.na(`Protein names`)]
  # message(paste0("\tPeptide removed due to missing Protein Names: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][is.na(`Gene Names`) | stri_isempty(`Gene Names`) | `Gene Names` == ""])
  input_files[["PEP"]] <- input_files[["PEP"]][!(is.na(`Gene Names`) | stri_isempty(`Gene Names`) | `Gene Names` == "")]
  message(paste0("\tPeptide removed due to missing Gene Names: ",to_remove))
  
  # Keep only first gene name
  input_files[["PEP"]][str_detect(`Gene Names`,";"), `Gene Names` := tstrsplit(`Gene Names`, ";", keep = 1)]
  input_files[["PEP"]][, EG.ModifiedSequence:=gsub(" ", "_", EG.ModifiedSequence)]
  #Made the matrix
  psm_sig_raw <- data.table("ID_peptide" = as.factor(paste(input_files[["PEP"]]$`Gene Names`, input_files[["PEP"]]$PEP.StrippedSequence, input_files[["PEP"]]$EG.ModifiedSequence, sep="_")), 
                            "Sample" = as.factor(input_files[["PEP"]][, ..sample_col][[1]]), 
                            "Intensity" = input_files[["PEP"]]$PEP.Quantity)
  suppressWarnings({
    psm_sig_raw <- psm_sig_raw[, .(Intensity = sum(Intensity)), by = c("ID_peptide", "Sample")]
    psm_sig_raw <- dcast(psm_sig_raw, formula = ID_peptide~Sample, value.var = "Intensity")
  })
  
  #Made the description of the mpeptides
  psm_peptide_table <- as.data.table(unique(input_files[["PEP"]][, c("Accession", "Protein names", "Gene Names", "PEP.StrippedSequence", "EG.ModifiedSequence")]))
  colnames(psm_peptide_table) <- c("Accession","Description","GeneName","Annotated_Sequence","Modifications")
  psm_peptide_table[, ID_peptide := paste(GeneName, Annotated_Sequence, Modifications, sep="_")]
  
  # Extract annotation
  c_anno<-input_files[["annotation"]]
  colnames(c_anno)<-tolower(colnames(c_anno))
  
  if(!("color" %in% colnames(c_anno))){
    message("Color column not found! Setting default color")
    c_anno<-merge.data.table(c_anno,
                             data.table("color"=colour_vec[1:length(unique(c_anno$condition))],
                                        "condition"=unique(c_anno$condition)),
                             by = "condition")
  }
  colour_vec<-na.omit(c_anno$color)
  names(colour_vec)<-na.omit(c_anno$sample)
  
  # Add chunk peptide annotation
  psm_anno_raw <- data.table("Accession" = psm_peptide_table$Accession,
                             "ID_peptide"=psm_peptide_table$ID_peptide,
                             "symbol"=psm_peptide_table$GeneName,
                             "sequence"=psm_peptide_table$Annotated_Sequence,
                             "modifications"=psm_peptide_table$Modifications)
  
  psm_anno_raw$old_id<-psm_anno_raw$ID_peptide
  
  # create peptide names linked to symbols
  psm_anno_raw[, row := seq(1,nrow(psm_anno_raw))]
  psm_anno_raw <- psm_anno_raw %>% group_by(symbol) %>% mutate("min"=min(row),"card"=n()) %>% ungroup()
  psm_anno_raw$rank<- psm_anno_raw$row - psm_anno_raw$min +1
  psm_anno_raw$id<- paste(psm_anno_raw$symbol,psm_anno_raw$rank,psm_anno_raw$card,sep="_")
  psm_anno_raw <- as.data.table(psm_anno_raw)
  
  
  # Convert to data.table
  psm_sig_prot_raw <- as.data.table(psm_sig_raw)
  psm_sig_pet_raw <- as.data.table(psm_sig_raw)
  
  # Preprocess Protein intensities
  psm_sig_prot_raw[psm_sig_prot_raw == 0] <- NA  # Transform 0s into NAs
  sig_thr <- filt_absent_value  # NA threshold
  
  psm_long_dt <- melt(psm_sig_prot_raw, id.vars = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_prot_df <- psm_sig_prot_raw[ID_peptide %in% filter_ID_peptide]
  psm_anno_df <- as.data.table(psm_anno_raw)[ID_peptide %in% filter_ID_peptide]
  
  # Log2 transformation
  psm_log_prot_df <- copy(psm_sig_prot_df)
  psm_log_prot_df[, (setdiff(names(psm_log_prot_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_prot_df), "ID_peptide")]
  
  # Filter proteins with only 1 peptID_peptidee
  psm_anno_df[, N_peptide := .N, by = symbol]
  filter_df_single_pep <- psm_anno_df[N_peptide >= min_peptide_protein, ID_peptide]
  psm_log_prot_df <- psm_log_prot_df[ID_peptide %in% filter_df_single_pep]
  psm_anno_df <- psm_anno_df[ID_peptide %in% filter_df_single_pep]
  
  # Preprocess peptID_peptidee intensities
  psm_sig_pet_raw[psm_sig_pet_raw == 0] <- NA  # Transform 0s into NAs
  
  psm_long_dt <- melt(psm_sig_pet_raw, id = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_pet_df <- psm_sig_pet_raw[ID_peptide %in% filter_ID_peptide]
  psm_peptide_table <- as.data.table(psm_peptide_table)[ID_peptide %in% filter_ID_peptide]
  
  # Determine tryptic condition
  peptides_df <- psm_peptide_table[, .(Accession, Annotated_Sequence)]
  peptides_df[, preAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,1], regex("\\[\\w+\\]")),
                                 2,
                                 str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                        n = 3)[,1], regex("\\[\\w+\\]")))-1)]
  peptides_df[, endAA := substr(Annotated_Sequence, nchar(Annotated_Sequence)-1, nchar(Annotated_Sequence)-1)]
  peptides_df[, postAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,3], regex("\\[\\w+\\]")),
                                  2,
                                  str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                         n = 3)[,3], regex("\\[\\w+\\]")))-1)]
  
  peptides_df[, fully_TRI := preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, NSEMI_TRI := preAA %in% c("K","R") & !endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, CSEMI_TRI := !preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, non_TRI := !fully_TRI & !NSEMI_TRI & !CSEMI_TRI]
  
  peptides_df[, tryptic_cond := fifelse(fully_TRI, "fully tryptic", fifelse(NSEMI_TRI, "N-semi tryptic", fifelse(CSEMI_TRI, "C-semi tryptic", "non tryptic")))]
  psm_peptide_table[, tryptic_cond := peptides_df$tryptic_cond]
  
  # Log2 transformation for peptID_peptidees
  psm_log_pet_df <- copy(psm_sig_pet_df)
  psm_log_pet_df[, (setdiff(names(psm_log_pet_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_pet_df), "ID_peptide")]
  
  
  psm_log_pet_df <- psm_log_pet_df[ID_peptide %in% filter_df_single_pep]
  psm_peptide_table <- psm_peptide_table[ID_peptide %in% filter_df_single_pep]
  
  message(paste0("N Proteins (after filter): ", uniqueN(psm_peptide_table$Accession)))
  message(paste0("N Peptides (after filter): ", uniqueN(psm_peptide_table$ID_peptide)))
  
  return(list("c_anno" = c_anno,
              "psm_anno_df" = psm_anno_df,
              "psm_log_prot_df" = psm_log_prot_df,
              "psm_peptide_table" = psm_peptide_table,
              "psm_log_pet_df" = psm_log_pet_df,
              "colour_vec" = colour_vec))
}


# Read FragPipe files
read_FragPipe_files <- function(anno_filename, pep_filename, 
                                condition_col="Condition", sample_col="Sample", 
                                color_col="Color", batch_corr_exe = FALSE, batch_col="batch",
                                filt_absent_value = 0, min_peptide_protein = 0){
  input_files <- list()
  
  message("Reading files...")
  # Read annotation file
  input_files[["annotation"]] <- tryCatch({
    as.data.table(read_xlsx(anno_filename))
  }, error=function(cond){
    stop(paste0("Missing file. The file \'ANNOTATION\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  
  # Read peptide file
  input_files[["PEP"]] <- tryCatch({
    fread(pep_filename)
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PEPTIDE\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  message("File read.")
  
  message("Starting preprocessing...")
  #Clean files and merge
  colToKeep<-intersect(colnames(input_files[["annotation"]]), c(condition_col, sample_col, color_col, batch_col))
  if(!(condition_col %in% colToKeep)){
    stop("\'Condition\' column  in \'ANNOTATION\' file.")
  }
  if(!(sample_col %in% colToKeep)){
    stop("\'Sample\' column  in \'ANNOTATION\' file.")
  }
  if(batch_corr_exe & !(batch_col %in% colToKeep)){ 
    stop(paste0("\'",batch_col,"\' column  in \'ANNOTATION\' file with batch correction activated"))
  }
  input_files[["annotation"]] <- input_files[["annotation"]][, ..colToKeep]
  setnames(input_files[["annotation"]], old = condition_col, new = "Condition")
  setnames(input_files[["annotation"]], old = sample_col, new = "Sample")
  
  # Cast to string if all condition are numeric
  input_files[["annotation"]][, Condition := as.character(Condition)]
  # Check if conditions start with number. It generate problem later in DEqMS
  if(any(str_starts(input_files[["annotation"]]$Condition, "[0-9]"))){
    input_files[["annotation"]][(str_starts(Condition, "[0-9]")), Condition := str_c("X.", Condition)]
  }
  
  
  # Manage Peptide file
  message("Cleaning data...")
  initial_peptide = nrow(input_files[["PEP"]])
  message(paste0("\tRaw number of peptide: ",initial_peptide))
  
  to_remove <- nrow(input_files[["PEP"]][grepl("Keratin|keratin", `Protein Description`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!grepl("Keratin|keratin", `Protein Description`)]
  message(paste0("\tKeratin peptide removed: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][grepl("CON_",`Protein ID`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!grepl("CON_",`Protein ID`)]
  message(paste0("\tCONTAMINANT peptide removed: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][is.na(`Gene`) | stri_isempty(`Gene`) | `Gene` == ""])
  input_files[["PEP"]] <- input_files[["PEP"]][!(is.na(`Gene`) | stri_isempty(`Gene`) | `Gene` == "")]
  message(paste0("\tPeptide removed due to missing Gene Names: ",to_remove))
  
  # Keep only first gene name
  input_files[["PEP"]][str_detect(`Gene`,";"), `Gene` := tstrsplit(`Gene`, ";", keep = 1)]
  input_files[["PEP"]][, `Assigned Modifications`:=gsub(" ", "_", `Assigned Modifications`)]
  
  peptides_df <- input_files[["PEP"]][, c("Gene", "Peptide Sequence", "Assigned Modifications", "Prev AA", "Next AA")]
  
  suppressWarnings({
    input_files[["PEP"]] <- melt(input_files[["PEP"]], 
                                 id.vars = c("Protein ID", "Protein Description", "Gene", "Peptide Sequence", "Assigned Modifications"),
                                 measure.vars = paste0(input_files$annotation$Sample, " Intensity"), 
                                 variable.name = "Sample", 
                                 value.name = "Intensity")
  })
  input_files[["PEP"]][, `Sample` := stri_replace(`Sample`, regex = " Intensity$", replacement = "")]
  
  to_remove <- nrow(input_files[["PEP"]][!(`Sample` %in% input_files[["annotation"]]$Sample)])
  input_files[["PEP"]] <- input_files[["PEP"]][(`Sample` %in% input_files[["annotation"]]$Sample)]
  input_files[["annotation"]] <- input_files[["annotation"]][Sample %in% input_files[["PEP"]]$`Sample`]
  message(paste0("\tPeptide removed since sample not in sample annotation: ",to_remove))
  
  #Made the matrix
  psm_sig_raw <- data.table("ID_peptide" = as.factor(paste(input_files[["PEP"]]$`Gene`, input_files[["PEP"]]$`Peptide Sequence`, input_files[["PEP"]]$`Assigned Modifications`, sep="_")), 
                            "Sample" = as.factor(input_files[["PEP"]]$`Sample`), 
                            "Intensity" = input_files[["PEP"]]$Intensity)
  suppressWarnings({
    psm_sig_raw <- psm_sig_raw[, .(Intensity = sum(Intensity)), by = c("ID_peptide", "Sample")]
    psm_sig_raw <- dcast(psm_sig_raw, formula = ID_peptide~Sample, value.var = "Intensity")
  })
  
  #Made the description of the mpeptides
  psm_peptide_table <- as.data.table(unique(input_files[["PEP"]][, c("Protein ID", "Protein Description", "Gene", "Peptide Sequence", "Assigned Modifications")]))
  colnames(psm_peptide_table) <- c("Accession","Description","GeneName","Annotated_Sequence","Modifications")
  psm_peptide_table[, ID_peptide := paste(GeneName, Annotated_Sequence, Modifications, sep="_")]
  
  # Extract annotation
  c_anno<-input_files[["annotation"]]
  colnames(c_anno)<-tolower(colnames(c_anno))
  
  if(!("color" %in% colnames(c_anno))){
    message("Color column not found! Setting default color")
    c_anno<-merge.data.table(c_anno,
                             data.table("color"=colour_vec[1:length(unique(c_anno$condition))],
                                        "condition"=unique(c_anno$condition)),
                             by = "condition")
  }
  colour_vec<-na.omit(c_anno$color)
  names(colour_vec)<-na.omit(c_anno$sample)
  
  # Add chunk peptide annotation
  psm_anno_raw <- data.table("Accession" = psm_peptide_table$Accession,
                             "ID_peptide"=psm_peptide_table$ID_peptide,
                             "symbol"=psm_peptide_table$GeneName,
                             "sequence"=psm_peptide_table$Annotated_Sequence,
                             "modifications"=psm_peptide_table$Modifications)
  
  psm_anno_raw$old_id<-psm_anno_raw$ID_peptide
  
  # create peptide names linked to symbols
  psm_anno_raw[, row := seq(1,nrow(psm_anno_raw))]
  psm_anno_raw <- psm_anno_raw %>% group_by(symbol) %>% mutate("min"=min(row),"card"=n()) %>% ungroup()
  psm_anno_raw$rank<- psm_anno_raw$row - psm_anno_raw$min +1
  psm_anno_raw$id<- paste(psm_anno_raw$symbol,psm_anno_raw$rank,psm_anno_raw$card,sep="_")
  psm_anno_raw <- as.data.table(psm_anno_raw)
  
  
  # Convert to data.table
  psm_sig_prot_raw <- as.data.table(psm_sig_raw)
  psm_sig_pet_raw <- as.data.table(psm_sig_raw)
  
  # Preprocess Protein intensities
  psm_sig_prot_raw[psm_sig_prot_raw == 0] <- NA  # Transform 0s into NAs
  sig_thr <- filt_absent_value  # NA threshold
  
  psm_long_dt <- melt(psm_sig_prot_raw, id.vars = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_prot_df <- psm_sig_prot_raw[ID_peptide %in% filter_ID_peptide]
  psm_anno_df <- as.data.table(psm_anno_raw)[ID_peptide %in% filter_ID_peptide]
  
  # Log2 transformation
  psm_log_prot_df <- copy(psm_sig_prot_df)
  psm_log_prot_df[, (setdiff(names(psm_log_prot_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_prot_df), "ID_peptide")]
  
  # Filter proteins with only 1 peptID_peptidee
  psm_anno_df[, N_peptide := .N, by = symbol]
  filter_df_single_pep <- psm_anno_df[N_peptide >= min_peptide_protein, ID_peptide]
  psm_log_prot_df <- psm_log_prot_df[ID_peptide %in% filter_df_single_pep]
  psm_anno_df <- psm_anno_df[ID_peptide %in% filter_df_single_pep]
  
  # Preprocess peptID_peptidee intensities
  psm_sig_pet_raw[psm_sig_pet_raw == 0] <- NA  # Transform 0s into NAs
  
  psm_long_dt <- melt(psm_sig_pet_raw, id = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_pet_df <- psm_sig_pet_raw[ID_peptide %in% filter_ID_peptide]
  psm_peptide_table <- as.data.table(psm_peptide_table)[ID_peptide %in% filter_ID_peptide]
  
  # Determine tryptic condition
  peptides_df[, ID_peptide := paste(`Gene`, `Peptide Sequence`, `Assigned Modifications`, sep="_")]
  peptides_df[, endAA := substr(`Peptide Sequence`, nchar(`Peptide Sequence`)-1, nchar(`Peptide Sequence`)-1)]
  
  peptides_df[, fully_TRI := `Prev AA` %in% c("K","R") & endAA %in% c("K","R") & (!`Next AA` %in% "P" | is.na(`Next AA`))]
  peptides_df[, NSEMI_TRI := `Prev AA` %in% c("K","R") & !endAA %in% c("K","R") & (!`Next AA` %in% "P" | is.na(`Next AA`))]
  peptides_df[, CSEMI_TRI := !`Prev AA` %in% c("K","R") & endAA %in% c("K","R") & (!`Next AA` %in% "P" | is.na(`Next AA`))]
  peptides_df[, non_TRI := !fully_TRI & !NSEMI_TRI & !CSEMI_TRI]
  
  peptides_df[, tryptic_cond := fifelse(fully_TRI, "fully tryptic", fifelse(NSEMI_TRI, "N-semi tryptic", fifelse(CSEMI_TRI, "C-semi tryptic", "non tryptic")))]
  peptides_df <- peptides_df[, c("ID_peptide", "tryptic_cond")]
  
  psm_peptide_table <- merge.data.table(psm_peptide_table, peptides_df, by = "ID_peptide")
  
  # Log2 transformation for peptID_peptidees
  psm_log_pet_df <- copy(psm_sig_pet_df)
  psm_log_pet_df[, (setdiff(names(psm_log_pet_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_pet_df), "ID_peptide")]
  
  
  psm_log_pet_df <- psm_log_pet_df[ID_peptide %in% filter_df_single_pep]
  psm_peptide_table <- psm_peptide_table[ID_peptide %in% filter_df_single_pep]
  
  message(paste0("N Proteins (after filter): ", uniqueN(psm_peptide_table$Accession)))
  message(paste0("N Peptides (after filter): ", uniqueN(psm_peptide_table$ID_peptide)))
  
  return(list("c_anno" = c_anno,
              "psm_anno_df" = psm_anno_df,
              "psm_log_prot_df" = psm_log_prot_df,
              "psm_peptide_table" = psm_peptide_table,
              "psm_log_pet_df" = psm_log_pet_df,
              "colour_vec" = colour_vec))
}


#' Read Phosphoproteomics Data
#'
#' This function reads phosphoproteomics data from files based on the specified software (Proteome Discoverer or MaxQuant). 
#' It supports reading annotation, peptide, protein group, and PSM data, with options for batch correction, filtering absent values, 
#' and keeping only phospho-modified peptides above a specified threshold.
#'
#' @param software Character; the proteomics software used, either "PD" for Proteome Discoverer or "MQ" for MaxQuant.
#' @param folder Character; the folder containing the data files.
#' @param keep_only_phosphomodification Logical; whether to keep only phospho-modified peptides. Default is TRUE.
#' @param phospho_thr Numeric; the threshold for phospho-modified peptides. Default is 0.75.
#' @param annotation_filename Character; the name or pattern of the annotation file. Default is "annotation".
#' @param peptide_filename Character; the name or pattern of the peptide file. Default is "pep".
#' @param proteinGroup_filename Character; the name or pattern of the protein group file. Default is "prot".
#' @param psm_filename Character; the name or pattern of the PSM file. Default is "psm".
#' @param condition_col Character; the column name representing condition information in the annotation file. Default is "Condition".
#' @param sample_col Character; the column name representing sample information in the annotation file. Default is "Sample".
#' @param color_col Character; the column name representing color information for the plot in the annotation file. Default is "Color".
#' @param batch_corr_exe Logical; whether batch correction should be applied. Default is FALSE.
#' @param batch_col Character; the column name representing batch information in the annotation file. Default is "batch".
#' @param filt_absent_value Numeric; the value used to filter out absent data. Default is 0.
#' @param min_peptide_protein Numeric; the value used to filter out protein by N peptide. Default is 0.
#'
#' @return A list containing phosphoproteomics data, including annotation and peptide data, 
#'         with batch correction applied if specified and only phospho-modified peptides kept if requested.
#'
#' @examples
#' \dontrun{
#' proteome_data <- read_phosphoproteomics(software = "MQ", folder = "data_folder")
#' }
#'
#' @importFrom dplyr ungroup mutate filter group_by n
#' @import data.table
#' @import readxl
#' @export
read_phosphoproteomics <- function(software, folder, 
                                   keep_only_phosphomodification = T, phospho_thr = 0.75, 
                                   annotation_filename = "annotation",peptide_filename = "pep", 
                                   proteinGroup_filename = "prot", psm_filename = "psm",
                                   condition_col="Condition", sample_col="Sample", 
                                   color_col="Color", batch_corr_exe = FALSE, batch_col="batch",
                                   filt_absent_value = 0, min_peptide_protein=0){
  
  if(!(software %in% c("PD","MQ"))){
    stop("Valid software is required. Write PD or MQ.",
         "\tPD: Proteome Discoverer",
         "\tMQ: MaxQuant")
  }
  
  # TODO: controllo che ci siano i file
  # TODO: batch correction verificare input
  
  anno_filename = list.files(folder, pattern = annotation_filename, full.names = T, ignore.case = T)
  if(is.null(anno_filename) | length(anno_filename) > 1){
    stop("Missing file annotation or wrong annotation filename parameter or multiple files with same pattern")
  }
  pep_filename = list.files(folder, pattern = peptide_filename, full.names = T, ignore.case = T)
  if(is.null(pep_filename) | length(pep_filename) > 1){
    stop("Missing file peptide or wrong peptide filename parameter or multiple files with same pattern")
  }
  if(software == "PD"){
    prot_filename = list.files(folder, pattern = proteinGroup_filename, full.names = T, ignore.case = T)
    if(is.null(prot_filename) | length(prot_filename) > 1){
      stop("Missing file proteinGroup or wrong proteinGroup filename parameter or multiple files with same pattern")
    }
  } else{
    prot_filename = NULL
  }
  if(software == "PD"){
    psm_file_filename = list.files(folder, pattern = psm_filename, full.names = T, ignore.case = T)
    if(is.null(psm_file_filename) | length(psm_file_filename) > 1){
      stop("Missing file PSM or wrong PSM filename parameter or multiple files with same pattern")
    }
  } else{
    psm_file_filename = NULL
  }
  
  proteome_data = NULL
  if(software == "PD"){
    proteome_data = read_phospho_PD_files(anno_filename, pep_filename, prot_filename, psm_file_filename, 
                                          keep_only_phosphomodification = keep_only_phosphomodification,
                                          batch_corr_exe= batch_corr_exe, 
                                          condition_col=condition_col, 
                                          sample_col=sample_col, 
                                          color_col=color_col,
                                          filt_absent_value = filt_absent_value, 
                                          phospho_thr = phospho_thr, 
                                          batch_col = batch_col,
                                          min_peptide_protein = min_peptide_protein)
  } else if(software == "MQ"){
    proteome_data = read_phospho_MQ_files(anno_filename, pep_filename, 
                                          keep_only_phosphomodification = keep_only_phosphomodification,
                                          condition_col=condition_col, 
                                          sample_col=sample_col, 
                                          color_col=color_col,
                                          batch_corr_exe = batch_corr_exe, 
                                          filt_absent_value = filt_absent_value, 
                                          phospho_thr = phospho_thr, 
                                          batch_col = batch_col,
                                          min_peptide_protein = min_peptide_protein)
  }
  
  return(proteome_data)
}


# Read MaxQuant files
read_phospho_MQ_files <- function(anno_filename, pep_filename, keep_only_phosphomodification = T,
                                  phospho_thr = 0.75, condition_col="Condition", sample_col="Sample", 
                                  color_col="Color", batch_corr_exe = FALSE, batch_col="batch",
                                  filt_absent_value = 0, min_peptide_protein = 0){
  input_files <- list()
  
  message("Reading files...")
  # Read annotation file
  input_files[["annotation"]] <- tryCatch({
    as.data.table(read_xlsx(anno_filename))
  }, error=function(cond){
    stop(paste0("Missing file. The file \'ANNOTATION\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  
  # Read peptide file
  input_files[["PEP"]] <- tryCatch({
    fread(pep_filename)
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PEPTIDE\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  message("File read.")
  
  message("Starting preprocessing...")
  #Clean files and merge
  colToKeep<-intersect(colnames(input_files[["annotation"]]), c(condition_col, sample_col, color_col, batch_col))
  if(!(condition_col %in% colToKeep)){
    stop("\'Condition\' column  in \'ANNOTATION\' file.")
  }
  if(!(sample_col %in% colToKeep)){
    stop("\'Sample\' column  in \'ANNOTATION\' file.")
  }
  if(batch_corr_exe & !(batch_col %in% colToKeep)){ 
    stop(paste0("\'",batch_col,"\' column  in \'ANNOTATION\' file with batch correction activated"))
  }
  input_files[["annotation"]] <- input_files[["annotation"]][, ..colToKeep]
  setnames(input_files[["annotation"]], old = condition_col, new = "Condition")
  setnames(input_files[["annotation"]], old = sample_col, new = "Sample")
  
  # Cast to string if all condition are numeric
  input_files[["annotation"]][, Condition := as.character(Condition)]
  # Check if conditions start with number. It generate problem later in DEqMS
  if(any(str_starts(input_files[["annotation"]]$Condition, "[0-9]"))){
    input_files[["annotation"]][(str_starts(Condition, "[0-9]")), Condition := str_c("X.", Condition)]
  }
  
  # Check if "Gene Names" and "Protein names" columnare present. Otherwise add them
  # TODO: check if true for all cases
  if(!("Gene names" %in% colnames(input_files[["PEP"]]))){
    input_files[["PEP"]][, `Gene names` := tstrsplit(tstrsplit(`Leading razor protein`, "\\|", keep = 3)[[1]], "_", keep = 1)[[1]]]
    input_files[["PEP"]][, `Leading razor protein` := tstrsplit(`Leading razor protein`, "\\|", keep = 2)[[1]]]
  }
  if(!("Protein names" %in% colnames(input_files[["PEP"]]))){
    data("anno_uniprot", envir = environment())
    anno_uniprot[, `Gene Names`:=NULL]
    if(length(intersect(anno_uniprot$`Leading razor protein`, input_files[["PEP"]]$`Leading razor protein`)) == 0){
      input_files[["PEP"]][, `Protein names` := `Leading razor protein`]
    }else{
      input_files[["PEP"]] <- merge.data.table(input_files[["PEP"]], unique(anno_uniprot, by="Leading razor protein", fromLast=FALSE), 
                                               by="Leading razor protein", 
                                               all.x = TRUE)
    }
  }
  if(nrow(input_files[["PEP"]][is.na(`Protein names`)])>0){
    data("anno_uniprot", envir = environment())
    anno_uniprot[, `Gene Names`:=NULL]
    unique_anno_uniprot <- unique(anno_uniprot, by="Leading razor protein", fromLast=FALSE)
    setnames(unique_anno_uniprot, old=c("Protein names"), new = c("desc"))
    input_files[["PEP"]][, `Protein names` := as.character(`Protein names`)]
    input_files[["PEP"]][is.na(`Protein names`), `Protein names` := unique_anno_uniprot[.SD, on=.(`Leading razor protein`), x.desc]]
  }
  
  # Manage Peptide file
  message("Cleaning data...")
  initial_peptide = nrow(input_files[["PEP"]])
  message(paste0("\tRaw number of peptide: ",initial_peptide))
  
  to_remove <- nrow(input_files[["PEP"]][grepl("Keratin|keratin", `Protein names`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!grepl("Keratin|keratin", `Protein names`)]
  message(paste0("\tKeratin peptide removed: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][grepl("CON_",`Leading proteins`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!grepl("CON_",`Leading proteins`)]
  message(paste0("\tCONTAMINANT peptide removed: ",to_remove))
  
  # to_remove <- nrow(input_files[["PEP"]][is.na(`Protein names`)])
  # input_files[["PEP"]] <- input_files[["PEP"]][!is.na(`Protein names`)]
  # message(paste0("\tPeptide removed due to missing Protein Names: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][is.na(`Gene names`) | stri_isempty(`Gene names`) | `Gene names` == ""])
  input_files[["PEP"]] <- input_files[["PEP"]][!(is.na(`Gene names`) | stri_isempty(`Gene names`) | `Gene names` == "")]
  message(paste0("\tPeptide removed due to missing Gene Names: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][!(`Raw file` %in% input_files[["annotation"]]$Sample)])
  input_files[["PEP"]] <- input_files[["PEP"]][(`Raw file` %in% input_files[["annotation"]]$Sample)]
  input_files[["annotation"]] <- input_files[["annotation"]][Sample %in% input_files[["PEP"]]$`Raw file`]
  message(paste0("\tPeptide removed since sample not in sample annotation: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][(Type == "MSMS")])
  input_files[["PEP"]] <- input_files[["PEP"]][(Type != "MSMS")]
  message(paste0("\tPeptide removed since sample not in sample annotation: ",to_remove))
  
  # Filter for phosphorylated modifications
  if(keep_only_phosphomodification){
    to_remove <- nrow(input_files[["PEP"]][!(grepl("Phospho|phospho", Modifications))])
    input_files$PEP <- input_files$PEP[grepl("Phospho|phospho", Modifications)]
    message(paste0("\tPeptide removed due to missing phosphorilation: ",to_remove))
  }
  
  # Remove phosphorylation with low score
  to_remove <- nrow(input_files[["PEP"]][!(sapply(str_extract_all(`Phospho (STY) Probabilities`, "1|(0\\.\\d+)"), 
                                                  function(x) any(as.numeric(x) > (phospho_thr))))])
  input_files$PEP <- input_files$PEP[sapply(str_extract_all(`Phospho (STY) Probabilities`, "1|(0\\.\\d+)"), 
                                            function(x) any(as.numeric(x) > (phospho_thr)))]
  message(paste0("\tPeptide removed due to low phosphorilation level (thr > ",phospho_thr,"): ",to_remove))
  
  # Modify sequence to detect only the phosphorylated sites
  input_files$PEP[, Sequence := sapply(1:.N, function(x) {
    pattern <- str_extract_all(`Phospho (STY) Probabilities`[x], "1|(0\\.\\d+)")
    replacement <- as.integer(as.numeric(unlist(pattern)) > (phospho_thr))
    stri_replace_all_regex(`Phospho (STY) Probabilities`[x], unlist(pattern), as.character(replacement), vectorize_all = FALSE)
  })]
  
  # Keep only first gene name
  input_files[["PEP"]][str_detect(`Gene names`,";"), `Gene names` := tstrsplit(`Gene names`, ";", keep = 1)]
  input_files[["PEP"]][, Modifications:=gsub(" ", "_", Modifications)]
  #Made the matrix
  psm_sig_raw <- data.table("ID_peptide" = as.factor(paste(input_files[["PEP"]]$`Gene names`, input_files[["PEP"]]$Sequence, input_files[["PEP"]]$Modifications, sep="_")), 
                            "Sample" = as.factor(input_files[["PEP"]]$`Raw file`), 
                            "Intensity" = input_files[["PEP"]]$Intensity)
  suppressWarnings({
    psm_sig_raw <- psm_sig_raw[, .(Intensity = sum(Intensity)), by = c("ID_peptide", "Sample")]
    psm_sig_raw <- dcast(psm_sig_raw, formula = ID_peptide~Sample, value.var = "Intensity")
  })
  colnames(psm_sig_raw)[-1] <- input_files[["annotation"]][match(colnames(psm_sig_raw)[-1], input_files[["annotation"]]$Sample)]$Sample
  
  #Made the description of the mpeptides
  psm_peptide_table <- as.data.table(unique(input_files[["PEP"]][, c("Leading razor protein", "Protein names", "Gene names", "Sequence", "Modifications", "Phospho (STY) Probabilities")]))
  colnames(psm_peptide_table) <- c("Accession","Description","GeneName","Annotated_Sequence","Modifications","Phospho_%")
  psm_peptide_table[, ID_peptide := paste(GeneName, Annotated_Sequence, Modifications, sep="_")]
  
  # Extract annotation
  c_anno<-input_files[["annotation"]]
  colnames(c_anno)<-tolower(colnames(c_anno))
  
  if(!("color" %in% colnames(c_anno))){
    message("Color column not found! Setting default color")
    c_anno<-merge.data.table(c_anno,
                             data.table("color"=colour_vec[1:length(unique(c_anno$condition))],
                                        "condition"=unique(c_anno$condition)),
                             by = "condition")
  }
  colour_vec<-na.omit(c_anno$color)
  names(colour_vec)<-na.omit(c_anno$sample)
  
  # Add chunk peptide annotation
  psm_anno_raw <- data.table("Accession" = psm_peptide_table$Accession,
                             "ID_peptide"=psm_peptide_table$ID_peptide,
                             "symbol"=psm_peptide_table$GeneName,
                             "sequence"=psm_peptide_table$Annotated_Sequence,
                             "modifications"=psm_peptide_table$Modifications)
  
  psm_anno_raw$old_id<-psm_anno_raw$ID_peptide
  
  # create peptide names linked to symbols
  psm_anno_raw[, row := seq(1,nrow(psm_anno_raw))]
  psm_anno_raw <- psm_anno_raw %>% group_by(symbol) %>% mutate("min"=min(row),"card"=n()) %>% ungroup()
  psm_anno_raw$rank<- psm_anno_raw$row - psm_anno_raw$min +1
  psm_anno_raw$id<- paste(psm_anno_raw$symbol,psm_anno_raw$rank,psm_anno_raw$card,sep="_")
  psm_anno_raw <- as.data.table(psm_anno_raw)
  
  
  # Convert to data.table
  psm_sig_prot_raw <- as.data.table(psm_sig_raw)
  psm_sig_pet_raw <- as.data.table(psm_sig_raw)
  
  # Preprocess Protein intensities
  psm_sig_prot_raw[psm_sig_prot_raw == 0] <- NA  # Transform 0s into NAs
  sig_thr <- filt_absent_value  # NA threshold
  
  psm_long_dt <- melt(psm_sig_prot_raw, id.vars = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_prot_df <- psm_sig_prot_raw[ID_peptide %in% filter_ID_peptide]
  psm_anno_df <- as.data.table(psm_anno_raw)[ID_peptide %in% filter_ID_peptide]
  
  # Log2 transformation
  psm_log_prot_df <- copy(psm_sig_prot_df)
  psm_log_prot_df[, (setdiff(names(psm_log_prot_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_prot_df), "ID_peptide")]
  
  # Filter proteins with only 1 peptID_peptidee
  psm_anno_df[, N_peptide := .N, by = symbol]
  filter_df_single_pep <- psm_anno_df[N_peptide >= min_peptide_protein, ID_peptide]
  psm_log_prot_df <- psm_log_prot_df[ID_peptide %in% filter_df_single_pep]
  psm_anno_df <- psm_anno_df[ID_peptide %in% filter_df_single_pep]
  
  # Preprocess peptID_peptidee intensities
  psm_sig_pet_raw[psm_sig_pet_raw == 0] <- NA  # Transform 0s into NAs
  
  psm_long_dt <- melt(psm_sig_pet_raw, id = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_pet_df <- psm_sig_pet_raw[ID_peptide %in% filter_ID_peptide]
  psm_peptide_table <- as.data.table(psm_peptide_table)[ID_peptide %in% filter_ID_peptide]
  
  # Determine tryptic condition
  # TODO: check this triptic test
  peptides_df <- psm_peptide_table[, .(Accession, Annotated_Sequence)]
  peptides_df[, preAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,1], regex("\\[\\w+\\]")),
                                 2,
                                 str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                        n = 3)[,1], regex("\\[\\w+\\]")))-1)]
  peptides_df[, endAA := substr(Annotated_Sequence, nchar(Annotated_Sequence)-1, nchar(Annotated_Sequence)-1)]
  peptides_df[, postAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,3], regex("\\[\\w+\\]")),
                                  2,
                                  str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                         n = 3)[,3], regex("\\[\\w+\\]")))-1)]
  
  peptides_df[, fully_TRI := preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, NSEMI_TRI := preAA %in% c("K","R") & !endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, CSEMI_TRI := !preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, non_TRI := !fully_TRI & !NSEMI_TRI & !CSEMI_TRI]
  
  peptides_df[, tryptic_cond := fifelse(fully_TRI, "fully tryptic", fifelse(NSEMI_TRI, "N-semi tryptic", fifelse(CSEMI_TRI, "C-semi tryptic", "non tryptic")))]
  psm_peptide_table[, tryptic_cond := peptides_df$tryptic_cond]
  
  # Log2 transformation for peptID_peptidees
  psm_log_pet_df <- copy(psm_sig_pet_df)
  psm_log_pet_df[, (setdiff(names(psm_log_pet_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_pet_df), "ID_peptide")]
  
  psm_log_pet_df <- psm_log_pet_df[ID_peptide %in% filter_df_single_pep]
  psm_peptide_table <- psm_peptide_table[ID_peptide %in% filter_df_single_pep]
  
  message(paste0("N Proteins (after filter): ", uniqueN(psm_peptide_table$Accession)))
  message(paste0("N Peptides (after filter): ", uniqueN(psm_peptide_table$`Phospho_%`)))
  
  return(list("c_anno" = c_anno,
              "psm_anno_df" = psm_anno_df,
              "psm_log_prot_df" = psm_log_prot_df,
              "psm_peptide_table" = psm_peptide_table,
              "psm_log_pet_df" = psm_log_pet_df,
              "colour_vec" = colour_vec))
}

# Read Proteome Discoverer
read_phospho_PD_files <- function(anno_filename, pep_filename, prot_filename, psm_file_filename,
                                  keep_only_phosphomodification = T, phospho_thr = 0.75,
                                  condition_col="Condition", sample_col="Sample", 
                                  color_col="Color", batch_corr_exe = FALSE, batch_col="batch",
                                  filt_absent_value = 0, min_peptide_protein = 0){
  input_files <- list()
  
  message("Reading files...")
  # Read annotation file
  input_files[["annotation"]] <- tryCatch({
    as.data.table(read_xlsx(anno_filename))
  }, error=function(cond){
    stop(paste0("Missing file. The file \'ANNOTATION\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  
  # Read peptide file
  input_files[["PEP"]] <- tryCatch({
    as.data.table(read_xlsx(pep_filename))
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PEPTIDE\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  # Read protein file
  input_files[["PROT"]] <- tryCatch({
    as.data.table(read_xlsx(prot_filename))
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PROTEIN\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  # Read protein file
  input_files[["PSM"]] <- tryCatch({
    as.data.table(read_xlsx(psm_file_filename, guess_max = 10000))
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PSM\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  message("File read.")
  
  
  message("Starting preprocessing...")
  #Clean files and merge
  colToKeep<-intersect(colnames(input_files[["annotation"]]), c("File ID",condition_col, sample_col, color_col, batch_col))
  if(!(condition_col %in% colToKeep)){
    stop("\'Condition\' column  in \'ANNOTATION\' file.")
  }
  if(!(sample_col %in% colToKeep)){
    stop("\'Sample\' column  in \'ANNOTATION\' file.")
  }
  if(batch_corr_exe & !(batch_col %in% colToKeep)){ 
    stop(paste0("\'",batch_col,"\' column  in \'ANNOTATION\' file with batch correction activated"))
  }
  input_files[["annotation"]] <- input_files[["annotation"]][, ..colToKeep]
  setnames(input_files[["annotation"]], old = "File ID", new = "File_ID")
  setnames(input_files[["annotation"]], old = condition_col, new = "Condition")
  setnames(input_files[["annotation"]], old = sample_col, new = "Sample")
  
  # Cast to string if all condition are numeric
  input_files[["annotation"]][, Condition := as.character(Condition)]
  
  # Check if conditions start with number. It generate problem later in DEqMS
  if(any(str_starts(input_files[["annotation"]]$Condition, "[0-9]"))){
    input_files[["annotation"]][(str_starts(Condition, "[0-9]")), Condition := str_c("X.", Condition)]
  }
  
  if(!("Modifications" %in% names(input_files$PEP))){
    setnames(input_files$PEP, grep("Modification",names(input_files$PEP), value = T), "Modifications")
  }
  
  # Prepare PROT table
  message("Cleaning data...")
  initial_peptide = nrow(input_files[["PROT"]])
  message(paste0("\tRaw number of protein: ",initial_peptide))
  
  to_remove <- nrow(input_files[["PROT"]][grepl("Keratin|keratin", Description)])
  input_files[["PROT"]] <- input_files[["PROT"]][!grepl("Keratin|keratin", Description)]
  message(paste0("\tKeratin protein removed: ",to_remove))
  
  to_remove <- nrow(input_files[["PROT"]][grepl("CON_",Accession) | grepl("CON_",Description)])
  input_files[["PROT"]] <- input_files[["PROT"]][!(grepl("CON_",Accession) | grepl("CON_",Description))]
  message(paste0("\tCONTAMINANT protein removed: ",to_remove))
  
  input_files$PROT <- na.omit(input_files$PROT[, .(Accession, Description)])
  
  input_files$PROT[, GeneName := sub(".*GN=(.+) PE=.*", "\\1", Description)]
  
  to_remove <- nrow(input_files[["PROT"]][is.na(GeneName) | stri_isempty(GeneName) | GeneName == ""])
  input_files[["PROT"]] <- input_files[["PROT"]][!(is.na(GeneName) | stri_isempty(GeneName) | GeneName == "")]
  message(paste0("\tProtein removed due to missing Gene Names: ",to_remove))
  
  # Maintain only the first UNIPROT code
  input_files$PEP[, `Master Protein Accessions` := sapply(`Master Protein Accessions`, function(x) strsplit(x, ";")[[1]][1])]
  
  # Filter for phosphorylated modifications
  if(keep_only_phosphomodification){
    to_remove <- nrow(input_files[["PEP"]][!(grepl("Phospho|phospho", Modifications))])
    input_files$PEP <- input_files$PEP[grepl("Phospho|phospho", Modifications)]
    message(paste0("\tPeptide removed due to missing phosphorilation: ",to_remove))
  }
  
  input_files[["PEP"]][, Modifications:=gsub(" ", "_", Modifications)]
  
  # Verify and clean PEP table
  abundance_cols <- grep("Abundance", names(input_files$PEP), value = TRUE)
  input_files$PEP <- input_files$PEP[, c("Annotated Sequence", "Modifications", 
                                         "Master Protein Accessions", "Positions in Master Proteins", abundance_cols), with = FALSE]
  
  # Modify sequence to detect only active phosphosite
  pattern_sub <- lapply(regmatches(input_files$PEP$Modifications, gregexpr("\\(\\d+\\.\\d+|\\(\\d+", input_files$PEP$Modifications)), function(x) gsub("\\(", "", x))
  replacement <- lapply(pattern_sub, function(x) as.integer(as.numeric(x) > (phospho_thr*100)))
  
  input_files$PEP[, new_seq := unlist(lapply(1:length(input_files$PEP$`Modifications`), 
                                             function(x){
                                               tryCatch({
                                                 stri_replace_all_regex(input_files$PEP$`Modifications`[x], 
                                                                        pattern_sub[[x]],
                                                                        as.character(replacement[[x]]),
                                                                        vectorize_all = F)
                                               }, error=function(cond){input_files$PEP$`Modifications`[x]})
                                             }))]
  
  input_files$PEP[, ID_peptide := paste(`Master Protein Accessions`, `Annotated Sequence`, new_seq, sep = "_")]
  
  # Collapse PEP table
  PEP_collapse <- input_files$PEP[, lapply(.SD, sum, na.rm = TRUE), by = ID_peptide, .SDcols = patterns("Abundance")]
  PEP_collapse <- merge(input_files$PEP[, .(`Annotated Sequence`, Modifications, `Master Protein Accessions`, `Positions in Master Proteins`, ID_peptide, new_seq)], PEP_collapse, by = "ID_peptide", all.y = TRUE)
  PEP_collapse <- PEP_collapse[!duplicated(ID_peptide)]
  
  # Merge by Uniprot
  uniprot_to_take <- merge(input_files$PROT, PEP_collapse, by.x = "Accession", by.y = "Master Protein Accessions")
  input_files$PD_PEP_matrix <- merge(input_files$PROT[Accession %in% uniprot_to_take$Accession], PEP_collapse, by.x = "Accession", by.y = "Master Protein Accessions")
  
  # Rename and clean column names
  input_files$PD_PEP_matrix[, GeneName := sub("-.*", "", GeneName)]
  colnames(input_files$PD_PEP_matrix) <- gsub("\\.", "_", colnames(input_files$PD_PEP_matrix))
  
  # Filter sample in file annotation
  colToKeep <- c(c("Accession","Description","GeneName","Annotated Sequence","Modifications",
                   "Positions in Master Proteins","ID_peptide","new_seq"),
                 lapply(input_files$annotation$File_ID, function(x){
                   colnames(input_files$PD_PEP_matrix)[grepl(str_c(x,":"), colnames(input_files$PD_PEP_matrix))]
                 }) %>% unlist())
  input_files$PD_PEP_matrix<-input_files$PD_PEP_matrix[, ..colToKeep]
  
  colnames(input_files$PD_PEP_matrix)<-c(c("Accession","Description","GeneName","Annotated Sequence","Modifications","Position in Master Proteins","ID_peptide","new_seq"),
                                         unlist(lapply(colnames(input_files$PD_PEP_matrix), function(x){
                                           lapply(input_files$annotation$File_ID, function(y){if(str_detect(x,paste0(y,":"))){
                                             input_files$annotation$Sample[which(input_files$annotation$File_ID %in% y)]
                                           }})
                                         })))
  
  colnames(input_files$PD_PEP_matrix) <- gsub(" ", "_", colnames(input_files$PD_PEP_matrix))
  
  #Made the description of the mpeptides
  psm_peptide_table <- as.data.table(unique(input_files[["PD_PEP_matrix"]][, c("ID_peptide","Accession","Description","GeneName","Annotated_Sequence","Modifications","Position_in_Master_Proteins")]))
  
  # Read PSM and filter Phospho Percentage
  if("ptmRS: Best Site Probabilities" %in% names(input_files[["PSM"]])){
    input_files[["PSM"]] <- input_files[["PSM"]][
      !(`Annotated Sequence` %in% input_files[["PD_PEP_matrix"]]$Annotated_Sequence) &
        grepl("Phospho|phospho", input_files[["PSM"]]$`ptmRS: Best Site Probabilities`) &
        !is.na(input_files[["PSM"]]$`Precursor Abundance`),]
    
    input_files[["PEP"]][, Modifications:=gsub(" ", "_", Modifications)]
    
    input_files[["PSM"]] <- input_files[["PSM"]][unlist(lapply(str_extract_all(input_files[["PSM"]]$`ptmRS: Best Site Probabilities`,
                                                                               "\\: \\d+\\.\\d+\\b|\\: \\d+\\b"),
                                                               function(x){any(as.double(unlist(str_remove_all(unlist(x), "\\: "))) > phospho_thr*100)})),]
    
    # Keep only the first Uniprot code
    input_files$PSM[, `Master Protein Accessions` := sapply(`Master Protein Accessions`, function(x) strsplit(x, ";")[[1]][1])]
    
    # Modify sequence to detect phosphorylated sites
    pattern_sub <- lapply(regmatches(input_files$PSM$`ptmRS: Best Site Probabilities`, gregexpr("\\: \\d+\\.\\d+\\b|\\: \\d+\\b", 
                                                                                                input_files$PSM$`ptmRS: Best Site Probabilities`)), 
                          function(x) gsub("\\: ", "", x))
    replacement <- lapply(pattern_sub, function(x) as.integer(as.numeric(x) > phospho_thr*100))
    
    input_files$PSM[, ID_peptide := unlist(lapply(1:length(input_files$PSM$`ptmRS: Best Site Probabilities`),
                                                  function(x){
                                                    stri_replace_all_regex(input_files$PSM$`ptmRS: Best Site Probabilities`[x],
                                                                           pattern_sub[[x]],
                                                                           as.character(replacement[[x]]),
                                                                           vectorize_all = F)}))]
  } else if("PhosphoRS: Best Site Probabilities" %in% names(input_files[["PSM"]])){
    input_files[["PSM"]] <- input_files[["PSM"]][
      !(`Annotated Sequence` %in% input_files[["PD_PEP_matrix"]]$Annotated_Sequence) &
        grepl("Phospho|phospho", input_files[["PSM"]]$`PhosphoRS: Best Site Probabilities`) &
        !is.na(input_files[["PSM"]]$`Precursor Abundance`),]
    
    input_files[["PEP"]][, Modifications:=gsub(" ", "_", Modifications)]
    
    input_files[["PSM"]] <- input_files[["PSM"]][unlist(lapply(str_extract_all(input_files[["PSM"]]$`PhosphoRS: Best Site Probabilities`,
                                                                               "\\: \\d+\\.\\d+\\b|\\: \\d+\\b"),
                                                               function(x){any(as.double(unlist(str_remove_all(unlist(x), "\\: "))) > phospho_thr*100)})),]
    
    # Keep only the first Uniprot code
    input_files$PSM[, `Master Protein Accessions` := sapply(`Master Protein Accessions`, function(x) strsplit(x, ";")[[1]][1])]
    
    # Modify sequence to detect phosphorylated sites
    pattern_sub <- lapply(regmatches(input_files$PSM$`PhosphoRS: Best Site Probabilities`, gregexpr("\\: \\d+\\.\\d+\\b|\\: \\d+\\b", 
                                                                                                input_files$PSM$`PhosphoRS: Best Site Probabilities`)), 
                          function(x) gsub("\\: ", "", x))
    replacement <- lapply(pattern_sub, function(x) as.integer(as.numeric(x) > phospho_thr*100))
    
    input_files$PSM[, ID_peptide := unlist(lapply(1:length(input_files$PSM$`PhosphoRS: Best Site Probabilities`),
                                                  function(x){
                                                    stri_replace_all_regex(input_files$PSM$`PhosphoRS: Best Site Probabilities`[x],
                                                                           pattern_sub[[x]],
                                                                           as.character(replacement[[x]]),
                                                                           vectorize_all = F)}))]
  }
  
  colnames(input_files$PSM) <- gsub(" ", "_", colnames(input_files$PSM))
  input_files$PSM[, Annotated_Sequence := toupper(Annotated_Sequence)]
  
  # Resolve ambiguous phosphosites
  id_pep_doubt <- psm_peptide_table[grepl("\\/", Modifications), .(Accession, Description, GeneName, Annotated_Sequence, Modifications, Position_in_Master_Proteins, ID_peptide)]
  
  pattern <- input_files$PSM[Annotated_Sequence %in% id_pep_doubt$Annotated_Sequence & Master_Protein_Accessions %in% id_pep_doubt$Accession, .(Annotated_Sequence, Master_Protein_Accessions, ID_peptide)]
  pattern[, modif_split := unlist(lapply(stri_split_regex(pattern$ID_peptide, ";"), 
                                         function(x){stri_c(unlist(lapply(stri_extract_all_regex(str_remove_all(x, "\\ "), "^.|.$"), 
                                                                          function(y){y[y[2] != "0"]})), collapse = "")}))]
  
  pattern2 <- pattern[id_pep_doubt, on = c("Annotated_Sequence"="Annotated_Sequence", "Master_Protein_Accessions"="Accession"), mult = "first"]
  pattern2[, phospo_PEP := unlist(lapply(stri_extract_all_regex(i.ID_peptide, "\\w\\d+\\(1|\\w\\/"), 
                                         function(x){stri_c(unlist(stri_extract_all_regex(x, "^.|.$")), collapse = "")}))]
  
  # Adjust doubtful phosphosites
  for (ids in unique(pattern2$i.ID_peptide)) {
    test <- pattern2[i.ID_peptide == ids & modif_split != phospo_PEP]
    if (nrow(test) > 0) {
      test$doubtSite<-lapply(1:nrow(test), function(x){str_remove(test$modif_split[x], test$phospo_PEP[x])})
      table_site <- table(unlist(test$doubtSite))
      site_to_take <- stri_extract_all_regex(names(table_site[which(table_site == max(table_site))])[1], "^.")[[1]]
      input_files$PD_PEP_matrix$Modifications[which(input_files$PD_PEP_matrix$ID_peptide == test$i.ID_peptide[1])]<- stri_replace_all_regex(test$Modifications[1], "\\w\\/\\w\\/\\w|\\w\\/\\w",paste0(site_to_take,"(100)"))
      input_files$PD_PEP_matrix$new_seq[which(input_files$PD_PEP_matrix$ID_peptide == test$i.ID_peptide[1])]<- stri_replace_all_regex(test$Modifications[1], "\\w\\/\\w\\/\\w|\\w\\/\\w",paste0(site_to_take,"(100)"))
      input_files$PD_PEP_matrix$ID_peptide[which(input_files$PD_PEP_matrix$ID_peptide == test$i.ID_peptide[1])]<- stri_replace_all_regex(test$i.ID_peptide[1], "\\w\\/\\w\\/\\w|\\w\\/\\w",paste0(site_to_take,"(100)"))
    }
  }
  
  
  #Maintain only the new complete matrix of abundance
  input_files$PD_PEP_matrix[, ID_peptide := as.factor(paste(input_files[["PD_PEP_matrix"]]$GeneName, input_files[["PD_PEP_matrix"]]$Annotated_Sequence, input_files[["PD_PEP_matrix"]]$new_seq, sep="_"))]
  colToKeep <- c("ID_peptide",input_files$annotation$Sample)
  psm_sig_raw<-input_files$PD_PEP_matrix[, ..colToKeep]
  
  #Made the description of the mpeptides
  psm_peptide_table <- as.data.table(unique(input_files[["PD_PEP_matrix"]][, c("ID_peptide","Accession","Description","GeneName","Annotated_Sequence","Modifications","Position_in_Master_Proteins")]))
  
  # Extract annotation
  c_anno<-input_files[["annotation"]]
  colnames(c_anno)<-tolower(colnames(c_anno))
  
  if(!("color" %in% colnames(c_anno))){
    message("Color column not found! Setting default color")
    c_anno<-merge.data.table(c_anno,
                             data.table("color"=colour_vec[1:length(unique(c_anno$condition))],
                                        "condition"=unique(c_anno$condition)),
                             by = "condition")
  }
  colour_vec<-na.omit(c_anno$color)
  names(colour_vec)<-na.omit(c_anno$sample)
  
  # Add chunk peptide annotation
  psm_anno_raw <- data.table("Accession" = psm_peptide_table$Accession,
                             "ID_peptide"=psm_peptide_table$ID_peptide,
                             "symbol"=psm_peptide_table$GeneName,
                             "sequence"=psm_peptide_table$Annotated_Sequence,
                             "modifications"=psm_peptide_table$Modifications)
  
  psm_anno_raw$old_id<-psm_anno_raw$ID_peptide
  
  # create peptide names linked to symbols
  psm_anno_raw[, row := seq(1,nrow(psm_anno_raw))]
  psm_anno_raw <- psm_anno_raw %>% group_by(symbol) %>% mutate("min"=min(row),"card"=n()) %>% ungroup()
  psm_anno_raw$rank<- psm_anno_raw$row - psm_anno_raw$min +1
  psm_anno_raw$id<- paste(psm_anno_raw$symbol,psm_anno_raw$rank,psm_anno_raw$card,sep="_")
  psm_anno_raw <- as.data.table(psm_anno_raw)
  
  
  # Convert to data.table
  psm_sig_prot_raw <- as.data.table(psm_sig_raw)
  psm_sig_pet_raw <- as.data.table(psm_sig_raw)
  
  # Preprocess Protein intensities
  psm_sig_prot_raw[psm_sig_prot_raw == 0] <- NA  # Transform 0s into NAs
  sig_thr <- filt_absent_value  # NA threshold
  
  psm_long_dt <- melt(psm_sig_prot_raw, id.vars = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_prot_df <- psm_sig_prot_raw[ID_peptide %in% filter_ID_peptide]
  psm_anno_df <- as.data.table(psm_anno_raw)[ID_peptide %in% filter_ID_peptide]
  
  # Log2 transformation
  psm_log_prot_df <- copy(psm_sig_prot_df)
  psm_log_prot_df[, (setdiff(names(psm_log_prot_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_prot_df), "ID_peptide")]
  
  # Filter proteins with only 1 peptID_peptidee
  psm_anno_df[, N_peptide := .N, by = symbol]
  filter_df_single_pep <- psm_anno_df[N_peptide >= min_peptide_protein, ID_peptide]
  psm_log_prot_df <- psm_log_prot_df[ID_peptide %in% filter_df_single_pep]
  psm_anno_df <- psm_anno_df[ID_peptide %in% filter_df_single_pep]
  
  # Preprocess peptID_peptidee intensities
  psm_sig_pet_raw[psm_sig_pet_raw == 0] <- NA  # Transform 0s into NAs
  
  psm_long_dt <- melt(psm_sig_pet_raw, id = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_long_dt[, is_na := is.na(counts)]
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is_na)), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, passes_c := min_c<=sig_thr]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(passes_c)), by = ID_peptide]
  
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_pet_df <- psm_sig_pet_raw[ID_peptide %in% filter_ID_peptide]
  psm_peptide_table <- as.data.table(psm_peptide_table)[ID_peptide %in% filter_ID_peptide]
  
  # Determine tryptic condition
  peptides_df <- psm_peptide_table[, .(Accession, Annotated_Sequence)]
  peptides_df[, preAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,1], regex("\\[\\w+\\]")),
                                 2,
                                 str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                        n = 3)[,1], regex("\\[\\w+\\]")))-1)]
  peptides_df[, endAA := substr(Annotated_Sequence, nchar(Annotated_Sequence)-1, nchar(Annotated_Sequence)-1)]
  peptides_df[, postAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,3], regex("\\[\\w+\\]")),
                                  2,
                                  str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                         n = 3)[,3], regex("\\[\\w+\\]")))-1)]
  
  peptides_df[, fully_TRI := preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, NSEMI_TRI := preAA %in% c("K","R") & !endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, CSEMI_TRI := !preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, non_TRI := !fully_TRI & !NSEMI_TRI & !CSEMI_TRI]
  
  peptides_df[, tryptic_cond := fifelse(fully_TRI, "fully tryptic", fifelse(NSEMI_TRI, "N-semi tryptic", fifelse(CSEMI_TRI, "C-semi tryptic", "non tryptic")))]
  psm_peptide_table[, tryptic_cond := peptides_df$tryptic_cond]
  
  # Log2 transformation for peptID_peptidees
  psm_log_pet_df <- copy(psm_sig_pet_df)
  psm_log_pet_df[, (setdiff(names(psm_log_pet_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_pet_df), "ID_peptide")]
  
  
  psm_log_pet_df <- psm_log_pet_df[ID_peptide %in% filter_df_single_pep]
  psm_peptide_table <- psm_peptide_table[ID_peptide %in% filter_df_single_pep]
  
  message(paste0("N Proteins (after filter): ", uniqueN(psm_peptide_table$Accession)))
  message(paste0("N Peptides (after filter): ", uniqueN(psm_peptide_table$ID_peptide)))
  
  return(list("c_anno" = c_anno,
              "psm_anno_df" = psm_anno_df,
              "psm_log_prot_df" = psm_log_prot_df,
              "psm_peptide_table" = psm_peptide_table,
              "psm_log_pet_df" = psm_log_pet_df,
              "colour_vec" = colour_vec))
}


#' Read Phospho-Proteome Data
#'
#' This function reads both proteome and phosphoproteome data from files based on the specified software (Proteome Discoverer or MaxQuant). 
#' It supports reading annotation, peptide, protein group, and PSM data for both proteome and phosphoproteome, with options for batch correction, 
#' filtering absent values, and keeping only phospho-modified peptides above a specified threshold.
#'
#' @param software Character; the proteomics software used, either "PD" for Proteome Discoverer or "MQ" for MaxQuant.
#' @param folder_proteome Character; the folder containing proteome data files.
#' @param folder_phospho Character; the folder containing phosphoproteome data files.
#' @param keep_only_phosphomodification Logical; whether to keep only phospho-modified peptides. Default is TRUE.
#' @param phospho_thr Numeric; the threshold for phospho-modified peptides. Default is 0.75.
#' @param annotation_proteome_filename Character; the name or pattern of the proteome annotation file. Default is "annotation".
#' @param peptide_proteome_filename Character; the name or pattern of the proteome peptide file. Default is "pep".
#' @param proteinGroup_proteome_filename Character; the name or pattern of the proteome protein group file. Default is "prot".
#' @param condition_proteome_col Character; the column name representing condition information in the proteome annotation file. Default is "Condition".
#' @param sample_proteome_col Character; the column name representing sample information in the proteome annotation file. Default is "Sample".
#' @param color_proteome_col Character; the column name representing color information for the plot in the proteome annotation file. Default is "Color".
#' @param annotation_phospho_filename Character; the name or pattern of the phosphoproteome annotation file. Default is "annotation".
#' @param peptide_phospho_filename Character; the name or pattern of the phosphoproteome peptide file. Default is "pep".
#' @param proteinGroup_phospho_filename Character; the name or pattern of the phosphoproteome protein group file. Default is "prot".
#' @param psm_phospho_filename Character; the name or pattern of the phosphoproteome PSM file. Default is "psm".
#' @param condition_phospho_col Character; the column name representing condition information in the phosphoproteome annotation file. Default is "Condition".
#' @param sample_phospho_col Character; the column name representing sample information in the phosphoproteome annotation file. Default is "Sample".
#' @param color_phospho_col Character; the column name representing color information for the plot in the phosphoproteome annotation file. Default is "Color".
#' @param batch_corr_exe Logical; whether batch correction should be applied. Default is FALSE.
#' @param batch_col Character; the column name representing batch information in the annotation file. Default is "batch".
#' @param filt_absent_value Numeric; the value used to filter out absent data. Default is 0.
#' @param min_peptide_protein Numeric; the value used to filter out protein by N peptide. Default is 0.
#'
#' @return A list containing both proteome and phosphoproteome data, including annotation and peptide data, 
#'         with batch correction applied if specified and only phospho-modified peptides kept if requested.
#'
#' @examples
#' \dontrun{
#' phospho_proteome_data <- read_phospho_proteome_proteomics(software = "MQ", 
#'                                                         folder_proteome = "proteome_folder", 
#'                                                         folder_phospho = "phospho_folder")
#' }
#'
#' @importFrom dplyr ungroup mutate filter group_by n
#' @import data.table
#' @import readxl
#' @export
read_phospho_proteome_proteomics <- function(software, 
                                             folder_proteome, folder_phospho, 
                                             keep_only_phosphomodification = T, phospho_thr = 0.75, 
                                             annotation_proteome_filename = "annotation",peptide_proteome_filename = "pep", 
                                             proteinGroup_proteome_filename = "prot",
                                             condition_proteome_col="Condition", sample_proteome_col="Sample", color_proteome_col="Color", 
                                             annotation_phospho_filename = "annotation",peptide_phospho_filename = "pep", 
                                             proteinGroup_phospho_filename = "prot", psm_phospho_filename = "psm",
                                             condition_phospho_col="Condition", sample_phospho_col="Sample", color_phospho_col="Color", 
                                             batch_corr_exe = FALSE, batch_col="batch",
                                             filt_absent_value = 0, min_peptide_protein = 0){
  
  if(!(software %in% c("PD","MQ"))){
    stop("Valid software is required. Write PD or MQ.",
         "\tPD: Proteome Discoverer",
         "\tMQ: MaxQuant")
  }
  
  # TODO: controllo che ci siano i file
  # TODO: batch correction verificare input
  
  # Verify file proteome
  anno_proteome_filename = list.files(folder_proteome, pattern = annotation_proteome_filename, full.names = T, ignore.case = T)
  if(is.null(anno_proteome_filename) | length(anno_proteome_filename) > 1){
    stop("Missing file proteome annotation or wrong annotation filename parameter or multiple files with same pattern")
  }
  pep_proteome_filename = list.files(folder_proteome, pattern = peptide_proteome_filename, full.names = T, ignore.case = T)
  if(is.null(pep_proteome_filename) | length(pep_proteome_filename) > 1){
    stop("Missing file proteome peptide or wrong proteome peptide filename parameter or multiple files with same pattern")
  }
  if(software == "PD"){
    prot_proteome_filename = list.files(folder_proteome, pattern = proteinGroup_proteome_filename, full.names = T, ignore.case = T)
    if(is.null(prot_proteome_filename) | length(prot_proteome_filename) > 1){
      stop("Missing file proteome proteinGroup or wrong proteome proteinGroup filename parameter or multiple files with same pattern")
    }
  } else{
    prot_proteome_filename = NULL
  }
  
  # Verify file phospho
  anno_phospho_filename = list.files(folder_phospho, pattern = annotation_phospho_filename, full.names = T, ignore.case = T)
  if(is.null(anno_phospho_filename) | length(anno_phospho_filename) > 1){
    stop("Missing file phospho annotation or wrong annotation filename parameter or multiple files with same pattern")
  }
  pep_phospho_filename = list.files(folder_phospho, pattern = peptide_phospho_filename, full.names = T, ignore.case = T)
  if(is.null(pep_phospho_filename) | length(pep_phospho_filename) > 1){
    stop("Missing file phospho peptide or wrong phospho peptide filename parameter or multiple files with same pattern")
  }
  if(software == "PD"){
    prot_phospho_filename = list.files(folder_phospho, pattern = proteinGroup_phospho_filename, full.names = T, ignore.case = T)
    if(is.null(prot_phospho_filename) | length(prot_phospho_filename) > 1){
      stop("Missing file phospho proteinGroup or wrong phospho proteinGroup filename parameter or multiple files with same pattern")
    }
  } else{
    prot_phospho_filename = NULL
  }
  if(software == "PD"){
    psm_phospho_file_filename = list.files(folder_phospho, pattern = psm_phospho_filename, full.names = T, ignore.case = T)
    if(is.null(psm_phospho_file_filename) | length(psm_phospho_file_filename) > 1){
      stop("Missing file PSM or wrong PSM filename parameter or multiple files with same pattern")
    }
  } else{
    psm_phospho_file_filename = NULL
  }
  
  # Read proteome data
  proteome_data = NULL
  message("Reading Proteome data...")
  if(software == "PD"){
    proteome_data = read_PD_files(anno_proteome_filename, pep_proteome_filename, prot_proteome_filename, 
                                  condition_col=condition_proteome_col, 
                                  sample_col=sample_proteome_col, 
                                  color_col=color_proteome_col,
                                  batch_corr_exe= batch_corr_exe, 
                                  filt_absent_value = filt_absent_value, 
                                  batch_col = batch_col,
                                  min_peptide_protein = min_peptide_protein)
  } else if(software == "MQ"){
    proteome_data = read_MQ_files(anno_proteome_filename, pep_proteome_filename, 
                                  condition_col=condition_proteome_col, 
                                  sample_col=sample_proteome_col, 
                                  color_col=color_proteome_col,
                                  batch_corr_exe = batch_corr_exe, 
                                  filt_absent_value = filt_absent_value, 
                                  batch_col = batch_col,
                                  min_peptide_protein = min_peptide_protein)
  }
  message("Read Proteome data!")
  
  # Read phospho data
  phospho_data = NULL
  message("Reading Phosphoproteomic data...")
  if(software == "PD"){
    phospho_data = read_phospho_PD_files(anno_phospho_filename, pep_phospho_filename, prot_phospho_filename, psm_phospho_file_filename, 
                                         keep_only_phosphomodification = keep_only_phosphomodification,
                                         condition_col=condition_phospho_col, 
                                         sample_col=sample_phospho_col, 
                                         color_col=color_phospho_col,
                                         batch_corr_exe= batch_corr_exe, 
                                         filt_absent_value = filt_absent_value, 
                                         phospho_thr = phospho_thr, 
                                         batch_col = batch_col,
                                         min_peptide_protein = min_peptide_protein)
  } else if(software == "MQ"){
    phospho_data = read_phospho_MQ_files(anno_phospho_filename, pep_phospho_filename, 
                                         keep_only_phosphomodification = keep_only_phosphomodification,
                                         condition_col=condition_phospho_col, 
                                         sample_col=sample_phospho_col, 
                                         color_col=color_phospho_col,
                                         batch_corr_exe = batch_corr_exe, 
                                         filt_absent_value = filt_absent_value,
                                         phospho_thr = phospho_thr,  
                                         batch_col = batch_col,
                                         min_peptide_protein = min_peptide_protein)
  }
  message("Reading Phosphoproteomic data!")
  
  phospho_proteome_data <- list("c_anno_proteome" = proteome_data$c_anno,
                                "c_anno_phospho" = phospho_data$c_anno,
                                "psm_anno_df" = proteome_data$psm_anno_df,
                                "psm_log_prot_df" = proteome_data$psm_log_prot_df,
                                "psm_peptide_table" = phospho_data$psm_peptide_table,
                                "psm_log_pet_df" = phospho_data$psm_log_pet_df,
                                "colour_vec_proteome" = proteome_data$colour_vec,
                                "colour_vec_phospho" = phospho_data$colour_vec)
  return(phospho_proteome_data)
}


#' Read Spatial Proteomics Data
#'
#' This function reads Spatial proteomics data from files based on the specified software (Proteome Discoverer or MaxQuant). 
#' It supports reading annotation, peptide, and protein group data, with options for batch correction and filtering absent values.
#'
#' @param software Character; the proteomics software used, either "PD" for Proteome Discoverer or "MQ" for MaxQuant.
#' @param folder Character; the folder containing the data files.
#' @param use_proteinGroups_MQ Logical; Only for MaxQuant. If FALSE it expect the evidence.txt file, if TRUE require peptide.txt and proteinGroups.txt files. Default is FALSE-
#' @param peptide_filename Character; the name or pattern of the peptide file. Default is "pep".
#' @param annotation_filename Character; the name or pattern of the annotation file. Default is "annotation".
#' @param sample_col Character; the column name representing sample information in the annotation file. Default is "Sample".
#' @param color_col Character; the column name representing color information for the plot in the annotation file. Default is "Color".
#' @param row_plate_col Character; the column name representing y coordinate of the sample in the annotation file. Default is "y".
#' @param column_plate_col Character; the column name representing y coordinate of the sample in the annotation file. Default is "x".
#' @param proteinGroup_filename Character; the name or pattern of the protein group file. Default is "prot".
#' @param filt_absent_value Numeric; the value used to filter out absent data. Default is 0.
#' @param min_peptide_protein Numeric; the value used to filter out protein by N peptide. Default is 0.
#'
#' @return A list containing Spatial proteomics data, including annotation and peptide data, 
#'         with batch correction applied if specified.
#'
#' @examples
#' \dontrun{
#' proteome_data <- read_spatial_proteomics(software = "MQ", folder = "data_folder")
#' }
#'
#' @importFrom dplyr ungroup mutate filter group_by n
#' @import data.table
#' @import stringi
#' @import stringr
#' @export
read_spatial_proteomics <- function(software, folder, use_proteinGroups_MQ = TRUE, 
                                    annotation_filename = "annotation", peptide_filename = "pep", proteinGroup_filename = "prot",
                                    sample_col="Sample", color_col="Color", row_plate_col = "y", column_plate_col = "x",
                                    filt_absent_value = 0,min_peptide_protein = 0){
  if(!(software %in% c("MQ"))){
    stop("Valid software is required. Write MQ.",
         "\tMQ: MaxQuant")
  }
  
  anno_filename = list.files(folder, pattern = annotation_filename, full.names = T, ignore.case = T)
  if(length(anno_filename) > 1){
    stop("Multiple annotation files with same pattern")
  } 
  pep_filename = list.files(folder, pattern = peptide_filename, full.names = T, ignore.case = T)
  if(is.null(pep_filename) | length(pep_filename) > 1){
    stop("Missing file peptide or wrong peptide filename parameter or multiple files with same pattern")
  }
  if(software == "PD"){
    prot_filename = list.files(folder, pattern = proteinGroup_filename, full.names = T, ignore.case = T)
    if(is.null(prot_filename) | length(prot_filename) > 1){
      stop("Missing file proteinGroup or wrong proteinGroup filename parameter or multiple files with same pattern")
    }
  } else if(use_proteinGroups_MQ){
    prot_filename = list.files(folder, pattern = proteinGroup_filename, full.names = T, ignore.case = T)
    if(is.null(prot_filename) | length(prot_filename) > 1){
      stop("Missing file proteinGroup or wrong proteinGroup filename parameter or multiple files with same pattern")
    }
  }else{
    prot_filename = NULL
  }
  
  proteome_data = NULL
  if(software == "MQ" & !use_proteinGroups_MQ){
    message("Reading MaxQuant files (evidence.txt)...")
    proteome_data = read_spatial_MQ_files(anno_filename, pep_filename,
                                          filt_absent_value = filt_absent_value,
                                          min_peptide_protein = min_peptide_protein)
  } else if(software == "MQ" & use_proteinGroups_MQ){
    message("Reading MaxQuant files (peptides.txt & proteinGroups.txt)...")
    proteome_data = read_spatial_MQ_prot_peptide_files(anno_filename, pep_filename, prot_filename,
                                                       filt_absent_value = filt_absent_value, 
                                                       min_peptide_protein = min_peptide_protein,
                                                       sample_col = sample_col,
                                                       color_col = color_col,
                                                       row_plate_col = row_plate_col,
                                                       column_plate_col = column_plate_col)
  }
  
  return(proteome_data)
}


# Read MaxQuant files Spatial proteomics
read_spatial_MQ_prot_peptide_files <- function(anno_filename, pep_filename, prot_filename,
                                               sample_col, color_col, row_plate_col, column_plate_col, 
                                       filt_absent_value = 0, min_peptide_protein = 0){
  input_files <- list()
  
  message("Reading files...")
  
  if(!is.null(anno_filename)){
    input_files[["annotation"]] <- tryCatch({
      as.data.table(read_xlsx(anno_filename))
    }, error=function(cond){
      stop(paste0("Missing file. The file \'ANNOTATION\' is missing or not have the pattern in the filename or there are duplicates files."))
    })
  }
  
  # Read peptide file
  input_files[["PEP"]] <- tryCatch({
    fread(pep_filename)
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PEPTIDE\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  
  # Read peptide file
  input_files[["PROT"]] <- tryCatch({
    fread(prot_filename)
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PROTEIN_GROUPS\' is missing or not have the pattern in the filename or there are duplicates files."))
  })
  message("File read.")
  message("Starting preprocessing...")
  
  #Clean files and merge
  colToKeep<-intersect(colnames(input_files[["annotation"]]), c(sample_col, color_col, row_plate_col, column_plate_col))
  if(!(sample_col %in% colToKeep)){
    stop("\'Sample\' column  in \'ANNOTATION\' file.")
  }
  if(!(row_plate_col %in% colToKeep)){
    stop("\'Row plate position\' column  in \'ANNOTATION\' file.")
  }
  if(!(column_plate_col %in% colToKeep)){
    stop("\'Column plate position\' column  in \'ANNOTATION\' file.")
  }
  input_files[["annotation"]] <- input_files[["annotation"]][, ..colToKeep]
  setnames(input_files[["annotation"]], old = sample_col, new = "Sample")
  setnames(input_files[["annotation"]], old = row_plate_col, new = "y")
  setnames(input_files[["annotation"]], old = column_plate_col, new = "x")
  
  
  # Extract ID cell and annotation
  intensity_cols <- grep("LFQ intensity", names(input_files[["PROT"]]), value = T)
  samples <- stri_replace_all(intensity_cols, regex = "LFQ intensity ", replacement = "")
  message(paste0("Detected ",nrow(input_files[["annotation"]])," cells"))
  
  # Manage Peptide file
  message("Cleaning data...")
  initial_peptide <- nrow(input_files$PEP)
  initial_protein <- nrow(input_files$PROT)
  message(paste0("\tRaw number of protein: ",initial_protein))
  
  to_remove <- nrow(input_files[["PROT"]][grepl("Keratin|keratin", `Protein names`)])
  input_files[["PROT"]] <- input_files[["PROT"]][!grepl("Keratin|keratin", `Protein names`)]
  message(paste0("\tKeratin protein removed: ",to_remove))
  
  to_remove <- nrow(input_files[["PROT"]][grepl("CON_",`Majority protein IDs`)])
  input_files[["PROT"]] <- input_files[["PROT"]][!grepl("CON_",`Majority protein IDs`)]
  message(paste0("\tCONTAMINANT protein removed: ",to_remove))
  
  # to_remove <- nrow(input_files[["PROT"]][is.na(`Protein names`)])
  # input_files[["PROT"]] <- input_files[["PROT"]][!is.na(`Protein names`)]
  # message(paste0("\tProtein removed due to missing Protein Names: ",to_remove))
  
  to_remove <- nrow(input_files[["PROT"]][is.na(`Gene names`) | stri_isempty(`Gene names`) | `Gene names` == ""])
  input_files[["PROT"]] <- input_files[["PROT"]][!(is.na(`Gene names`) | stri_isempty(`Gene names`) | `Gene names` == "")]
  message(paste0("\tPeptide removed due to missing Gene Names: ",to_remove))
  
  # Keep only first gene name
  input_files[["PROT"]][str_detect(`Gene names`,";"), `Gene names` := tstrsplit(`Gene names`, ";", keep = 1)]
  input_files[["PROT"]][str_detect(`Majority protein IDs`,";"), `Majority protein IDs` := tstrsplit(`Majority protein IDs`, ";", keep = 1)]
  
  # Column Extraction (LFQ intensity Plate... & protein information) + merging in data
  intensity_cols <- grep("LFQ intensity", names(input_files[["PROT"]]), value = T)
  id_cols <- c("Majority protein IDs", "Protein names", "Gene names", "Peptides", "Unique peptides", "Peptide IDs",intensity_cols)
  input_files[["PROT"]] <- input_files[["PROT"]][, ..id_cols]
  setnames(input_files[["PROT"]], 
           old=c("Majority protein IDs", "Protein names", "Gene names", "Peptides", "Unique peptides","Peptide IDs"),
           new = c("Accession","Description","GeneName","N_peptides","N_unique_peptides","Peptide_IDs"))

  # Filter peptides
  message(paste0("\tRaw number of peptides: ",initial_peptide))
  pep_ids <- paste0(input_files[["PROT"]]$Peptide_IDs, collapse = ";")
  pep_ids_split <- unique(stri_split(pep_ids, regex = ";")[[1]])
  to_remove <- nrow(input_files[["PEP"]][!(id %in% pep_ids_split)])
  input_files[["PEP"]] <- input_files[["PEP"]][id %in% pep_ids_split]
  message(paste0("\tPeptide removed due to filtered proteins: ",to_remove))
  
  
  # Check if "Gene Names" and "Protein names" columnare present. Otherwise add them
  # TODO: check if true for all cases
  if(!("Gene names" %in% colnames(input_files[["PEP"]]))){
    input_files[["PEP"]][, `Gene names` := tstrsplit(tstrsplit(`Leading razor protein`, "\\|", keep = 3)[[1]], "_", keep = 1)[[1]]]
    input_files[["PEP"]][, `Leading razor protein` := tstrsplit(`Leading razor protein`, "\\|", keep = 2)[[1]]]
  }
  if(!("Protein names" %in% colnames(input_files[["PEP"]]))){
    data("anno_uniprot", envir = environment())
    anno_uniprot[, `Gene Names`:=NULL]
    if(length(intersect(anno_uniprot$`Leading razor protein`, input_files[["PEP"]]$`Leading razor protein`)) == 0){
      input_files[["PEP"]][, `Protein names` := `Leading razor protein`]
    }else{
      input_files[["PEP"]] <- merge.data.table(input_files[["PEP"]], unique(anno_uniprot, by="Leading razor protein", fromLast=FALSE), 
                                               by="Leading razor protein", 
                                               all.x = TRUE)
    }
  } 
  if(nrow(input_files[["PEP"]][is.na(`Protein names`)])>0){
    data("anno_uniprot", envir = environment())
    anno_uniprot[, `Gene Names`:=NULL]
    unique_anno_uniprot <- unique(anno_uniprot, by="Leading razor protein", fromLast=FALSE)
    setnames(unique_anno_uniprot, old=c("Protein names"), new = c("desc"))
    input_files[["PEP"]][, `Protein names` := as.character(`Protein names`)]
    input_files[["PEP"]][is.na(`Protein names`), `Protein names` := unique_anno_uniprot[.SD, on=.(`Leading razor protein`), x.desc]]
  }
  
  to_remove <- nrow(input_files[["PEP"]][grepl("Keratin|keratin", `Protein names`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!grepl("Keratin|keratin", `Protein names`)]
  message(paste0("\tKeratin peptide removed: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][grepl("CON_",`Leading razor protein`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!grepl("CON_",`Leading razor protein`)]
  message(paste0("\tCONTAMINANT peptide removed: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][is.na(`Gene names`) | stri_isempty(`Gene names`) | `Gene names` == ""])
  input_files[["PEP"]] <- input_files[["PEP"]][!(is.na(`Gene names`) | stri_isempty(`Gene names`) | `Gene names` == "")]
  message(paste0("\tPeptide removed due to missing Gene Names: ",to_remove))
  
  # Keep only first gene name
  input_files[["PEP"]][str_detect(`Gene names`,";"), `Gene names` := tstrsplit(`Gene names`, ";", keep = 1)]
  input_files[["PEP"]][, Modifications:=NA]
  
  # Column Extraction (LFQ intensity Plate... & protein information) + merging in data
  intensity_cols <- grep("Intensity ", names(input_files[["PEP"]]), value = T)
  id_cols <- c("Leading razor protein", "Protein names", "Gene names", "Sequence", "Modifications",intensity_cols)
  input_files[["PEP"]] <- input_files[["PEP"]][, ..id_cols]
  setnames(input_files[["PEP"]], 
           old=c("Leading razor protein", "Protein names", "Gene names", "Sequence"),
           new = c("Accession","Description","GeneName","Annotated_Sequence"))
  
  #Merge PROT PEP
  ids_prot <- unique(input_files$PROT$Accession)
  ids_pep <- unique(input_files$PEP$Accession)
  intersect_id <- intersect(ids_prot, ids_pep)
  
  to_remove <- nrow(input_files[["PROT"]][!(Accession %in% intersect_id)])
  input_files[["PROT"]] <- input_files[["PROT"]][Accession %in% intersect_id]

  to_remove <- nrow(input_files[["PEP"]][!(Accession %in% intersect_id)])
  input_files[["PEP"]] <- input_files[["PEP"]][Accession %in% intersect_id]

  # Collapse duplicated peptide
  suppressWarnings({
    input_files[["PEP_l"]] <- melt(input_files[["PEP"]], 
                                      id.vars = c("Accession","Description","GeneName","Annotated_Sequence", "Modifications"),
                                      measure.vars = grep("Intensity ", names(input_files$PEP), value = TRUE), 
                                      variable.name = "Cell_plate", 
                                      value.name = "Intensity")
  })
  input_files[["PEP_l"]][, `Cell_plate` := stri_replace(`Cell_plate`, regex = "Intensity ", replacement = "")]
  
  #Made the matrix
  psm_sig_raw <- data.table("ID_peptide" = as.factor(paste(input_files[["PEP_l"]]$GeneName, input_files[["PEP_l"]]$Annotated_Sequence, input_files[["PEP_l"]]$Modifications, sep="_")), 
                            "Sample" = as.factor(input_files[["PEP_l"]]$Cell_plate), 
                            "Intensity" = input_files[["PEP_l"]]$Intensity)
  suppressWarnings({
    psm_sig_raw <- psm_sig_raw[, .(Intensity = sum(Intensity)), by = c("ID_peptide", "Sample")]
    psm_sig_raw <- dcast(psm_sig_raw, formula = ID_peptide~Sample, value.var = "Intensity")
  })
  
  # if(any(duplicated(input_files[["PROT"]]$GeneName))){
  # Collapse duplicated peptide
  suppressWarnings({
    input_files[["PROT_l"]] <- melt(input_files[["PROT"]], 
                                    id.vars = c("Accession","Description","GeneName","N_peptides","N_unique_peptides","Peptide_IDs"),
                                    measure.vars = grep("LFQ intensity ", names(input_files$PROT), value = TRUE), 
                                    variable.name = "Cell_plate", 
                                    value.name = "Intensity")
  })
  input_files[["PROT_l"]][, `Cell_plate` := stri_replace(`Cell_plate`, regex = "LFQ intensity ", replacement = "")]
  
  #Made the matrix
  psm_sig_raw_prot <- data.table("ID_peptide" = as.factor(input_files[["PROT_l"]]$GeneName), 
                                 "Sample" = as.factor(input_files[["PROT_l"]]$Cell_plate), 
                                 "Intensity" = input_files[["PROT_l"]]$Intensity)
  suppressWarnings({
    psm_sig_raw_prot <- psm_sig_raw_prot[, .(Intensity = sum(Intensity)), by = c("ID_peptide", "Sample")]
    psm_sig_raw_prot <- dcast(psm_sig_raw_prot, formula = ID_peptide~Sample, value.var = "Intensity")
  })
  # }
  
  #Made the description of the mpeptides
  psm_peptide_table <- as.data.table(unique(input_files[["PEP"]][, c("Accession","Description","GeneName","Annotated_Sequence", "Modifications")]))
  protein_table <- as.data.table(unique(input_files[["PROT"]][, c("Accession","Description","GeneName","N_peptides","N_unique_peptides","Peptide_IDs")]))
  psm_peptide_table <- merge.data.table(psm_peptide_table, protein_table, by=c("Accession","Description","GeneName"), all.x = T)
  psm_peptide_table[, ID_peptide := paste(GeneName, Annotated_Sequence, Modifications, sep="_")]
  
  # Extract annotation
  c_anno<-input_files[["annotation"]]
  colnames(c_anno)<-tolower(colnames(c_anno))
  
  if(!("color" %in% colnames(c_anno))){
    message("Color column not found! Setting default color")
    c_anno<-merge.data.table(c_anno,
                             data.table("color"="#0078AEAA",
                                        "sample"=unique(c_anno$sample)),
                             by = "sample")
  }
  colour_vec<-na.omit(c_anno$color)
  names(colour_vec)<-na.omit(c_anno$sample)
  
  # Add chunk peptide annotation
  psm_anno_raw <- data.table("Accession" = psm_peptide_table$Accession,
                             "ID_peptide"=psm_peptide_table$ID_peptide,
                             "symbol"=psm_peptide_table$GeneName,
                             "sequence"=psm_peptide_table$Annotated_Sequence,
                             "modifications"=psm_peptide_table$Modifications)
  
  psm_anno_raw$old_id<-psm_anno_raw$ID_peptide
  
  # create peptide names linked to symbols
  psm_anno_raw[, row := seq(1,nrow(psm_anno_raw))]
  psm_anno_raw <- psm_anno_raw %>% group_by(symbol) %>% mutate("min"=min(row),"card"=n()) %>% ungroup()
  psm_anno_raw$rank<- psm_anno_raw$row - psm_anno_raw$min +1
  psm_anno_raw$id<- paste(psm_anno_raw$symbol,psm_anno_raw$rank,psm_anno_raw$card,sep="_")
  psm_anno_raw <- as.data.table(psm_anno_raw)
  
  # Convert to data.table
  psm_sig_prot_raw <- as.data.table(psm_sig_raw_prot)
  psm_sig_pet_raw <- as.data.table(psm_sig_raw)
  
  # Preprocess Protein intensities
  psm_sig_prot_raw[psm_sig_prot_raw == 0] <- NA  # Transform 0s into NAs

  psm_sig_prot_df <- psm_sig_prot_raw
  psm_anno_df <- as.data.table(psm_anno_raw)
  
  # Log2 transformation
  psm_log_prot_df <- copy(psm_sig_prot_df)
  psm_log_prot_df[, (setdiff(names(psm_log_prot_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_prot_df), "ID_peptide")]
  
  # Preprocess peptID_peptidee intensities
  psm_sig_pet_raw[psm_sig_pet_raw == 0] <- NA  # Transform 0s into NAs
  
  psm_sig_pet_df <- psm_sig_pet_raw
  psm_peptide_table <- as.data.table(psm_peptide_table)
  
  # Determine tryptic condition
  peptides_df <- psm_peptide_table[, .(Accession, Annotated_Sequence)]
  peptides_df[, preAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,1], regex("\\[\\w+\\]")),
                                 2,
                                 str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                        n = 3)[,1], regex("\\[\\w+\\]")))-1)]
  peptides_df[, endAA := substr(Annotated_Sequence, nchar(Annotated_Sequence)-1, nchar(Annotated_Sequence)-1)]
  peptides_df[, postAA := str_sub(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."), n = 3)[,3], regex("\\[\\w+\\]")),
                                  2,
                                  str_length(str_extract(str_split_fixed(Annotated_Sequence,regex("\\."),
                                                                         n = 3)[,3], regex("\\[\\w+\\]")))-1)]
  
  peptides_df[, fully_TRI := preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, NSEMI_TRI := preAA %in% c("K","R") & !endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, CSEMI_TRI := !preAA %in% c("K","R") & endAA %in% c("K","R") & (!postAA %in% "P" | is.na(postAA))]
  peptides_df[, non_TRI := !fully_TRI & !NSEMI_TRI & !CSEMI_TRI]
  
  peptides_df[, tryptic_cond := fifelse(fully_TRI, "fully tryptic", fifelse(NSEMI_TRI, "N-semi tryptic", fifelse(CSEMI_TRI, "C-semi tryptic", "non tryptic")))]
  psm_peptide_table[, tryptic_cond := peptides_df$tryptic_cond]
  
  # Log2 transformation for peptID_peptidees
  psm_log_pet_df <- copy(psm_sig_pet_df)
  psm_log_pet_df[, (setdiff(names(psm_log_pet_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_pet_df), "ID_peptide")]
  
  message(paste0("N Proteins (after filter): ", uniqueN(psm_log_prot_df$ID_peptide)))
  message(paste0("N Peptides (after filter): ", uniqueN(psm_peptide_table$ID_peptide)))
  
  return(list("c_anno" = c_anno,
              "psm_anno_df" = psm_anno_df,
              "psm_log_prot_df" = psm_log_prot_df,
              "psm_peptide_table" = psm_peptide_table,
              "psm_log_pet_df" = psm_log_pet_df,
              "colour_vec" = colour_vec))
}