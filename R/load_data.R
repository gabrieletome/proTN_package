#' From Proteome Discoverer or MaxQuant files to lists of data tables.
#'
#' This function reads proteomics data.
#'
#' @param software Must be \strong{PD} for Proteome Discoverer or \strong{MQ} for MaxQuant
#' @param folder Path with files
#' @param annotation_filename Must 
#' @param peptide_filename Must 
#' @param proteinGroup_filename Must
#' @param condition_col Must b
#' @param sample_col Must 
#' @param color_col Must
#' @param batch_corr_exe description
#' @param batch_col description
#' @param filt_absent_value description
#' @return A list of data tables.
#' @details \strong{ProTN}
#' @examples 
#' ## ## Example:
#' ## example
#' ## example2
#' @import data.table
#' @export
read_proteomics <- function(software, folder, annotation_filename = "annotation",
                            peptide_filename = "pep", proteinGroup_filename = "prot", 
                            condition_col="Condition", sample_col="Sample", 
                            color_col="Color", batch_corr_exe = FALSE, batch_col="batch",
                            filt_absent_value = 0){
  
  source("./R/functions.R")
  
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
  
  proteome_data = NULL
  if(software == "PD"){
    proteome_data = read_PD_files(anno_filename, pep_filename, prot_filename, 
                                  batch_corr_exe= batch_corr_exe, filt_absent_value = filt_absent_value)
  } else if(software == "MQ"){
    proteome_data = read_MQ_files(anno_filename, pep_filename, 
                                  batch_corr_exe = batch_corr_exe, filt_absent_value = filt_absent_value)
  }
  
  return(proteome_data)
}


# Read MaxQuant files
read_MQ_files <- function(anno_filename, pep_filename, 
                          condition_col="Condition", sample_col="Sample", 
                          color_col="Color", batch_corr_exe = FALSE, batch_col="batch",
                          filt_absent_value = 0){
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
    stop("\'Condition\' column missin in \'ANNOTATION\' file.")
  }
  if(!(sample_col %in% colToKeep)){
    stop("\'Sample\' column missin in \'ANNOTATION\' file.")
  }
  if(batch_corr_exe & !(batch_col %in% colToKeep)){ 
    stop(paste0("\'",batch_col,"\' column missin in \'ANNOTATION\' file with batch correction activated"))
  }
  input_files[["annotation"]] <- input_files[["annotation"]][, ..colToKeep]
  setnames(input_files[["annotation"]], old = condition_col, new = "Condition")
  setnames(input_files[["annotation"]], old = sample_col, new = "Sample")
  
  # Check if conditions start with number. It generate problem later in DEqMS
  if(any(str_starts(input_files[["annotation"]]$Condition, "[0-9]"))){
    input_files[["annotation"]][(str_starts(Condition, "[0-9]")), Condition := str_c("X.", Condition)]
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
  
  to_remove <- nrow(input_files[["PEP"]][is.na(`Protein names`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!is.na(`Protein names`)]
  message(paste0("\tPeptide removed due to missing Protein Names: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][is.na(`Gene names`) | stri_isempty(`Gene names`) | `Gene names` == ""])
  input_files[["PEP"]] <- input_files[["PEP"]][!(is.na(`Gene names`) | stri_isempty(`Gene names`) | `Gene names` == "")]
  message(paste0("\tPeptide removed due to missing Gene Names: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][!(`Raw file` %in% input_files[["annotation"]]$Sample)])
  input_files[["PEP"]] <- input_files[["PEP"]][(`Raw file` %in% input_files[["annotation"]]$Sample)]
  message(paste0("\tPeptide removed since sample not in sample annotation: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][(Type == "MSMS")])
  input_files[["PEP"]] <- input_files[["PEP"]][(Type != "MSMS")]
  message(paste0("\tPeptide removed since sample not in sample annotation: ",to_remove))
  
  # Keep only first gene name
  input_files[["PEP"]][str_detect(`Gene names`,";"), `Gene names` := tstrsplit(`Gene names`, ";", keep = 1)]
  input_files[["PEP"]][, Modifications:=gsub(" ", "_", Modifications)]
  #Made the matrix
  psm_sig_raw <- data.table("ID_peptide" = as.factor(paste(input_files[["PEP"]]$`Gene names`, input_files[["PEP"]]$Sequence, input_files[["PEP"]]$Modifications, sep="_")), 
                            "Sample" = as.factor(input_files[["PEP"]]$`Raw file`), 
                            "Intensity" = input_files[["PEP"]]$Intensity)
  psm_sig_raw <- dcast(psm_sig_raw, formula = ID_peptide~Sample, value.var = "Intensity", fun.aggregate = sum)
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
  psm_anno_raw <- data.table("ID_peptide"=psm_peptide_table$ID_peptide,
                             "symbol"=psm_peptide_table$GeneName,
                             "sequence"=psm_peptide_table$Annotated_Sequence,
                             "modifications"=psm_peptide_table$Modifications)
  
  psm_anno_raw$old_id<-psm_anno_raw$ID_peptide
  
  # create peptide names linked to symbols
  psm_anno_raw[, row := seq(1,nrow(psm_anno_raw))]
  psm_anno_raw <- psm_anno_raw %>% dplyr::group_by(symbol) %>% dplyr::mutate("min"=min(row),"card"=n()) %>% ungroup()
  psm_anno_raw$rank<- psm_anno_raw$row - psm_anno_raw$min +1
  psm_anno_raw$id<- paste(psm_anno_raw$symbol,psm_anno_raw$rank,psm_anno_raw$card,sep="_")
  psm_anno_raw <- as.data.table(psm_anno_raw)
  
  
  # Convert to data.table
  psm_sig_prot_raw <- as.data.table(psm_sig_raw)
  psm_sig_pet_raw <- as.data.table(psm_sig_raw)
  
  n_prot_preFILT <- uniqueN(psm_anno_raw$symbol)
  n_pep_pre_FILT <- nrow(psm_peptide_table)
  
  # Preprocess Protein intensities
  psm_sig_prot_raw[psm_sig_prot_raw == 0] <- NA  # Transform 0s into NAs
  sig_thr <- filt_absent_value  # NA threshold
  
  psm_long_dt <- melt(psm_sig_prot_raw, id.vars = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is.na(counts))), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(min_c <= sig_thr)), by = ID_peptide]
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_prot_df <- psm_sig_prot_raw[ID_peptide %in% filter_ID_peptide]
  psm_anno_df <- as.data.table(psm_anno_raw)[ID_peptide %in% filter_ID_peptide]
  
  # Log2 transformation
  psm_log_prot_df <- copy(psm_sig_prot_df)
  psm_log_prot_df[, (setdiff(names(psm_log_prot_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_prot_df), "ID_peptide")]
  
  # Filter proteins with only 1 peptID_peptidee
  # TODO: change this filter
  # filter_df_single_pep <- psm_anno_df[, .N, by = symbol][N > 1]
  # filter_df_single_pep <- merge(filter_df_single_pep, psm_anno_df, by = "symbol")[, ID_peptide]
  # psm_log_prot_df <- psm_log_prot_df[ID_peptide %in% filter_df_single_pep]
  # psm_anno_df <- psm_anno_df[ID_peptide %in% filter_df_single_pep]
  
  # Preprocess peptID_peptidee intensities
  psm_sig_pet_raw[psm_sig_pet_raw == 0] <- NA  # Transform 0s into NAs
  
  psm_long_dt <- melt(psm_sig_pet_raw, id = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is.na(counts))), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(min_c <= sig_thr)), by = ID_peptide]
  
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
  
  # Filter proteins with only 1 peptID_peptidee
  filter_df_single_pep <- psm_anno_df[, .N, by = symbol][N > 0]
  filter_df_single_pep <- merge(filter_df_single_pep, psm_anno_df, by = "symbol")[, symbol]
  filter_df_single_pep <- psm_peptide_table[GeneName %in% filter_df_single_pep, ID_peptide]
  
  psm_log_pet_df <- psm_log_pet_df[ID_peptide %in% filter_df_single_pep]
  psm_peptide_table <- psm_peptide_table[ID_peptide %in% filter_df_single_pep]
  
  message(paste0("\nN° Proteins (raw): ", n_prot_preFILT))
  message(paste0("N° Proteins (after filter): ", uniqueN(psm_anno_df$symbol)))
  message(paste0("\nN° Peptides (raw): ", n_pep_pre_FILT))
  message(paste0("N° Peptides (after filter): ", uniqueN(psm_peptide_table$ID_peptide)))
  
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
                          filt_absent_value = 0){
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
    stop("\'Condition\' column missin in \'ANNOTATION\' file.")
  }
  if(!(sample_col %in% colToKeep)){
    stop("\'Sample\' column missin in \'ANNOTATION\' file.")
  }
  if(batch_corr_exe & !(batch_col %in% colToKeep)){ 
    stop(paste0("\'",batch_col,"\' column missin in \'ANNOTATION\' file with batch correction activated"))
  }
  input_files[["annotation"]] <- input_files[["annotation"]][, ..colToKeep]
  setnames(input_files[["annotation"]], old = "File ID", new = "File_ID")
  setnames(input_files[["annotation"]], old = condition_col, new = "Condition")
  setnames(input_files[["annotation"]], old = sample_col, new = "Sample")
  
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
  
  # Maintain only the first UNIPROT code
  input_files$PEP[, `Master Protein Accessions` := sapply(`Master Protein Accessions`, function(x) strsplit(x, ";")[[1]][1])]
  
  # Verify and clean PEP table
  abundance_cols <- grep("Abundance", names(input_files$PEP), value = TRUE)
  input_files$PEP <- input_files$PEP[, c("Annotated Sequence", if ("Modifications" %in% names(input_files$PEP)) "Modifications" else "Modifications in Master Proteins", 
                                         "Master Protein Accessions", "Positions in Master Proteins", abundance_cols), with = FALSE]
  # Merge PROT and PEP tables
  input_files$PD_PEP_matrix <- merge.data.table(input_files$PROT, input_files$PEP, by.x = "Accession", by.y = "Master Protein Accessions")
  
  # Filter sample in file annotation
  colToKeep <- c(c("Accession","Description","GeneName","Annotated Sequence",
                   if("Modifications" %in% colnames(input_files$PEP)){"Modifications"}else{"Modifications in Master Proteins"},
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
  psm_anno_raw <- data.table("ID_peptide"=psm_peptide_table$ID_peptide,
                             "symbol"=psm_peptide_table$GeneName,
                             "sequence"=psm_peptide_table$Annotated_Sequence,
                             "modifications"=psm_peptide_table$Modifications)
  
  psm_anno_raw$old_id<-psm_anno_raw$ID_peptide
  
  # create peptide names linked to symbols
  psm_anno_raw[, row := seq(1,nrow(psm_anno_raw))]
  psm_anno_raw <- psm_anno_raw %>% dplyr::group_by(symbol) %>% dplyr::mutate("min"=min(row),"card"=n()) %>% ungroup()
  psm_anno_raw$rank<- psm_anno_raw$row - psm_anno_raw$min +1
  psm_anno_raw$id<- paste(psm_anno_raw$symbol,psm_anno_raw$rank,psm_anno_raw$card,sep="_")
  psm_anno_raw <- as.data.table(psm_anno_raw)
  
  
  # Convert to data.table
  psm_sig_prot_raw <- as.data.table(psm_sig_raw)
  psm_sig_pet_raw <- as.data.table(psm_sig_raw)
  
  n_prot_preFILT <- uniqueN(psm_anno_raw$symbol)
  n_pep_pre_FILT <- nrow(psm_peptide_table)
  
  # Preprocess Protein intensities
  psm_sig_prot_raw[psm_sig_prot_raw == 0] <- NA  # Transform 0s into NAs
  sig_thr <- filt_absent_value  # NA threshold
  
  psm_long_dt <- melt(psm_sig_prot_raw, id.vars = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is.na(counts))), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(min_c <= sig_thr)), by = ID_peptide]
  
  filter_ID_peptide <- psm_filter_dt[passes_c > 0, ID_peptide]
  psm_sig_prot_df <- psm_sig_prot_raw[ID_peptide %in% filter_ID_peptide]
  psm_anno_df <- as.data.table(psm_anno_raw)[ID_peptide %in% filter_ID_peptide]
  
  # Log2 transformation
  psm_log_prot_df <- copy(psm_sig_prot_df)
  psm_log_prot_df[, (setdiff(names(psm_log_prot_df), "ID_peptide")) := lapply(.SD, log2), .SDcols = setdiff(names(psm_log_prot_df), "ID_peptide")]
  
  # Filter proteins with only 1 peptID_peptidee
  # TODO: change this filter
  # filter_df_single_pep <- psm_anno_df[, .N, by = symbol][N > 1]
  # filter_df_single_pep <- merge(filter_df_single_pep, psm_anno_df, by = "symbol")[, ID_peptide]
  # psm_log_prot_df <- psm_log_prot_df[ID_peptide %in% filter_df_single_pep]
  # psm_anno_df <- psm_anno_df[ID_peptide %in% filter_df_single_pep]
  
  # Preprocess peptID_peptidee intensities
  psm_sig_pet_raw[psm_sig_pet_raw == 0] <- NA  # Transform 0s into NAs
  
  psm_long_dt <- melt(psm_sig_pet_raw, id = "ID_peptide", variable.name = "sample", value.name = "counts")
  psm_long_dt <- merge(psm_long_dt, as.data.table(c_anno), by = "sample", all.x = TRUE)
  
  psm_filter_dt <- psm_long_dt[, .(min_c = sum(is.na(counts))), by = .(ID_peptide, condition)]
  psm_filter_dt <- psm_filter_dt[, .(passes_c = sum(min_c <= sig_thr)), by = ID_peptide]
  
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
  
  # Filter proteins with only 1 peptID_peptidee
  filter_df_single_pep <- psm_anno_df[, .N, by = symbol][N > 0]
  filter_df_single_pep <- merge(filter_df_single_pep, psm_anno_df, by = "symbol")[, symbol]
  filter_df_single_pep <- psm_peptide_table[GeneName %in% filter_df_single_pep, ID_peptide]
  
  psm_log_pet_df <- psm_log_pet_df[ID_peptide %in% filter_df_single_pep]
  psm_peptide_table <- psm_peptide_table[ID_peptide %in% filter_df_single_pep]
  
  message(paste0("\nN° Proteins (raw): ", n_prot_preFILT))
  message(paste0("N° Proteins (after filter): ", uniqueN(psm_anno_df$symbol)))
  message(paste0("\nN° Peptides (raw): ", n_pep_pre_FILT))
  message(paste0("N° Peptides (after filter): ", uniqueN(psm_peptide_table$ID_peptide)))
  
  return(list("c_anno" = c_anno,
              "psm_anno_df" = psm_anno_df,
              "psm_log_prot_df" = psm_log_prot_df,
              "psm_peptide_table" = psm_peptide_table,
              "psm_log_pet_df" = psm_log_pet_df,
              "colour_vec" = colour_vec))
}
