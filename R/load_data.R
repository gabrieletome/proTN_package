#' Extract case study
#'
#' @param path_phospho Character; path where extract case study
#' @param path_proteome Character; path where extract case study
#' @examples
#' \dontrun{
#' extract_example(path_proteome = tempdir())
#' }
#' @import utils
#' @import stringr
#' @importFrom grDevices cairo_pdf dev.off pdf
#' @importFrom stats cmdscale dist fisher.test model.matrix na.omit p.adjust prcomp rnorm sd setNames
#' @export
extract_example = function(path_phospho = NULL, 
                           path_proteome = NULL) {
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
  if(is.null(path_phospho) & is.null(path_proteome)){
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
#' @param annotation_filename Character; the name or pattern of the annotation file. Default is "annotation".
#' @param peptide_filename Character; the name or pattern of the peptide file. Default is "pep".
#' @param proteinGroup_filename Character; the name or pattern of the protein group file. Default is "prot".
#' @param condition_col Character; the column name representing condition information in the annotation file. Default is "Condition".
#' @param sample_col Character; the column name representing sample information in the annotation file. Default is "Sample".
#' @param color_col Character; the column name representing color information for the plot in the annotation file. Default is "Color".
#' @param batch_corr_exe Logical; whether batch correction should be applied. Default is FALSE.
#' @param batch_col Character; the column name representing batch information in the annotation file. Default is "batch".
#' @param filt_absent_value Numeric; the value used to filter out absent data. Default is 0.
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
read_proteomics <- function(software, folder, annotation_filename = "annotation",
                            peptide_filename = "pep", proteinGroup_filename = "prot", 
                            condition_col="Condition", sample_col="Sample", 
                            color_col="Color", batch_corr_exe = FALSE, batch_col="batch",
                            filt_absent_value = 0){
  
  # source("./R/functions.R")
  
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
  
  message(paste0("N Proteins (after filter): ", uniqueN(psm_anno_df$symbol)))
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
  
  message(paste0("N Proteins (after filter): ", uniqueN(psm_anno_df$symbol)))
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
                                   filt_absent_value = 0){
  
  # source("./R/functions.R")
  
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
                                          batch_corr_exe= batch_corr_exe, filt_absent_value = filt_absent_value, phospho_thr = phospho_thr)
  } else if(software == "MQ"){
    proteome_data = read_phospho_MQ_files(anno_filename, pep_filename, keep_only_phosphomodification = keep_only_phosphomodification,
                                          batch_corr_exe = batch_corr_exe, filt_absent_value = filt_absent_value)
  }
  
  return(proteome_data)
}


# Read MaxQuant files
read_phospho_MQ_files <- function(anno_filename, pep_filename, keep_only_phosphomodification = T,
                                  phospho_thr = 0.75, condition_col="Condition", sample_col="Sample", 
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
  message(paste0("\tPeptide removed due to low phosphorilation level: ",to_remove))
  
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
  psm_sig_raw <- dcast(psm_sig_raw, formula = ID_peptide~Sample, value.var = "Intensity", fun.aggregate = sum)
  colnames(psm_sig_raw)[-1] <- input_files[["annotation"]][match(colnames(psm_sig_raw)[-1], input_files[["annotation"]]$Sample)]$Sample
  
  #Made the description of the mpeptides
  psm_peptide_table <- as.data.table(unique(input_files[["PEP"]][, c("Leading razor protein", "Protein names", "Gene names", "Sequence", "Modifications", "Phospho (STY) Probabilities")]))
  colnames(psm_peptide_table) <- c("Accession","Description","GeneName","Annotated_Sequence","Modifications","Phospho_%")
  psm_peptide_table[, ID_peptide := paste(GeneName, Annotated_Sequence, Modifications, sep="_")]
  
  #TODO: verify if is essential
  # {
  #   input_files$PEP_info <- psm_peptide_table
  #   # Set row names
  #   input_files$PEP_info[, rowname := paste(Accession,`Phospho_%`,Modifications, sep = "_")]
  #   
  #   # Extract phosphorylation pattern
  #   pattern <- data.table(Phospho = input_files$PEP_info$`Phospho_%`)
  #   
  #   # Split phosphorylation percentages into columns
  #   split_values <- tstrsplit(pattern$Phospho, "\\(|\\)", type.convert = TRUE)
  #   
  #   # Convert to data.table
  #   pattern <- data.table(do.call(cbind, split_values))
  #   setnames(pattern, as.character(seq_len(ncol(pattern))))
  #   
  #   # Add grouping columns
  #   pattern[, `:=`(group = input_files$PEP_info$ID_peptide,
  #                  seq = input_files$PEP_info$Annotated_Sequence)]
  #   
  #   # Convert even-indexed columns to numeric
  #   num_cols <- na.omit(names(pattern)[seq(2, ncol(pattern)-2, by = 2)])
  #   pattern[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
  #   
  #   # Aggregate means for numeric columns
  #   pattern2 <- pattern[, lapply(.SD, mean, na.rm = TRUE), by = group, .SDcols = num_cols]
  #   
  #   # Round numeric values
  #   pattern2[, (num_cols) := lapply(.SD, round, 3), .SDcols = num_cols]
  #   
  #   # Merge original values back
  #   setkey(pattern, group)
  #   setkey(pattern2, group)
  #   pattern2 <- pattern[pattern2]
  #   
  #   # Convert to character and replace NA values
  #   pattern2 <- pattern2[, lapply(.SD, as.character)]
  #   pattern2[is.na(pattern2)] <- "tmpC"
  #   
  #   # Format new phosphorylation sequence
  #   base_string <- paste0(rep("%s(%s)", ncol(pattern)-2), collapse = "")
  #   pattern2[, newseq := gsub("\\(tmpC\\)|tmpC| ", "", do.call(sprintf, c(fmt = base_string, .SD))), .SDcols = names(pattern2)[1:ncol(pattern)-2]]
  #   
  #   # Update PEP_info with new phosphorylation percentages
  #   input_files$PEP_info <- unique(input_files$PEP[, .(Accession = `Leading razor protein`, 
  #                                                      Description = `Protein names`, 
  #                                                      GeneName = `Gene names`, 
  #                                                      `Annotated Sequence` = Sequence, 
  #                                                      Modifications)])
  #   setkey(input_files$PEP_info, Accession, Modifications, `Annotated Sequence`)
  #   setkey(pattern2, group)
  #   
  #   input_files$PEP_info[pattern2, Phospho_% := newseq, on = .(Accession, Modifications, `Annotated Sequence`)]
  #   
  # }
  
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
  
  message(paste0("N Proteins (after filter): ", uniqueN(psm_anno_df$symbol)))
  message(paste0("N Peptides (after filter): ", uniqueN(psm_peptide_table$ID_peptide)))
  
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
  # Read protein file
  input_files[["PSM"]] <- tryCatch({
    as.data.table(read_xlsx(psm_file_filename))
  }, error=function(cond){
    stop(paste0("Missing file. The file \'PSM\' is missing or not have the pattern in the filename or there are duplicates files."))
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
  psm_anno_raw <- data.table("ID_peptide"=psm_peptide_table$ID_peptide,
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
  
  message(paste0("N Proteins (after filter): ", uniqueN(psm_anno_df$symbol)))
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
                                             filt_absent_value = 0){
  # source("./R/functions.R")
  
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
                                  batch_corr_exe= batch_corr_exe, filt_absent_value = filt_absent_value)
  } else if(software == "MQ"){
    proteome_data = read_MQ_files(anno_proteome_filename, pep_proteome_filename, 
                                  batch_corr_exe = batch_corr_exe, filt_absent_value = filt_absent_value)
  }
  message("Read Proteome data!")
  
  # Read phospho data
  phospho_data = NULL
  message("Reading Phosphoproteomic data...")
  if(software == "PD"){
    phospho_data = read_phospho_PD_files(anno_phospho_filename, pep_phospho_filename, prot_phospho_filename, psm_phospho_file_filename, 
                                          keep_only_phosphomodification = keep_only_phosphomodification,
                                          batch_corr_exe= batch_corr_exe, filt_absent_value = filt_absent_value, phospho_thr = phospho_thr)
  } else if(software == "MQ"){
    phospho_data = read_phospho_MQ_files(anno_phospho_filename, pep_phospho_filename, keep_only_phosphomodification = keep_only_phosphomodification,
                                          batch_corr_exe = batch_corr_exe, filt_absent_value = filt_absent_value)
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
