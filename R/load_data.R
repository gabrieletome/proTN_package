#' From Proteome Discoverer or MaxQuant files to lists of data tables.
#'
#' This function reads proteomics data.
#'
#' @param software Must be \strong{PD} for Proteome Discoverer or \strong{MQ} for MaxQuant
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
                            batch_corr_exe = FALSE){
  
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
    proteome_data = read_PD_files(anno_filename, pep_filename, prot_filename, batch_corr_exe= batch_corr_exe)
  } else if(software == "MQ"){
    proteome_data = read_MQ_files(anno_filename, pep_filename, batch_corr_exe = batch_corr_exe)
  }
  
  return(proteome_data)
}


# Read MaxQuant files
# Output: psm_sig_raw, c_anno
read_MQ_files <- function(anno_filename, pep_filename, 
                          condition_col="Condition", sample_col="Sample", 
                          color_col="Color", batch_corr_exe = FALSE, batch_col="batch"){
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
  
  to_remove <- nrow(input_files[["PEP"]][is.na(`Gene names`)])
  input_files[["PEP"]] <- input_files[["PEP"]][!is.na(`Gene names`)]
  message(paste0("\tPeptide removed due to missing Gene Names: ",to_remove))
  
  to_remove <- nrow(input_files[["PEP"]][!(`Raw file` %in% input_files[["annotation"]]$Sample)])
  input_files[["PEP"]] <- input_files[["PEP"]][(`Raw file` %in% input_files[["annotation"]]$Sample)]
  message(paste0("\tPeptide removed since sample not in sample annotation: ",to_remove))
                                                
  to_remove <- nrow(input_files[["PEP"]][(Type == "MSMS")])
  input_files[["PEP"]] <- input_files[["PEP"]][(Type != "MSMS")]
  message(paste0("\tPeptide removed since sample not in sample annotation: ",to_remove))
     
  # Keep only first gene name
  input_files[["PEP"]][str_detect(`Gene names`,";"), `Gene names` := tstrsplit(`Gene names`, ";", keep = 1)]

  #Made the matrix
  psm_sig_raw <- data.table("ID_peptide" = as.factor(paste(input_files[["PEP"]]$`Leading razor protein`, input_files[["PEP"]]$Modifications, input_files[["PEP"]]$Sequence, sep="_")), 
                            "Sample" = as.factor(input_files[["PEP"]]$`Raw file`), 
                            "Intensity" = input_files[["PEP"]]$Intensity)
  psm_sig_raw <- dcast(psm_sig_raw, formula = ID_peptide~Sample, value.var = "Intensity", fun.aggregate = sum)
  colnames(psm_sig_raw)[-1] <- input_files[["annotation"]][match(colnames(psm_sig_raw)[-1], input_files[["annotation"]]$Sample)]$Sample
  
  #Made the description of the mpeptides
  psm_peptide_table <- as.data.table(unique(input_files[["PEP"]][, c("Leading razor protein", "Protein names", "Gene names", "Sequence", "Modifications")]))
  colnames(psm_peptide_table) <- c("Accession","Description","GeneName","Annotated_Sequence","Modifications")
  psm_peptide_table[, ID_peptide := paste(Accession, Modifications, Annotated_Sequence, sep="_")]
  
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
  
  # TODO: add chunk peptide annotation
  
  # TODO: add chunk preprocess_sig_raw
  
  
  return(list("annotation" = c_anno,
              "raw_dt" = psm_sig_raw,
              "peptide_annotation" = psm_peptide_table,
              "colour_vec" = colour_vec))
}
