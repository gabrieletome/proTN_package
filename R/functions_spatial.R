#' Parallel Moran's I Computation for Protein Spatial Data
#'
#' Computes Moran's I index of spatial autocorrelation for multiple proteins
#' in parallel. Supports multi-plate layouts, optional imputation of missing
#' values, and plotting of spatial maps.
#'
#' @param proteome_data A list containing proteome data, including gene intensity data and sample annotations.
#' @param geneNames A character vector of gene names to process.
#' @param nPixels Integer number of total pixels per gene.
#' @param nsim Number of Monte Carlo simulations to use for significance testing (default: 999).
#' @param impute Logical; if TRUE, imputes missing values using local means (3x3 kernel).
#' @param cutoff Minimum number of valid (non-NA) values required to compute Moran's I (default: 8).
#' @param ncores Numeric; Set the max available cores for multithreading (default: ncores-1)
#'
#' @return A `data.table` with columns:
#' \describe{
#'   \item{gene}{Gene name.}
#'   \item{moran}{Observed Moran's I statistic.}
#'   \item{p.value}{Monte Carlo p-value.}
#'   \item{valid.values}{Number of valid (non-NA) pixels.}
#' }
#' 
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom parallel detectCores makeCluster clusterExport parLapply stopCluster
#' @importFrom terra rast focal values
#' @importFrom raster rasterToPolygons
#' @importFrom spdep poly2nb nb2listw moran.mc
#' @importFrom viridis viridis
#' @importFrom dplyr filter
#' @importFrom data.table data.table rbindlist
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' results <- moran_index(protein.data, geneNames = c("TP53", "SMN1"), nPixels = 384)
#' }
moran_index <- function(proteome_data, geneNames = NULL, nPixels = NULL,
                        nsim = 999, impute = FALSE, cutoff = 8, ncores = NULL) {
  if(!("dat_gene" %in% names(proteome_data))){
    if(!("psm_log_prot_df" %in% names(proteome_data))){
      message("Normalized gene data.table not detected. Computing moran index on the raw values.")
      impute = TRUE
      protein.data <- proteome_data$psm_log_prot_df
    } else{
      stop("Data.table intensity not detected")
    }
  } else{
    protein.data <- proteome_data$dat_gene
    impute = FALSE
  }
  
  if(!("c_anno" %in% names(proteome_data))){
    stop("Missing annotation data.table")
  }
  c_anno <- proteome_data$c_anno
  
  if(is.null(geneNames)){
    message("Computing Moran Index for all genes, can took several minutes")
    geneNames <- unique(protein.data$GeneName)
  }
  
  if(is.null(nPixels)){
    message("Detecting number of pixel from annotation")
    nPixels <- nrow(c_anno)
  }
  
  # Helper function for imputation
  fill_na <- function(x, i = 5) {
    if (is.na(x[i])) round(mean(x, na.rm = TRUE), 0) else x[i]
  }
  
  # Detect and register parallel backend
  if(is.null(ncores)){
    available_cores <- max(1, detectCores() - 1)
  } else{
    available_cores <- min(ncores, max(1, detectCores() - 1))
  }
  message(paste0("Multicore using: ",available_cores," cores."))
  
  if(length(geneNames) >= 1000){
    
    
    cl <- makeCluster(available_cores)
    clusterEvalQ(cl, {
      library(terra)
      library(spdep)
      library(viridis)
      library(data.table)
      library(dplyr)
      library(raster)
    })
    clusterExport(cl, varlist = c("protein.data", "geneNames", "c_anno",
                                  "nsim", "impute", "cutoff", "plot", "fill_na"),
                  envir = environment())
    
    
    results_list <- parLapply(cl, geneNames, function(gene) {
      # results_list <- foreach(
      #   gene = geneNames,
      #   .packages = c("terra", "spdep", "viridis", "dplyr", "data.table")
      # ) %dopar% {
      
      vec <- melt(protein.data[GeneName == gene], id.vars = "GeneName", value.name = "Intensity", variable.name = "sample")
      
      if (all(is.na(vec$Intensity))) {
        return(data.table(gene, moran = NA_real_, p.value = NA_real_, valid.values = 0L))
      }
      
      vec <- merge.data.table(vec, c_anno, by = "sample")
      # Convert to matrix: rows = y, columns = x
      mat <- dcast(vec, y ~ x, value.var = "Intensity")
      intensity_matrix <- as.matrix(mat[, -"y"])
      rownames(intensity_matrix) <- mat$y
      rast_mat <- rast(intensity_matrix)
      
      # Impute missing values if requested
      if (impute) {
        rast_mat <- focal(rast_mat, w = matrix(1, 3, 3), fun = fill_na, pad = TRUE, na.rm = FALSE)
      }
      
      # Optional plotting
      # if (plot) {
      #   plot(
      #     as.matrix(rast_mat), xlab = "", ylab = "", xaxt = "n", na.col = "white",
      #     col = viridis(24),
      #     main = paste(gene, ifelse(impute, "(Imputed)", ""))
      #   )
      # }
      
      valid <- sum(!is.na(values(rast_mat)))
      if (valid < cutoff) {
        return(data.table(gene, moran = NA_real_, p.value = NA_real_, valid.values = valid))
      }
      
      # Compute neighborhood structure
      w <- tryCatch({
        poly2nb(rasterToPolygons(as(rast_mat, "Raster"), na.rm = FALSE))
      },
      error = function(e) {as.character(e)}
      )
      if (is.null(w)) {
        return(data.table(gene, moran = NA_real_, p.value = NA_real_, valid.values = valid))
      }
      
      # Compute Moran's I
      mi <- tryCatch({
        moran.mc(
          values(rast_mat),
          nb2listw(w, style = "B"),
          nsim = nsim, na.action = na.omit, zero.policy = TRUE
        )
      },
      error = function(e)  {as.character(e)}
      )
      
      if (is.null(mi)) {
        return(data.table(gene, moran = NA_real_, p.value = NA_real_, valid.values = valid))
      }
      
      data.table(gene, moran = mi$statistic, p.value = mi$p.value, valid.values = valid)
    })
    
    stopCluster(cl = cl)
    
  } else{
    registerDoParallel(cores = available_cores)
    
    results_list <- foreach(
      gene = geneNames,
      .packages = c("terra", "spdep", "viridis", "dplyr", "data.table")
    ) %dopar% {
      
      vec <- melt(protein.data[GeneName == gene], id.vars = "GeneName", value.name = "Intensity", variable.name = "sample")
      
      if (all(is.na(vec$Intensity))) {
        return(data.table(gene, moran = NA_real_, p.value = NA_real_, valid.values = 0L))
      }
      
      vec <- merge.data.table(vec, c_anno, by = "sample")
      # Convert to matrix: rows = y, columns = x
      mat <- dcast(vec, y ~ x, value.var = "Intensity")
      intensity_matrix <- as.matrix(mat[, -"y"])
      rownames(intensity_matrix) <- mat$y
      rast_mat <- rast(intensity_matrix)
      
      # Impute missing values if requested
      if (impute) {
        rast_mat <- focal(rast_mat, w = matrix(1, 3, 3), fun = fill_na, pad = TRUE, na.rm = FALSE)
      }
      
      # Optional plotting
      # if (plot) {
      #   plot(
      #     as.matrix(rast_mat), xlab = "", ylab = "", xaxt = "n", na.col = "white",
      #     col = viridis(24),
      #     main = paste(gene, ifelse(impute, "(Imputed)", ""))
      #   )
      # }
      
      valid <- sum(!is.na(values(rast_mat)))
      if (valid < cutoff) {
        return(data.table(gene, moran = NA_real_, p.value = NA_real_, valid.values = valid))
      }
      
      # Compute neighborhood structure
      w <- tryCatch({
        poly2nb(rasterToPolygons(as(rast_mat, "Raster"), na.rm = FALSE))
      },
      error = function(e) {as.character(e)}
      )
      if (is.null(w)) {
        return(data.table(gene, moran = NA_real_, p.value = NA_real_, valid.values = valid))
      }
      
      # Compute Moran's I
      mi <- tryCatch({
        moran.mc(
          values(rast_mat),
          nb2listw(w, style = "B"),
          nsim = nsim, na.action = na.omit, zero.policy = TRUE
        )
      },
      error = function(e)  {as.character(e)}
      )
      
      if (is.null(mi)) {
        return(data.table(gene, moran = NA_real_, p.value = NA_real_, valid.values = valid))
      }
      
      data.table(gene, moran = mi$statistic, p.value = mi$p.value, valid.values = valid)
    }
    
    stopImplicitCluster()
  }
  results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
  results$FDR <- p.adjust(results$p.value, method = "BH") # FDR (Benjamini-Hochberg)
  colnames(results) <- c("GeneName", "moran_index", "p_val", "valid_pixels", "p_adj")
  results <- results[order(p_val), ] # Results ordered by p_value
  
  return(results)
}


