# PhosR merging file
# All function of PhosR package are merged
# in a single file.
# ProTN download the code of PhosR from github
# since the package in CRAN/Bioconductor does not
# permit the change of species.
# Link: https://github.com/PYangLab/PhosR


#' @title RUV for phosphoproteomics data normalisation
#'
#' @description This is a wrapper implementation of RUVIII for phosphoproteomics
#' data normalisation. This function will call
#' tailImpute function to impute all the missing values (if there is any) in the
#' phosphoproteomics data for
#' applying RUVIII. It will then return the normalised values for quantified
#' phosphosites and remove imputed values.
#'
#' @param mat a matrix (or PhosphoExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples.
#' @param M is the design matrix as defined in RUVIII.
#' @param ctl is the stable phosphosites (or negative controls as defined in
#' RUVIII).
#' @param k is the number of unwanted factors as defined in RUVIII.
#' @param m a numeric number for controlling mean downshifting.
#' @param s a numeric number for controlling standard deviation of downshifted
#' sampling values.
#' @param keepImpute a boolean to keep the missing value in the returned matrix.
#' @param ... additional parameters that may be passed to RUVIII.
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#'
#' @return A normalised matrix.
#'
#' 
#' @importFrom ruv RUVIII
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay
#'
#' @aliases RUVproteome
#'
#' 
RUVphospho <- function(mat, M, ctl, k = NULL, m = 1.6, s = 0.6,
                       keepImpute = FALSE, assay = NULL, ...) {
  if (missing(mat))
    stop("Parameter mat is missing!")
  if (missing(M))
    stop("Parameter M is missing!")
  if (missing(ctl))
    stop("Parameter ctl is missing!")
  
  pe = FALSE
  if (methods::is(mat, "PhosphoExperiment")) {
    pe = TRUE
    mat.orig = mat
    if (is.null(assay)) {
      mat = as.matrix(SummarizedExperiment::assay(mat))
    } else {
      mat = as.matrix(SummarizedExperiment::assay(mat, assay))
    }
  }
  mat = RUV(mat = mat, M = M, ctl = ctl, k = k, m = m, s = s,
            keepImpute = keepImpute, ...)
  if (pe) {
    SummarizedExperiment::assay(mat.orig,"normalised",withDimnames=FALSE) = 
      mat
    mat = mat.orig
  }
  mat
}


#' @importFrom ruv RUVIII
#' 
RUVproteome <- function(mat, M, ctl, k = NULL, m = 1.8, s = 0.3,
                        keepImpute = FALSE, assay = NULL, ...) {
  pe = FALSE
  if (methods::is(mat, "PhosphoExperiment")) {
    mat.orig = mat
    if (is.null(assay)) {
      mat = as.matrix(SummarizedExperiment::assay(mat))
    } else {
      mat = as.matrix(SummarizedExperiment::assay(mat, assay))
    }
    pe = TRUE
  }
  
  mat = RUV(mat = mat, M = M, ctl = ctl, k = k, m = m, s = s,
            keepImpute = keepImpute, ...)
  
  if (pe) {
    SummarizedExperiment::assay(mat.orig,"normalised",withDimnames=FALSE) = 
      mat
    mat = mat.orig
  }
  mat
}

# An internal wrapper function for calling RUVIII
RUV <- function(mat, M, ctl, k = NULL, m = m, s = s, keepImpute = keepImpute,
                ...) {
  mat.complete <- c()
  if (sum(is.na(mat)) > 0) {
    mat.complete <- t(tImpute(mat, m = m, s = s))
  } else {
    mat.complete <- t(mat)
  }
  
  mat.ruv <- t(RUVIII(Y = mat.complete, M = M, k = k, ctl = ctl,
                      return.info = FALSE, ...))
  
  if ((sum(is.na(mat)) > 0) & (keepImpute == FALSE)) {
    mat.ruv[which(is.na(mat))] <- NA
    return(mat.ruv)
  } else {
    return(mat.ruv)
  }
}

#' Generate set of stable phosphoporylated sites
#'
#' @usage getSPS(phosData, assays, conds, num)
#'
#' @param phosData a list of users' PhosphoExperiment objects from which 
#' generate SPSs
#' @param assays an assay to use for each dataset in phosData
#' @param conds a list of vector contains the conditions labels for each sample 
#' in the phosphoExperiment objects
#' @param num the number of identified SPSs, by default is 100
#'
#' @return A vectors of stably phosphorylated sites
#'
#' 
#' 

getSPS <-function (phosData, assays="Quantification", conds, num = 100) {
  if (missing(phosData)) 
    stop("phosData is missing")
  if (missing(conds)) 
    stop("conds is missing")
  
  
  sites <- sites.unique <- mat.max <- list()
  n <- length(phosData)
  m <- length(conds)
  if (n < 2) {
    stop("Please use more than one dataset")
  }
  if (n != m) {
    stop("Please use the same number of datasets and conditions")
  }
  for (i in seq(n)) {
    if (!"PhosphoExperiment" %in% is(phosData[[i]])) {
      stop("Wrong phosData, need to be a PhosphoExperiment object")
    }
  }
  for (i in seq(n)) {
    sites[[i]] <- paste(toupper(GeneSymbol(phosData[[i]])), 
                        paste(Residue(phosData[[i]]), 
                              Site(phosData[[i]]), sep = ""), 
                        sep = ";")
    
    sites.unique[[i]] <- unique(sites[[i]])
    nrep <- ncol(phosData[[i]])/length(unique(conds[[i]]))
    if (nrep == 1) {
      mat.mean <- SummarizedExperiment::assay(phosData[[i]], assays)
    } else {
      grps <- conds[[i]]
      mat.mean <- PhosR::meanAbundance(
        SummarizedExperiment::assay(phosData[[i]], assays),grps)
    }
    fun_val = (rep(0, ncol(mat.mean)))
    names(fun_val) = colnames(mat.mean)
    sites.mean <- 
      t(
        vapply(split(as.data.frame(mat.mean), sites[[i]]),
               colMeans, FUN.VALUE = fun_val
        )
      )
    sites.max <- apply(sites.mean, 1, function(x) {
      x[which.max(abs(x))]
    })
    mat.max[[i]] <- sort(abs(sites.max), decreasing = TRUE)
  }
  o <- as.data.frame(table(unlist(sites.unique)))
  if (length(which(o$Freq > 1)) < 200) {
    stop("Fewer than 200 overlapped sites")
  }
  if (length(which(o$Freq == m)) > 1000) {
    top <- as.character(o[which(o$Freq == m), 1])
  }
  else {
    message("Warning: there aren't enough overlappling sites")
    top <- as.character(o[which(o$Freq > 1), 1])
  }
  Ts <- data.frame(mat.max[[1]][top])
  for (i in seq(2,n,1)) {
    Ts <- cbind(Ts, mat.max[[i]][top])
  }
  Tc <- (apply(-abs(Ts), 2, rank) - 0.5)/nrow(Ts)
  Tt4 <- pchisq(-2 * rowSums(log(Tc)), (n - 1) * 2, lower.tail = FALSE)
  names(Tt4) <- top
  sites.sorted <- names(sort(Tt4, decreasing = TRUE))
  sites.sorted = paste0(sites.sorted, ";")
  
  return(sites.sorted[seq(num)])
}

#' Create frequency matrix
#'
#' @usage createFrequencyMat(substrates.seq)
#'
#' @param substrates.seq A substrate sequence
#'
#' @return A frequency matrix of amino acid from substrates.seq.
#'
#' 
createFrequencyMat <- function(substrates.seq) {
  # substrates.seq.split <-
  # sapply(substrates.seq, strsplit, '')
  
  substrates.seq.split = Map(function(substrates.seq) {
    unlist(strsplit(substrates.seq, split = ""))
  }, substrates.seq)
  
  frequency.mat <- matrix(data = 0, nrow = 20,
                          ncol = length(substrates.seq.split[[1]]))
  rownames(frequency.mat) <- c("A", "R",
                               "N", "D", "C", "E", "Q", "G", "H",
                               "I", "L", "K", "M", "F", "P", "S",
                               "T", "W", "Y", "V")
  colnames(frequency.mat) <- paste("p",
                                   seq_len(length(substrates.seq.split[[1]])),
                                   sep = "")
  
  # calculate frequency
  frequency.list = lapply(seq(ncol(frequency.mat)), function(i) {
    aa = mapply(function(x, i) x[i],
                substrates.seq.split, MoreArgs = list(i = i))
    tmp_col = frequency.mat[,i]
    for (j in seq_len(length(aa))) {
      if (aa[j] == "_") {
        next
      }
      tmp_col[aa[j]] <- tmp_col[aa[j]] + 1
    }
    tmp_col
  })
  frequency.mat = matrix(unlist(frequency.list), ncol = ncol(frequency.mat))
  colnames(frequency.mat) = paste("p",
                                  seq_len(length(substrates.seq.split[[1]])),
                                  sep = "")
  rownames(frequency.mat) <- c("A", "R",
                               "N", "D", "C", "E", "Q", "G", "H",
                               "I", "L", "K", "M", "F", "P", "S",
                               "T", "W", "Y", "V")
  
  frequency.mat <- frequency.mat/length(substrates.seq)
  return(frequency.mat)
}

#' Frequency scoring
#'
#' @usage frequencyScoring(sequence.list, frequency.mat)
#'
#' @param sequence.list A vector list of sequences
#' @param frequency.mat A matrix output from `createFrequencyMat`
#'
#' @return A vector of frequency score
#'
#' 
frequencyScoring <- function(sequence.list, frequency.mat) {
  
  frequency.score <- c()
  
  frequency.score = lapply(seq(length(sequence.list)), function(idx) {
    if (sequence.list[idx] == "") {
      sequence.list[idx] = "_"
    }
    seqs = unlist(mapply(strsplit, sequence.list[idx],
                         MoreArgs = list(split = "")))
    score <- 0
    if (is.na(sequence.list[idx])) {
      frequency.score <- c(frequency.score,
                           score)
      return()
    }
    for (i in seq_len(length(seqs))) {
      aa <- c("A", "R", "N", "D", "C",
              "E", "Q", "G", "H", "I",
              "L", "K", "M", "F", "P",
              "S", "T", "W", "Y", "V")
      if (!seqs[i] %in% aa) {
        next
      }
      score <- frequency.mat[seqs[i],
                             i] + score
    }
    frequency.score <- c(frequency.score,
                         score)
    frequency.score
  })
  frequency.score = unlist(frequency.score)
  
  if (!is.null(names(sequence.list))) {
    names(frequency.score) <- names(sequence.list)
  }
  
  return(frequency.score)
}
#' Plot signalome map
#'
#' @usage plotSignalomeMap(signalomes, color)
#'
#' @param signalomes output from `Signalomes` function
#' @param color a string specifying the color vector for kinases
#'
#' @return a ggplot object
#'
#'
#' @importFrom tidyr spread
#' @importFrom dplyr count
#' 
#' 
plotSignalomeMap <- function(signalomes, color) {
  
  df <- stack(signalomes$kinaseSubstrates)
  modules <- signalomes$proteinModule
  names(modules) <- unlist(lapply(
    strsplit(
      as.character(names(signalomes$proteinModules)), 
      ";"
    ), "[[", 1))
  df$cluster <- modules[df$values]
  
  df_balloon <- df
  df_balloon <- na.omit(df_balloon) %>% dplyr::count(.data$cluster, .data$ind)
  df_balloon$ind <- as.factor(df_balloon$ind)
  df_balloon$cluster <- as.factor(df_balloon$cluster)
  df_balloon <- tidyr::spread(df_balloon, .data$ind, .data$n)[,-1]
  df_balloon[is.na(df_balloon)] <- 0
  
  df_balloon <- do.call(rbind, lapply(seq(nrow(df_balloon)), function(x) {
    
    res <- unlist(lapply(df_balloon[x,], 
                         function(y) y/sum(df_balloon[x,])*100))
    
  }))
  
  df_balloon <- melt(as.matrix(df_balloon))
  colnames(df_balloon) <- c("cluster", "ind", "n")
  
  g <- ggplot2::ggplot(df_balloon, aes(x = .data$ind, y = .data$cluster)) + 
    geom_point(aes(col=.data$ind, size=.data$n)) + 
    scale_color_manual(values=color) + 
    scale_size_continuous(range = c(2, 17)) + 
    theme_classic() + 
    theme(
      aspect.ratio=0.25, 
      legend.position = "bottom",
      axis.line = element_blank(),            
      axis.title = element_blank(),   
      panel.grid.major.x = element_blank(),  
      panel.grid.minor.x = element_blank()) 
  g
  
}

#' Plot kinase network
#'
#' @usage plotKinaseNetwork(KSR, predMatrix, threshold = 0.9, color, 
#' type = NULL, verbose = FALSE)
#'
#' @param KSR Kinase-substrate relationship scoring results
#' @param predMatrix Output of kinaseSubstratePred function
#' @param threshold Threshold used to select interconnected kinases for
#'  the expanded signalomes
#' @param color A string specifying the color vector for nodes
#' @param type A type (\code{graph} or \code{chord}) of plot. If NULL, network 
#' graph is plotted
#' @param verbose Default to \code{TRUE} to show messages during the progress.
#' All messages will be suppressed if set to \code{FALSE}
#'  
#' @return a graphical plot
#' 
#' @importFrom network network
#' @importFrom stats cor
#' @importFrom GGally ggnet2
#' @import circlize
#' 
#' 
#' 
plotKinaseNetwork <- function(KSR, predMatrix, threshold = 0.9, color=NULL, 
                              type = NULL, verbose = FALSE) {
  
  if (is.null(type)) { type = "chord" } else {
    type <- match.arg(type, c("graph", "chord"),
                      several.ok = FALSE)
  }
  
  cor_mat <- stats::cor(KSR$combinedScoreMatrix)
  diag(cor_mat) <- 0
  
  cor_mat_values = cor_mat
  
  cor_mat <- apply(cor_mat, 2, function(x) x  > threshold)
  cor_mat[cor_mat == FALSE] <- 0
  cor_mat[cor_mat == TRUE] <- 1
  
  if (type == "graph") {
    if (is.null(color)) {
      stop("Parameter color cannot be NULL.")
    }
    if (verbose) {
      print("Generating network graph...")
    }
    
    n = network::network(cor_mat, directed=FALSE)
    network::set.edge.value(n, "cor", cor_mat_values^2)
    my_col = RColorBrewer::brewer.pal(11, "Spectral")
    
    g <- GGally::ggnet2(n, 
                        node.size=10, 
                        node.color=color, 
                        edge.size = "cor", 
                        size = "degree",
                        size.cut=3,
                        label=colnames(cor_mat),
                        label.size=2,
                        mode="circle",
                        label.color="black")
    print(g)
    
  }
  
  if (type == "chord") {
    
    if (verbose) {
      print("Generating circular graph...")
    }
    if (is.null(color)) {
      stop("Parameter color cannot be NULL.")
    }
    grid.col <- color
    names(grid.col) <- rownames(cor_mat)
    
    n = length(grid.col)
    
    circos.clear()
    circos.par(start.degree = 180)
    circos.initialize(factors = "a", xlim = c(0, n))
    chordDiagram(cor_mat, transparency = 0.2,
                 order = rownames(cor_mat),
                 grid.col = grid.col,
                 annotationTrack = c("name", "grid"), scale = TRUE)
    title("Kinase network")
    
  }
}

#' phosphosite/Gene set over-representation analysis
#'
#' @description This function performes phosphosite (or gene) set 
#' over-representation analysis using Fisher's exact test.
#'
#' @param geneSet an array of gene or phosphosite IDs (IDs are gene symbols etc
#' that match to your pathway annotation list).
#' @param annotation a list of pathways with each element containing an array of
#' gene or phosphosite IDs.
#' @param universe the universe/backgrond of all genes or phosphosites in your
#' profiled dataset.
#' @param alter test for enrichment ('greater', default), depletion ('less'), or
#' 'two.sided'.
#'
#' @return A matrix of pathways and their associated substrates and p-values.
#'
#'
#' 
pathwayOverrepresent <- function(geneSet, annotation, universe,
                                 alter = "greater") {
  if (missing(geneSet)) {
    stop("Parameter geneSet is missing!")
  }
  if (missing(annotation)) {
    stop("Parameter annotation is missing!")
  }
  if (missing(universe)) {
    stop("Parameter universe is missing!")
  }
  
  fisherTest.mat <- matrix("NA", ncol = 3, nrow = length(annotation))
  colnames(fisherTest.mat) <- c("pvalue", "# of substrates",
                                "substrates")
  for (i in seq_len(length(annotation))) {
    di <- length(intersect(geneSet, annotation[[i]]))
    dn <- length(setdiff(geneSet, annotation[[i]]))
    ndi <- length(setdiff(annotation[[i]], geneSet))
    ndn <- length(setdiff(universe, union(geneSet, annotation[[i]])))
    p <- stats::fisher.test(rbind(c(di, ndi), c(dn, ndn)),
                            alternative = alter)$p.value
    substrates <- paste(intersect(geneSet, annotation[[i]]),
                        collapse = "|")
    fisherTest.mat[i, ] <- c(p, di, substrates)
  }
  rownames(fisherTest.mat) <- names(annotation)
  fisherTest.mat <- fisherTest.mat[order(as.numeric(fisherTest.mat[,
                                                                   1])), ]
  return(fisherTest.mat)
}

#' Phosphosite/Gene set enrichment analysis
#'
#' This function performes phosphosite (or gene) set enrichment analysis using 
#' Wilcoxon Rank Sum test.
#'
#' @param geneStats an array of statistics (e.g. log2 FC) of all quantified
#' genes or phosphosite with names of the array as gene or phosphosite IDs.
#' @param annotation a list of pathways with each element containing an array of
#'  gene IDs.
#' @param alter test for enrichment ('greater', default), depletion ('less'), or
#'  'two.sided'.
#'
#' @return A matrix of pathways and their associated substrates and p-values.
#'
#'
#' 
pathwayRankBasedEnrichment <- function(geneStats, annotation,
                                       alter = "greater") {
  if (missing(geneStats))
    stop("Parameter geneStats is missing!")
  if (missing(annotation))
    stop("Parameter annotation is missing!")
  
  wilcoxTest.mat <- matrix("NA", ncol = 3, nrow = length(annotation))
  colnames(wilcoxTest.mat) <- c("pvalue", "# of substrates",
                                "substrates")
  
  # perform wilcox sum rank test for pathway enrichment
  # analysis
  for (i in seq_len(length(annotation))) {
    p.in <- intersect(names(geneStats), annotation[[i]])
    p.out <- setdiff(names(geneStats), annotation[[i]])
    
    if (length(geneStats[p.in]) <= 3) {
      wilcoxTest.mat[i, seq_len(3)] <- NA
    } else {
      wilcoxTest.mat[i, 1] <- stats::wilcox.test(geneStats[p.in],
                                                 geneStats[p.out], alternative = alter)$p.value
      wilcoxTest.mat[i, 2] <- length(p.in)
      wilcoxTest.mat[i, 3] <- paste(p.in, collapse = ";")
    }
  }
  
  rownames(wilcoxTest.mat) <- names(annotation)
  wilcoxTest.mat <- wilcoxTest.mat[order(as.numeric(wilcoxTest.mat[,
                                                                   1])), ]
  return(wilcoxTest.mat)
}

#' @title Select by treatment groups (replicate block)
#'
#' @description Select phosphosites that have been quantified in a given
#' percentage of treatment groups (e.g. 0.75 as 3 out of 4 replicates)
#' in n groups.
#'
#' @author Pengyi Yang, Taiyun Kim
#'
#' @usage selectGrps(mat, grps, percent, n, assay)
#'
#' @param mat a matrix (PhosphoExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples in replicates for different 
#' treatments. 
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#' @param grps a string specifying the grouping (replicates).
#' @param percent a percent from 0 to 1, specifying the percentage of quantified
#'  values in any treatment group.
#' @param n an integer indicating n or more replicates pass the percentage
#' filtering for a phosphosite to be included.
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#'
#' @return a filtered matrix (or a PhosphoExperiment Oject) with at least 
#' 'percent' quantification in one or more conditions. If an input \code{mat} is 
#' a SummarizedExperiment object, filtered SummarizedExperiment object will be 
#' returned.
#'
#'
selectGrps <- function(mat, grps, percent, n = 1, assay = NULL) {
  if (missing(mat))
    stop("Parameter mat is missing!")
  if (missing(grps))
    stop("Parameter grps is missing!")
  if (missing(percent))
    stop("Parameter percent is missing!")
  
  if (length(grps) != ncol(mat))
    stop("Length of vector grps must be equal to number of columns in mat")
  if ((percent < 0) || (percent > 1))
    stop("Parameter percent must be a numeric value between 0 and 1")
  
  mat.filtered = mat
  
  if (methods::is(mat, "PhosphoExperiment")) {
    if (is.null(assay)) {
      mat = SummarizedExperiment::assay(mat)
    } else {
      mat = SummarizedExperiment::assay(mat, assay)
    }
  }
  
  
  # split the matrix by groups and organise them to a list
  tmp <- lapply(split(seq_len(ncol(mat)), grps), function(i) mat[,
                                                                 i])
  
  test <- do.call(cbind, lapply(tmp, function(x) {
    rowSums(!is.na(x))/ncol(x) >= percent
  }))
  
  sel = rowSums(test) >= n
  
  return(mat.filtered[sel,])
}


#' @title selectTimes
#'
#' @usage selectTimes(mat, timepoint, order, percent, w, assay)
#'
#' @param mat a matrix (or PhosphoExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples in replicates for different 
#' treatments.
#' @param timepoint a timepoint as factor with a length equal to the
#' number of columns of mat.
#' @param order a vector specifying the order of timepoints.
#' @param percent a percent (decimal) from 0 to 1, to filter phosphosites with
#' with missing value larger than percent per timepoint.
#' @param w a timepoint window for selection of phosphosites to remove.
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#'
#' @return a filtered matrix. If param \code{mat} is a SummarizedExperiment 
#' object, a SummarizedExperiment object will be returned.
#'
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' 
#' 
#'
selectTimes <- function(mat, timepoint, order, percent, w = 1, assay = NULL) {
  if (missing(mat))
    stop("Parameter mat is missing!")
  if (missing(timepoint))
    stop("Parameter timepoint is missing!")
  if (missing(order))
    stop("Parameter order is missing!")
  if (missing(percent))
    stop("Parameter percent is missing!")
  if ((percent < 0) || (percent > 1))
    stop("Parameter percent must be a numeric value between 0 and 1")
  
  
  mat.orig = mat
  if (methods::is(mat, "PhosphoExperiment")) {
    if (is.null(assay)) {
      mat = SummarizedExperiment::assay(mat)
    } else {
      mat = SummarizedExperiment::assay(mat, assay)
    }
  }
  
  # split the matrix by groups and organise them to a list
  tmp <- lapply(split(seq_len(ncol(mat)), timepoint), function(i) mat[,
                                                                      i,drop = FALSE])
  
  test <- do.call(cbind, lapply(tmp, function(x) {
    rowSums(!is.na(x))/ncol(x) >= percent
  }))[, order]
  
  if ((w == 1) | (w == length(order))) {
    sel = rowSums(test) >= w
    mat.filtered <- mat.orig[sel, ]
  } else if (w > length(order)) {
    stop("w is greater than the number of time points.
Please try a smaller w.")
  } else {
    idx <- length(order) - w + 1
    sel = Reduce(union, lapply(seq_len(idx), function(i) {
      which(rowSums(test[, seq(i, i+w-1),drop = FALSE]) == w)
    }))
    mat.filtered <- mat.orig[sel, ]
  }
  
  
  return(mat.filtered)
}



#' @title Select phosphosite by percentage of quantification
#'
#' @description Select phosphosites that have been quantified in more than a
#' given percentage of samples
#'
#' @usage selectOverallPercent(mat, percent, n, assay)
#'
#' @param mat a matrix (or PhosphoExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples in replicates for different 
#' treatments.
#' @param percent a percent from 0 to 1, specifying the percentage of quantified
#' values in across all samples for retaining a phosphosite for subsequent
#' analysis.
#' @param n an integer indicating n or more quantified values required for
#' retaining a phosphosite for subsequent analysis.
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#'
#' @return a filtered matrix
#'
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' 
#' 
#'
selectOverallPercent <- function(mat, percent = NULL, n = NULL, assay = NULL) { 
  
  
  if (missing(mat))
    stop("Parameter mat is missing!")
  if ((!is.null(percent)) && ((percent < 0) || (percent > 1)))
    stop("Parameter percent must be a numeric value between 0 and 1")
  if (is.null(percent) & is.null(n))
    stop("specify either percentage of number of quantified values for a
given phosphosite to be retained.")
  
  mat.orig = mat
  if (methods::is(mat, "PhosphoExperiment")) {
    if (is.null(assay)) {
      mat = SummarizedExperiment::assay(mat)
    } else {
      mat = SummarizedExperiment::assay(mat, assay)
    }
  }
  
  if (!is.null(percent) & is.null(n)) {
    sel = rowSums(!is.na(mat))/ncol(mat) >= percent
    mat.filtered <- mat.orig[sel, ]
  } else if (is.null(percent) & !is.null(n)) {
    sel = rowSums(!is.na(mat)) >= n
    mat.filtered <- mat.orig[sel, ]
  } else {
    sel = (rowSums(!is.na(mat))/ncol(mat) >=
             percent) & (rowSums(!is.na(mat)) >= n)
    mat.filtered <- mat.orig[sel, ]
  }
  
  
  return(mat.filtered)
}


#' @title Select phosphosites by localisation score
#'
#' @description Select phosphosites with a localisation score higher than the 
#' pre-defined probability score (default score = 0.75)
#'
#' @usage selectLocalisedSites(mat, loc=NULL, prob = 0.75)
#'
#' @param mat a matrix (or PhosphoExperiment object) with rows corresponding to 
#' phosphosites and columns corresponding to samples in replicates for different 
#' treatments.
#' @param loc a vector of localisation scores
#' @param prob a percent from 0 to 1, specifying the localisation probability 
#' of quantified values in across all samples for retaining a phosphosite for 
#' subsequent analysis.
#'
#' @return a filtered matrix
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' 
#' 
#'
#'

selectLocalisedSites <- function(mat, loc=NULL, prob = 0.75) { 
  
  if (missing(mat))
    stop("Parameter mat is missing!")
  if ((!is.null(prob)) && ((prob < 0) || (prob > 1)))
    stop("Parameter prob must be a numeric value between 0 and 1")
  if ((is.null(loc)) && is.null(Localisation(mat)) && 
      is.null(Localisation(mat)) && sum(is.na(Localisation(mat))) > 0)
    stop("Some or all localisation scores are missing")
  if ((!is.null(loc)) && (length(loc) != nrow(mat)))
    stop("Length of loc should equal to the number of rows in mat")
  
  mat.orig = mat
  if (methods::is(mat, "PhosphoExperiment")) {
    loc = Localisation(mat)
  } else {
    if (is.null(loc))
      stop("Localisation probabilty scores missing")
  }
  
  sel = loc >= prob
  mat.filtered <- mat.orig[sel,]
  
  
  return(mat.filtered)
}

#' Site- and condition-specific (sc) impute
#'
#' Impute the missing values for a phosphosite across replicates within a single
#' condition (or treatment)
#' if there are n or more quantified values of that phosphosite in that
#' condition.
#'
#' @usage scImpute(mat, percent, grps, assay)
#'
#' @param mat a matrix (or PhosphoExperiment object) with rows correspond to 
#' phosphosites and columns correspond to replicates within a condition.
#' @param percent a percent from 0 to 1, specifying the percentage of quantified
#' values in any treatment group.
#' @param grps a string specifying the grouping (replciates).
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#'
#' @return An imputed matrix. If param \code{mat} is a PhosphoExperiment 
#' object, a PhosphoExperiment object will be returned.
#'
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' 
#' 
scImpute <- function(mat, percent, grps, assay = NULL) {
  if (missing(mat)) {
    stop("Parameter mat is missing!")
  }
  if (missing(percent)) {
    stop("Parameter percent is missing!")
  }
  if (missing(grps)) {
    stop("Parameter grps is missing!")
  }
  if (length(grps) != ncol(mat)) {
    stop("Length of vector grps must be equal to number of columns in mat")
  }
  if ((percent < 0) || (percent > 1)) {
    stop("Parameter percent must be a numeric value between 0 and 1")
  }
  
  pe = FALSE
  if (methods::is(mat, "PhosphoExperiment")) {
    pe = TRUE
    mat.orig = mat
    if (is.null(assay)) {
      mat = SummarizedExperiment::assay(mat)
    } else {
      mat = SummarizedExperiment::assay(mat, assay)
    }
  }
  tmp <- lapply(split(seq_len(ncol(mat)), grps), function(i) mat[,
                                                                 i])
  mat.imputed <- do.call(cbind, lapply(tmp, stImp, percent = percent))[,
                                                                       colnames(mat)]
  
  if (pe) {
    SummarizedExperiment::assay(mat.orig, "imputed") = mat.imputed
    mat.imputed = mat.orig
  }
  
  return(mat.imputed)
}


stImp <- function(mat, percent) {
  for (i in seq_len(nrow(mat))) {
    idx <- which(!is.na(mat[i, ]))
    if (length(idx)/ncol(mat) >= percent) {
      ms <- mean(mat[i, ], na.rm = TRUE)
      sds <- stats::sd(mat[i, ], na.rm = TRUE)
      nid <- which(is.na(mat[i, ]))
      mat[i, nid] <- stats::rnorm(length(nid), mean = ms,
                                  sd = sds)
    }
  }
  return(mat)
}

#' Tail-based impute
#'
#' Tail-based imputation approach as implemented in Perseus.
#'
#' @usage tImpute(mat, m, s, assay)
#'
#' @param mat a matrix (or PhosphoExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples.
#' @param m a numeric number for controlling mean downshifting.
#' @param s a numeric number for controlling standard deviation of downshifted
#' sampling values.
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#'
#' @return An imputed matrix. If param \code{mat} is a SummarizedExperiment 
#' object, a SummarizedExperiment object will be returned.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom methods is
#'
#' 
tImpute <- function(mat, m = 1.6, s = 0.6, assay = NULL) {
  if (missing(mat)) {
    stop("Paramter mat is missing!")
  }
  
  mat.orig = mat
  pe = FALSE
  if (methods::is(mat, "PhosphoExperiment")) {
    if (is.null(assay)) {
      mat = SummarizedExperiment::assay(mat)
    } else {
      mat = SummarizedExperiment::assay(mat, assay)
    }
    pe = TRUE
  }
  
  ms <- colMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 2, stats::sd, na.rm = TRUE)
  mat.impute <- mat
  for (i in seq_len(ncol(mat.impute))) {
    r <- stats::rnorm(n = sum(is.na(mat.impute[, i])), mean = (ms[i] -
                                                                 sds[i] * m), sd = (sds[i] * s))
    mat.impute[which(is.na(mat.impute[, i])), i] <- r
  }
  
  if (pe) {
    SummarizedExperiment::assay(mat.orig, "imputed") = mat.impute
    mat.impute = mat.orig
  }
  
  return(mat.impute)
}



#' Paired-tail (pt) based impute
#'
#' Impute the missing values for mat2 using tail imputation approach if mat1 has
#' more than percent1 (percentage) of quantified values
#' and mat2 has less than percent2 (percentage) quantified values,
#' and vice versa if paired is set to be true. That is if mat2 has percentage of
#' quantified values more than percent1 and mat1 has percentage quantified
#' values less than percent2.
#'
#' @usage ptImpute(
#'     mat1, 
#'     mat2, 
#'     percent1, 
#'     percent2, 
#'     m = 1.6, 
#'     s = 0.6, 
#'     paired = TRUE, 
#'     verbose = TRUE,
#'     assay
#' )
#'
#' @param mat1 a matrix (or PhosphoExperiment object) with rows correspond to
#'  phosphosites and columns correspond to replicates within treatment1.
#' @param mat2 a matrix (or PhosphoExperiment object) with rows correspond to
#'  phosphosites and columns correspond to replicates within treatment2.
#' @param percent1 a percent indicating minimum quantified percentages required
#' for considering for imputation.
#' @param percent2 a percent indicating minimum quantified percentages required
#' for considering for imputation.
#' @param m a numeric number of for controlling mean downshifting.
#' @param s a numeric number of for controlling standard deviation of
#' downshifted sampling values.
#' @param paired a flag indicating whether to impute for both treatment1 and
#' treatment2 (default) or treatment2 only (if paired=FALSE).
#' @param verbose Default to \code{TRUE} to show messages during the progress.
#' All messages will be suppressed if set to \code{FALSE}
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#' 
#' 
#' @return An imputed matrix
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assay rowData colData
#' @importFrom methods is
#' 
#' 
#'
ptImpute <- function(mat1, mat2, percent1, percent2, m = 1.6,
                     s = 0.6, paired = TRUE, verbose = TRUE, assay = NULL) {
  if (missing(mat1))
    stop("Paramter mat1 is missing!")
  if (missing(mat2))
    stop("Paramter mat2 is missing!")
  if (missing(percent1))
    stop("Paramter percent1 is missing!")
  if (missing(percent2))
    stop("Paramter percent2 is missing!")
  
  pe = FALSE
  mat1.orig = mat1
  mat2.orig = mat2
  
  if (methods::is(mat1, "PhosphoExperiment") && 
      methods::is(mat2, "PhosphoExperiment")) {
    pe = TRUE
    if (is.null(assay)) {
      mat1 = SummarizedExperiment::assay(mat1)
      mat2 = SummarizedExperiment::assay(mat2)
    } else {
      mat1 = SummarizedExperiment::assay(mat1, assay)
      mat2 = SummarizedExperiment::assay(mat2, assay)
    }
  }
  # matrix before imputation
  mat1.raw = mat1
  mat2.raw = mat2
  
  # impute for mat2
  idx1 <- which((rowSums(!is.na(mat1))/ncol(mat1)) >= percent1 &
                  (rowSums(!is.na(mat2))/ncol(mat2)) <= percent2)
  
  if (verbose)
    message(paste("idx1:", length(idx1)))
  
  ms <- colMeans(mat2, na.rm = TRUE)
  sds <- apply(mat2, 2, stats::sd, na.rm = TRUE)
  if (length(idx1) > 0) {
    for (i in seq_len(ncol(mat2))) {
      mat2[idx1, i] <- stats::rnorm(length(idx1), mean = (ms[i] -
                                                            sds[i] * m), sd = (sds[i] * s))
    }
  }
  
  if (paired == TRUE) {
    # impute for mat1
    # If paired = TRUE, estimate idx2 before imputing mat2
    idx2 <- which((rowSums(!is.na(mat2.raw))/ncol(mat2.raw)) >= percent1 &
                    (rowSums(!is.na(mat1.raw))/ncol(mat1.raw)) <= percent2)
    if (verbose)
      message(paste("idx2:", length(idx2)))
    
    ms <- colMeans(mat1, na.rm = TRUE)
    sds <- apply(mat1, 2, stats::sd, na.rm = TRUE)
    if (length(idx2) > 0) {
      for (i in seq_len(ncol(mat1))) {
        mat1[idx2, i] <- stats::rnorm(length(idx2), mean = (ms[i] -
                                                              sds[i] * m), sd = (sds[i] * s))
      }
    }
    if (pe) {
      mat = cbind(mat1.orig, mat2.orig)
      if (is.null(assay)) {
        SummarizedExperiment::assay(mat, withDimnames = FALSE) = 
          cbind(mat1, mat2)
      } else {
        SummarizedExperiment::assay(mat, assay, withDimnames = FALSE) = 
          cbind(mat1, mat2)
      }
      
    } else {
      mat = cbind(mat1, mat2)
    }
    
    return(mat)
  } else {
    if (pe) {
      if (is.null(assay)) {
        SummarizedExperiment::assay(mat2.orig, withDimnames = FALSE) = 
          mat2
      } else {
        SummarizedExperiment::assay(mat2.orig, assay, 
                                    withDimnames = FALSE) = mat2            
      }
      mat2 = mat2.orig
    }
    return(mat2)
  }
}

#' @title PhosR Signalomes
#'
#' @description A function to generate signalomes
#'
#' @usage Signalomes(KSR, predMatrix, exprsMat, KOI, threskinaseNetwork=0.9,
#' signalomeCutoff=0.5, module_res = NULL, filter = FALSE, verbose = TRUE)
#'
#' @param KSR kinase-substrate relationship scoring results
#' @param predMatrix output of kinaseSubstratePred function
#' @param exprsMat a matrix with rows corresponding to phosphosites and columns
#' corresponding to samples
#' @param KOI a character vector that contains kinases of interest for which
#' expanded signalomes will be generated
#' @param threskinaseNetwork threshold used to select interconnected kinases for
#'  the expanded signalomes
#' @param signalomeCutoff threshold used to filter kinase-substrate
#' relationships
#' @param module_res parameter to select number of final modules
#' @param filter parameter to filter modules with only few proteins
#' @param verbose Default to \code{TRUE} to show messages during the progress.
#' All messages will be suppressed if set to \code{FALSE}
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom grDevices colorRampPalette
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr count
#' @importFrom dplyr n
#' @importFrom graphics title
#' @importFrom graphics par
#' @importFrom graphics barplot
#' @importFrom stats hclust
#' @importFrom rlang .data
#'
#' @return A list of 3 elements.
#'  \code{Signalomes}, \code{proteinModules} and \code{kinaseSubstrates}
#'
#' 
#' 

Signalomes <- function(KSR, 
                       predMatrix, 
                       exprsMat, 
                       KOI, 
                       threskinaseNetwork = 0.9,
                       signalomeCutoff = 0.5, 
                       module_res = NULL, 
                       filter = FALSE, 
                       verbose = TRUE) {
  
  if (!is.null(module_res)) {
    if (module_res < 20) {
      module_res = as.integer(module_res)
    } else {
      stop("module resolution should be an integer lower than 20")
    }
  } 
  
  ############## generate objects required for signalome function
  protein_assignment = mapply("[[",
                              strsplit(rownames(KSR$combinedScoreMatrix),";"),
                              MoreArgs = list(1))
  # KinaseFamily = PhosR::KinaseFamily
  utils::data("KinaseFamily", envir = environment())
  kinaseGroup <- KinaseFamily[, "kinase_group"]
  names(kinaseGroup) <- KinaseFamily[, "gene_symbol"]
  
  ############## set color palette
  my_color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,
                                                                           "Accent"))
  kinase_all_color <- my_color_palette(ncol(KSR$combinedScoreMatrix))
  names(kinase_all_color) <- colnames(KSR$combinedScoreMatrix)
  kinase_signalome_color <- kinase_all_color[colnames(predMatrix)]
  my_color_palette_kinaseGroup <-
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "Set2"))
  kinaseGroup_color <-
    my_color_palette_kinaseGroup(length(unique(kinaseGroup)))
  names(kinaseGroup_color) <- unique(kinaseGroup)
  
  ############## interconnected kinases
  resKinaseNetwork <- .kinaseNetwork(predMatrix, KSR, threskinaseNetwork,
                                     kinase_signalome_color)
  
  ############## cluster phosphosites
  substrate_clusters <- .phosphositeClusters(KSR, verbose)
  
  ############## generate coassignment
  cluster_assignment <- as.factor(substrate_clusters)
  dat.long <- data.frame(table(cluster_assignment, protein_assignment))
  dftoHeatmap <- tidyr::pivot_wider(dat.long, names_from = 
                                      .data$protein_assignment, values_from = .data$Freq)[,-1]
  dftoHeatmap[is.na(dftoHeatmap)] <- 0
  dftoHeatmap[dftoHeatmap > 0] <- 1
  hclust_res <- stats::hclust(stats::dist(t(dftoHeatmap)),
                              method = "ward.D")
  tree_height <- as.numeric(names(table(hclust_res$height)))
  branching <- as.numeric(table(hclust_res$height))
  
  tree_height_calc = unlist(lapply(seq(2,length(tree_height), 1), 
                                   function(x) {
                                     h <- tree_height[[x]]
                                     m <- stats::cutree(hclust_res, h = h)
                                     return(length(table(m)))
                                   }))
  
  if (!is.null(module_res)) {
    hcutree = which(tree_height_calc <= module_res) + 1
  } else {
    hcutree <- min(tree_height[tree_height > 0])
  }
  modules <- stats::cutree(hclust_res, h = tree_height[[hcutree[[1]]]])
  
  if (filter) {
    filter_modules = modules %in% which(table(modules) < 10)
    modules[filter_modules] = "noModule"
  } 
  
  ############## generate signalomes
  signalomeSubstrates <- .phosRsignalome(predMatrix, signalomeCutoff,
                                         kinase_signalome_color, modules)
  ############## generate kinase-specific signalomes (extended signalome)
  signalomes_of_KOI <- .getSignalomes(predMatrix, exprsMat,
                                      KOI, resKinaseNetwork, signalomeSubstrates, modules,
                                      kinaseGroup, verbose = verbose)
  signalome_res <- list(Signalomes = signalomes_of_KOI,
                        proteinModules = modules,
                        kinaseSubstrates = signalomeSubstrates)
  return(signalome_res)
}

.kinaseNetwork <- function(predMatrix, KSR, threskinaseNetwork,
                           kinase_signalome_color) {
  
  kinase_cor <- stats::cor(KSR$combinedScoreMatrix)
  
  cor_kinase_mat <- kinase_cor
  diag(cor_kinase_mat) <- 0
  kinase_network <- lapply(seq_len(ncol(cor_kinase_mat)), function(x)
    names(which(cor_kinase_mat[, x] > threskinaseNetwork)))
  names(kinase_network) <- colnames(cor_kinase_mat)
  
  cor_kinase_mat <- apply(cor_kinase_mat, 2,
                          function(x) x > threskinaseNetwork)
  cor_kinase_mat[cor_kinase_mat == FALSE] <- 0
  cor_kinase_mat[cor_kinase_mat == TRUE] <- 1
  
  network <- igraph::graph_from_adjacency_matrix(cor_kinase_mat,
                                                 mode = "undirected", diag = FALSE)
  
  kinaseNetwork.res <- list(kinaseNetwork = kinase_network,
                            kinaseCor = cor_kinase_mat)
  
  return(kinaseNetwork.res)
  
}

.phosphositeClusters <- function(KSR, verbose = TRUE) {
  substrate_cor <- stats::cor(t(KSR$combinedScoreMatrix))
  substrate_hclust <- stats::hclust(stats::dist(KSR$combinedScoreMatrix),
                                    method = "ward.D")
  if (verbose)
    message("calculating optimal number of clusters...")
  res <- lapply(seq(2,10,1), function(x) {
    substrate_clusters <- stats::cutree(substrate_hclust, k = x)
    cor.res <- lapply(seq_len(x), function(y) {
      substrate_cor = substrate_cor[substrate_clusters == y,
                                    substrate_clusters == y]
      diag(substrate_cor) <- 0
      return(substrate_cor)
    })
    cor.res <- lapply(cor.res, function(x) median(x))
    cor.res <- unlist(cor.res)
    cor.logic <- sum(cor.res >= 0.5) == x
    if (cor.logic) {
      cluster = mean(cor.res)
      names(cluster) = x
      return(cluster)
    }
  })
  if (isTRUE(is.null(unlist(res)))) {
    res <- lapply(seq(2,10,1), function(x) {
      substrate_clusters <- cutree(substrate_hclust, k = x)
      cor.res <- lapply(seq_len(x), function(y) {
        substrate_cor = substrate_cor[substrate_clusters == y,
                                      substrate_clusters == y]
        diag(substrate_cor) <- 0
        return(substrate_cor)
      })
      cor.res <- unlist(lapply(cor.res, function(x) median(x)))
      cor.logic <- sum(cor.res >= 0.1) == x
      if (cor.logic) {
        cluster = median(cor.res)
        names(cluster) = x
        return(cluster)
      }
    })
    res <- as.numeric(names(which(unlist(res) == max(unlist(res)))[1]))
  } else {
    res <- as.numeric(names(which(unlist(res) == max(unlist(res)))[1]))
  }
  if (verbose)
    message(paste0("optimal number of clusters = ", res))
  substrate_clusters <- stats::cutree(substrate_hclust, k = res)
  return(substrate_clusters)
}

#' @import circlize
#' @importFrom rlang .data
.phosRsignalome <- function(predMatrix, signalomeCutoff, kinase_signalome_color,
                            modules) {
  
  signalomeKinase <- colnames(predMatrix)
  
  signalomeSubstrates <- list()
  for (i in seq_len(length(signalomeKinase))) {
    signalomeSubstrates[[i]] =
      mapply(function(x) x[1],
             strsplit(names(which(
               predMatrix[,signalomeKinase[[i]]] > signalomeCutoff)),
               ";"))
  }
  names(signalomeSubstrates) <- signalomeKinase
  
  ############## generate circlize plot
  dftoPlot_signalome <- stack(signalomeSubstrates)
  dftoPlot_signalome$modules <- modules[dftoPlot_signalome$values]
  #adjacencyData <- with(dftoPlot_signalome, table(ind, modules))
  d = table(dftoPlot_signalome$ind, dftoPlot_signalome$modules)
  adjacencyData <- matrix(table(dftoPlot_signalome$ind, 
                                dftoPlot_signalome$modules), nrow = nrow(d), ncol = ncol(d))
  rownames(adjacencyData) = rownames(d)
  colnames(adjacencyData) = colnames(d)
  
  m = sort(as.integer(unique(dftoPlot_signalome$modules)))
  grid.col <- c(kinase_signalome_color, rep("grey", length(unique(m))))
  names(grid.col) <- c(rownames(adjacencyData), as.character(unique(m)))
  adjacencyData = adjacencyData[,!grepl("noModule",colnames(adjacencyData))]
  
  n = length(grid.col)
  circos.clear()
  circos.par(start.degree = 180)
  circos.initialize(factors = "a", xlim = c(0, n))
  chordDiagram(adjacencyData, transparency = 0.2,
               order = c(rownames(adjacencyData),
                         rev(unique(m))),
               grid.col = grid.col, #big.gap = 15,
               annotationTrack = c("name", "grid"), scale = TRUE)
  title("Signalomes")
  
  return(signalomeSubstrates)
}

#' @importFrom dplyr count %>%
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
.getSignalomes <- function(predMatrix, exprsMat, KOI, resKinaseNetwork,
                           signalomeSubstrates, modules, kinaseGroup, verbose = TRUE) {
  
  protein = mapply("[[", strsplit(rownames(exprsMat), ";"),
                   MoreArgs = list(1))
  KinaseSubstrateList <- resKinaseNetwork$kinaseNetwork
  
  ############## annotate kinase-substrate relationship
  kinaseAnnot = annoKinaseSubstrateRelation(predMatrix)
  
  ############## proportion of regulation by kinase
  signalomeMatrix <- stack(signalomeSubstrates)
  signalomeMatrix$cluster <- modules[signalomeMatrix$values]
  
  balloon_bycluster <- signalomeMatrix
  balloon_bycluster <- na.omit(balloon_bycluster) %>%
    dplyr::count(.data$cluster, .data$ind)
  balloon_bycluster$ind <- as.factor(balloon_bycluster$ind)
  balloon_bycluster$cluster <- as.factor(balloon_bycluster$cluster)
  balloon_bycluster <- as.data.frame(
    tidyr::pivot_wider(balloon_bycluster, 
                       names_from = .data$ind, values_from = .data$n))
  rownames(balloon_bycluster) = balloon_bycluster$cluster
  balloon_bycluster = balloon_bycluster[,-1]
  
  balloon_bycluster[is.na(balloon_bycluster)] <- 0
  balloon <- do.call(rbind, lapply(seq_len(nrow(balloon_bycluster)),
                                   function(x) {
                                     res = mapply(function(y, balloon_bycluster, x) {
                                       y/sum(balloon_bycluster[x, ]) * 100
                                     }, balloon_bycluster[x, ],
                                     MoreArgs = list(balloon_bycluster = balloon_bycluster, x = x))
                                   }))
  rownames(balloon) = rownames(balloon_bycluster)
  colnames(balloon) = colnames(balloon_bycluster)
  kinaseProportions <- round(balloon, 3)
  
  ############## generate kinase specific signalomes
  res = generateSignalome(kinaseAnnot, kinaseGroup, predMatrix, KOI,
                          KinaseSubstrateList, kinaseProportions,
                          signalomeSubstrates, exprsMat, protein, modules, 
                          verbose = verbose)
  return(res)
}

annoKinaseSubstrateRelation = function(predMatrix) {
  kinaseAnnot <- lapply(seq_len(nrow(predMatrix)), function(x) {
    scores <- predMatrix[x, ]
    siteName <- rownames(predMatrix)[[x]]
    kinaseTop <- names(which(scores == max(scores)))
    score <- as.numeric(max(scores))
    res <- c(siteName, kinaseTop, score)
    return(res)
  })
  kinaseAnnot <- do.call(rbind, kinaseAnnot)
  rownames(kinaseAnnot) <- kinaseAnnot[, 1]
  kinaseAnnot <- data.frame(kinaseAnnot[, -c(1)], stringsAsFactors = FALSE)
  colnames(kinaseAnnot) <- c("kinase", "score")
  kinaseAnnot$score <- as.numeric(kinaseAnnot$score)
  kinaseAnnot
}

generateSignalome = function(kinaseAnnot, kinaseGroup, predMatrix, KOI,
                             KinaseSubstrateList, kinaseProportions,
                             signalomeSubstrates, exprsMat, protein, modules,
                             verbose = TRUE) {
  annotation <- data.frame(kinase = as.factor(kinaseAnnot$kinase),
                           kinaseFamily = as.factor(kinaseGroup[kinaseAnnot$kinase]),
                           score = as.numeric(kinaseAnnot$score))
  rownames(annotation) <- rownames(predMatrix)
  
  kinaseProp = kinaseProportions[!grepl("noModule",
                                        rownames(kinaseProportions)),]
  m = modules[!grepl("noModule", modules)]
  
  res <- lapply(KOI, function(x) {
    if (x %in% names(KinaseSubstrateList)) {
      regModule <- which(kinaseProp[, x] > 1)
      kinases <- unique(c(x, KinaseSubstrateList[[x]]))
      substrates <- unique(unlist(
        lapply(kinases, function(x) signalomeSubstrates[[x]])))
      exprs_dat <- lapply(regModule, function(y) {
        exprsMat[protein %in%
                   names(m[m == y]) &
                   protein %in% substrates, ]
      })
      exprs_dat <- do.call(rbind, exprs_dat)
      annotation_dat <- annotation[rownames(annotation) %in%
                                     rownames(exprs_dat), ]
      kinaseSignalome <- list(exprs = exprs_dat,
                              annotation = annotation_dat)
      return(kinaseSignalome)
    } else {
      if (verbose)
        message(paste0(x, " is not found"))
    }
  })
  names(res) <- KOI
  res
}

#' Median centering and scaling
#'
#' Median centering and scaling of an input numeric matrix
#'
#'
#' @param mat a matrix with rows correspond to phosphosites and columns
#' correspond to samples.
#' @param scale a boolean flag indicating whether to scale the samples.
#' @param grps a string or factor specifying the grouping (replciates).
#' @param reorder To reorder the columns by group (\code{grps}).
#' By default (\code{reorder=FALSE}), original column order is maintained.
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#'
#' @return A median scaled matrix
#'
#'
medianScaling <- function(mat, scale = FALSE, grps = NULL, reorder = FALSE, 
                          assay = NULL) {
  pe = FALSE
  if (methods::is(mat, "PhosphoExperiment")) {
    pe = TRUE
    mat.orig = mat
    if (is.null(assay)) {
      mat = SummarizedExperiment::assay(mat)
    } else {
      mat = SummarizedExperiment::assay(mat, assay)
    }
  }
  mat.medianScaled <- NULL
  if (!is.null(grps)) {
    tmp <- lapply(split(seq_len(ncol(mat)), grps), function(i) mat[,
                                                                   i])
    mat.medianScaled <- do.call(cbind, lapply(tmp, medianScale,
                                              scale = scale))
    
    if (!reorder) {
      mat.medianScaled <- mat.medianScaled[, colnames(mat)]
    }
  } else {
    mat.medianScaled <- medianScale(mat, scale = scale)
  }
  
  if (pe) {
    SummarizedExperiment::assay(mat.orig, "scaled") = mat.medianScaled
    mat.medianScaled = mat.orig
  }
  
  return(mat.medianScaled)
}

medianScale <- function(mat, scale) {
  if (scale) {
    normcont <- stats::median(apply(mat, 2, stats::median,
                                    na.rm = TRUE))
    adjval <- apply(mat, 2, stats::median, na.rm = TRUE) -
      normcont
    mat.scaled <- normalizeMedianAbsValues(sweep(mat, 2,
                                                 adjval, "-"))
    return(mat.scaled)
  } else {
    normcont <- stats::median(apply(mat, 2, stats::median,
                                    na.rm = TRUE))
    adjval <- apply(mat, 2, stats::median, na.rm = TRUE) -
      normcont
    mat.unscaled <- sweep(mat, 2, adjval, "-")
    return(mat.unscaled)
  }
}


#' Standardisation
#'
#' Standardisation by z-score transformation.
#'
#' @usage standardise(mat)
#'
#' @param mat a matrix (or a PhosphoExperiment object) with rows correspond 
#' to phosphosites and columns correspond to samples.
#'
#' @return A standardised matrix
#'
#'
#' 
#'
standardise <- function(mat) {
  
  pe = FALSE
  if (methods::is(mat, "PhosphoExperiment")) {
    pe = TRUE
    mat.orig = mat
    if (is.null(assay)) {
      mat = SummarizedExperiment::assay(mat)
    } else {
      mat = SummarizedExperiment::assay(mat, assay)
    }
  }
  
  means <- apply(mat, 1, mean)
  sds <- apply(mat, 1, stats::sd)
  
  X2 <- sweep(mat, MARGIN = 1, STATS = means, FUN = "-")
  zX <- sweep(X2, MARGIN = 1, STATS = sds, FUN = "/")
  
  if (pe) {
    SummarizedExperiment::assay(mat.orig, "standardise") = zX
    zX = mat.orig
  }
  
  return(zX)
}

#' @title Multi-intersection, union
#'
#' @description A recusive loop for intersecting multiple sets.
#'
#' @aliases mUnion
#'
#' @usage
#' mIntersect(x, y, ...)
#' mUnion(x, y, ...)
#'
#' @param x,y,... objects to find intersection/union.
#'
#' @return An intersection/union of input parameters
#'
#' 
#'
mIntersect <- function(x, y, ...) {
  if (missing(...))
    intersect(x, y) else intersect(x, mIntersect(y, ...))
}


#' 
#'
mUnion <- function(x, y, ...) {
  if (missing(...))
    union(x, y) else union(x, mUnion(y, ...))
}



#' Minmax scaling
#'
#' Perform a minmax standardisation to scale data into 0 to 1 range
#'
#' @usage minmax(mat)
#'
#' @param mat a matrix with rows correspond to phosphosites and columns
#' correspond to condition
#'
#' @return Minmax standardised matrix
#'
#' 
#'
minmax <- function(mat) {
  # minmax normalise
  apply(mat, 2, function(x) {
    (x - min(x))/(max(x) - min(x))
  })
}

#' Summarising phosphosites to proteins
#'
#' Summarising phosphosite-level information to proteins for performing
#' downstream gene-centric analyses.
#'
#' @usage phosCollapse(mat, id, stat, by='min')
#'
#' @param mat a matrix with rows correspond to phosphosites and columns
#' correspond to samples.
#' @param id an array indicating the groupping of phosphosites etc.
#' @param stat an array containing statistics of phosphosite such as
#' phosphorylation levels.
#' @param by how to summarise phosphosites using their statistics. Either by
#' 'min' (default), 'max', or 'mid'.
#'
#' @return A matrix summarised to protein level
#'
#' 
#'
phosCollapse <- function(mat, id, stat, by = "min") {
  listMat <- split(data.frame(mat, stat), id)
  
  matNew <- c()
  if (by == "min") {
    matNew <- as.matrix(do.call(rbind, lapply(listMat, function(x) {
      if (nrow(x) == 1) {
        as.numeric(x[1, seq_len(ncol(x) - 1)])
      } else {
        x[order(as.numeric(x[, ncol(x)]))[1], seq_len(ncol(x) - 1)]
      }
    })))
  } else if (by == "max") {
    matNew <- as.matrix(do.call(rbind, lapply(listMat, function(x) {
      if (nrow(x) == 1) {
        as.numeric(x[1, seq_len(ncol(x) - 1)])
      } else {
        x[order(as.numeric(x[, ncol(x)]), decreasing = TRUE)[1],
          seq_len(ncol(x) - 1)]
      }
    })))
  } else if (by == "mid") {
    matNew <- as.matrix(do.call(rbind, lapply(listMat, function(x) {
      if (nrow(x) == 1) {
        as.numeric(x[1, seq_len(ncol(x) - 1)])
      } else {
        mid <- round(nrow(x)/2)
        x[order(as.numeric(x[, ncol(x)]))[mid], seq_len(ncol(x) - 1)]
      }
    })))
  } else {
    stop("Unrecognised way of collapsing the data!")
  }
  
  return(matNew)
}

#' ANOVA test
#'
#' @description Performs an ANOVA test and returns its adjusted p-value
#'
#'
#' @usage matANOVA(mat, grps)
#'
#' @param mat An p by n matrix where p is the number of phosphosites and n is
#' the number of samples
#' @param grps A vector of length n, with group or time point information of the
#'  samples
#'
#' @return A vector of multiple testing adjusted p-values
#'
#' 
matANOVA <- function(mat, grps) {
  ps <- apply(mat, 1, function(x) {
    summary(stats::aov(as.numeric(x) ~ grps))[[1]][["Pr(>F)"]][1]
  })
  
  # adjust for multiple testing
  ps.adj <- stats::p.adjust(ps, method = "fdr")
  return(ps.adj)
}

#' @title Obtain average expression from replicates
#'
#' @usage meanAbundance(mat, grps)
#'
#' @param mat a matrix with rows correspond to phosphosites and columns
#' correspond to samples.
#' @param grps a string specifying the grouping (replciates).
#'
#' @return a matrix with mean expression from replicates
#'
#'
#' 
meanAbundance <- function(mat, grps) {
  # meanMat <- sapply(split(seq_len(ncol(mat)), grps),
  # function(i) rowMeans(mat[,i]))[,unique(grps)]
  meanMat = mapply(function(i, mat) {
    rowMeans(mat[, i], na.rm = TRUE)
  }, split(seq_len(ncol(mat)), grps), MoreArgs = list(mat = mat))[,
                                                                  unique(grps)]
  return(meanMat)
}

############################### Global kinase annotatiion ##

#' @title Kinase substrate scoring
#'
#' @description This function generates substrate scores for kinases that pass
#' filtering based on both motifs and dynamic profiles
#'
#' @param substrate.list A list of kinases with each element containing an array
#'  of substrates.
#' @param mat A matrix with rows correspond to phosphosites and columns
#' correspond to samples.
#' @param seqs An array containing aa sequences surrounding each of all
#' phosphosites.
#' Each sequence has length of 15 (-7, p, +7).
#' @param numMotif Minimum number of sequences used for compiling motif for
#' each kinase.
#' Default is 5.
#' @param numSub Minimum number of phosphosites used for compiling
#' phosphorylation
#' profile for each kinase. Default is 1.
#' @param species Motif list species to be used. Currently there are 
#' \code{mouse} (default), \code{human} and \code{rat}.
#' @param verbose Default to \code{TRUE} to show messages during the progress.
#' All messages will be suppressed if set to \code{FALSE}
#'
#' @importFrom preprocessCore normalize.quantiles
#'
#' @return A list of 4 elements.
#' \code{motifScoreMatrix}, \code{profileScoreMatrix},
#' \code{combinedScoreMatrix}, \code{ksActivityMatrix} (kinase activity matrix)
#' and their \code{weights}.
#'
#' 
kinaseSubstrateScore_local <- function(substrate.list, mat, seqs, numMotif = 5, 
                                       numSub = 1, species = "mouse", verbose = TRUE) {
  
  ks.profile.list <- kinaseSubstrateProfile(substrate.list, mat)
  # motif.mouse.list = PhosR::motif.mouse.list
  data("KinaseMotifs", envir = environment())
  if (!(species %in% c("mouse", "human", "rat"))) {
    stop("Parameter 'species' must be one of 'mouse', 'human' or 'rat'")
  }
  
  if (any(is.na(mat))) {
    stop("Phosphosite quantification matrix contains NAs. Please remove NAs in mat.")
  }
  
  
  if (species == "mouse") {
    motif.list <- motif.mouse.list
  } else if (species == "human") {
    motif.list <- motif.human.list
  } else if (species == "rat") {
    motif.list <- motif.rat.list
  }
  
  if (verbose) {
    message(paste("Number of kinases passed motif size filtering:", 
                  sum(motif.list$NumInputSeq >= numMotif)))
    message(paste("Number of kinases passed profile size filtering:",
                  sum(ks.profile.list$NumSub >= numSub)))
  }
  
  motif.list.filtered <- motif.list[
    which(motif.list$NumInputSeq >= numMotif) ]
  ks.profile.list.filtered <- ks.profile.list[
    which(ks.profile.list$NumSub >= numSub) ]
  # scoring all phosphosites against all motifs
  motifScoreMatrix =
    scorePhosphositesMotifs(mat, motif.list.filtered, seqs, verbose)
  if (verbose) {
    message("done.")
    # scoring all phosphosites against all profiles
    message("Scoring phosphosites against kinase-substrate profiles:")
  }
  profileScoreMatrix = scorePhosphositeProfile(mat, ks.profile.list.filtered)
  if (verbose) {
    message("done.")
    ### prioritisation by integrating the two parts
    message("Generating combined scores for phosphosites
by motifs and phospho profiles:")
  }
  o <- intersect(colnames(motifScoreMatrix), colnames(profileScoreMatrix))
  combinedScoreMatrix <- matrix(NA, nrow = nrow(motifScoreMatrix),
                                ncol = length(o))
  colnames(combinedScoreMatrix) <- o
  rownames(combinedScoreMatrix) <- rownames(motifScoreMatrix)
  # normalising weights for the two parts
  w1 <- log(rank(motif.list$NumInputSeq[o]) + 1)
  w2 <- log(rank(ks.profile.list$NumSub[o]) + 1)
  w3 <- w1 + w2
  for (i in seq_len(length(o))) {
    # weight the two parts by the number of
    # motifs and quantified known substrates
    combinedScoreMatrix[, i] <- (w1[i]/(w1[i] +
                                          w2[i]) * motifScoreMatrix[, o[i]]) +
      (w2[i]/(w1[i] + w2[i]) * profileScoreMatrix[, o[i]])
  }
  if (verbose)
    message("done.")
  # visualise
  message("Preparing table.")
  ksActivityMatrix <- do.call(rbind, ks.profile.list.filtered)[o, ]
  phosScoringMatrices <- list(motifScoreMatrix = motifScoreMatrix,
                              profileScoreMatrix = profileScoreMatrix,
                              combinedScoreMatrix = combinedScoreMatrix,
                              ksActivityMatrix = ksActivityMatrix, weights = w3)
  message("Done.")
  # kinaseSubstrateHeatmap_local(phosScoringMatrices)
  return(phosScoringMatrices)
}

scorePhosphositesMotifs = function(mat, motif.mouse.list.filtered, seqs, 
                                   verbose = TRUE) {
  motifScoreMatrix <- matrix(NA, nrow = nrow(mat),
                             ncol = length(motif.mouse.list.filtered))
  rownames(motifScoreMatrix) <- rownames(mat)
  colnames(motifScoreMatrix) <- names(motif.mouse.list.filtered)
  # extracting flanking sequences
  seqWin = mapply(function(x) {
    mid <- (nchar(x) + 1)/2
    substr(x, start = (mid - 7), stop = (mid + 7))
  }, seqs)
  
  if (verbose)
    message("Scoring phosphosites against kinase motifs:")
  for (i in seq_len(length(motif.mouse.list.filtered))) {
    motifScoreMatrix[, i] <- frequencyScoring(seqWin,
                                              motif.mouse.list.filtered[[i]])
    if (verbose)
      message(paste(i, ".", sep = ""))
  }
  motifScoreMatrix <- minmax(motifScoreMatrix)
  motifScoreMatrix
}

scorePhosphositeProfile = function(mat, ks.profile.list.filtered) {
  profileScoreMatrix <- (t(apply(mat, 1,
                                 cor, t(do.call(rbind, ks.profile.list.filtered)))) + 1)/2
  rownames(profileScoreMatrix) <- rownames(mat)
  colnames(profileScoreMatrix) <- names(ks.profile.list.filtered)
  profileScoreMatrix
}

#' @title Kinase substrate profiling
#'
#' @description This function generates substrate profiles for kinases that have
#' one or more substrates quantified in the phosphoproteome data.
#'
#'
#' @param substrate.list a list of kinases with each element containing an array
#' of substrates.
#' @param mat a matrix with rows correspond to phosphosites and columns
#' correspond to samples.
#' @return Kinase profile list.
#'
#' 
kinaseSubstrateProfile <- function(substrate.list, mat) {
  
  if (any(is.na(mat))) {
    stop("Phosphosite quantification matrix contains NAs. Please remove NAs in mat.")
  }
  
  # generate kinase substrate profile list
  ks.profile.list <- lapply(substrate.list,
                            function(x) {
                              ns <- intersect(x, rownames(mat))
                              m <- c()
                              if (length(ns) == 1) {
                                m <- mat[ns, ]
                              } else if (length(ns) > 1) {
                                m <- apply(mat[ns, ], 2,
                                           median)
                              } else {
                                m <- NA
                              }
                              return(m)
                            })
  
  # ks.profile.list$NumSub <-
  # sapply(substrate.list, function(x){
  # sum(x %in% rownames(mat)) })
  ks.profile.list$NumSub <- mapply(function(x,
                                            mat) {
    sum(x %in% rownames(mat))
  }, substrate.list, MoreArgs = list(mat = mat))
  
  return(ks.profile.list)
}


kinaseActivityHeatmap <- function(ksProfileMatrix) {
  # KinaseFamily = PhosR::KinaseFamily
  utils::data("KinaseFamily", envir = environment())
  o <- intersect(rownames(ksProfileMatrix),
                 rownames(KinaseFamily))
  annotation_row = data.frame(group = KinaseFamily[o,
                                                   "kinase_group"], family = KinaseFamily[o,
                                                                                          "kinase_family"])
  rownames(annotation_row) <- o
  
  pheatmap(ksProfileMatrix, annotation_row = annotation_row,
           fontsize = 7, main = paste("Kinase substrate profiles"))
  
}

#' @title Kinase-substrate annotation prioritisation heatmap
#'
#' @param phosScoringMatrices a matrix returned from kinaseSubstrateScore_local
#' @param top the number of top ranked phosphosites for each kinase to be
#' included in the heatmap. Default is 1.
#' @param printPlot indicate whether the plot should be saved as a PDF
#' in the specified directory. Default is NULL, otherwise specify TRUE.
#' @param filePath path name to save the plot as a PDF file. 
#' Default saves in the working directory.
#' @param width width of PDF.
#' @param height height of PDF.
#'
#' @return a pheatmap object.
#'
#' @import pheatmap
#' @importFrom grDevices dev.off pdf
#'
#' 
kinaseSubstrateHeatmap_local <- function(phosScoringMatrices, top = 3, printPlot=NULL, 
                                         filePath="./kinaseSubstrateHeatmap.pdf", width=20, height=20) {
  # KinaseFamily = PhosR::KinaseFamily
  utils::data("KinaseFamily", envir = environment())
  ####### heatmap 1
  sites <- c()
  for (i in seq_len(ncol(phosScoringMatrices$combinedScoreMatrix))) {
    sites <- union(sites, names(
      sort(phosScoringMatrices$combinedScoreMatrix[,i],
           decreasing = TRUE)[seq_len(top)]))
  }
  o <- intersect(colnames(phosScoringMatrices$combinedScoreMatrix),
                 rownames(KinaseFamily))
  annotation_col = data.frame(group = KinaseFamily[o,
                                                   "kinase_group"], family = KinaseFamily[o,
                                                                                          "kinase_family"])
  rownames(annotation_col) <- o
  if (is.null(printPlot)==TRUE) {
    pheatmap(phosScoringMatrices$combinedScoreMatrix[sites,
    ], annotation_col = annotation_col,
    cluster_rows = TRUE, cluster_cols = TRUE,
    fontsize = 7, main = paste("Top",
                               top, "phosphosite(s) for each kinase"))
    
  } else {
    
    pdf(file=filePath, width=width, height=height)
    
    pheatmap(phosScoringMatrices$combinedScoreMatrix[sites,
    ], annotation_col = annotation_col,
    cluster_rows = TRUE, cluster_cols = TRUE,
    fontsize = 7, main = paste("Top",
                               top, "phosphosite(s) for each kinase"))
    
    dev.off()
  }
}

#' @title Phosphosite annotation
#'
#' @description This function plots the combined scores of each of all kinases
#' for a given phosphosites
#'
#'
#' @param site site the ID of a phosphosite
#' @param phosScoringMatrices output from function kinaseSubstrateScore_local()
#' @param predMatrix a prediction matrix from kinaseSubstratePred()
#'
#' @return A graphical plot
#'
#' 
siteAnnotate <- function(site, phosScoringMatrices,
                         predMatrix) {
  od <- order(predMatrix[site, ], decreasing = TRUE)
  kinases <- colnames(predMatrix)[od]
  # cols <- rainbow(length(kinases))
  
  par(mfrow = c(4, 1))
  barplot(predMatrix[site, kinases], las = 2,
          ylab = "Prediction score", col = "red3",
          main = paste("Site =", site), ylim = c(0,
                                                 1))
  barplot(phosScoringMatrices$combinedScoreMatrix[site,
                                                  kinases], las = 2, ylab = "Combined score",
          col = "orange2", ylim = c(0, 1))
  barplot(phosScoringMatrices$motifScoreMatrix[site,
                                               kinases], las = 2, ylab = "Motif score",
          col = "green4", ylim = c(0, 1))
  barplot(phosScoringMatrices$profileScoreMatrix[site,
                                                 kinases], las = 2, ylab = "Profile score",
          col = "lightblue3", ylim = c(0, 1))
}


#' kinaseSubstratePred
#'
#' A machine learning approach for predicting specific kinase for a given
#' substrate. This prediction framework utilise adaptive sampling.
#'
#' @usage
#' kinaseSubstratePred(
#'     phosScoringMatrices,
#'     ensembleSize = 10,
#'     top = 50,
#'     cs = 0.8,
#'     inclusion = 20,
#'     iter = 5,
#'     verbose = TRUE
#' )
#'
#' @param phosScoringMatrices An output of kinaseSubstrateScore_local
#' @param ensembleSize An ensemble size.
#' @param top a number to select top kinase substrates.
#' @param cs Score threshold.
#' @param inclusion A minimal number of substrates required for a kinase to be
#' selected.
#' @param iter A number of iterations for adaSampling.
#' @param verbose Default to \code{TRUE} to show messages during the progress.
#' All messages will be suppressed if set to \code{FALSE}
#'
#' @return Kinase prediction matrix
#'
#' 
#'
kinaseSubstratePred <- function(phosScoringMatrices,
                                ensembleSize = 10, top = 50, cs = 0.8,
                                inclusion = 20, iter = 5, verbose = TRUE) {
  # create the list of kinase-substrates for prediction
  substrate.list = substrateList(phosScoringMatrices, top, cs, inclusion)
  
  # building the positive training set
  if (verbose)
    message("Predicting kinases for phosphosites:")
  featureMat <- phosScoringMatrices$combinedScoreMatrix
  predMatrix <- matrix(0, nrow = nrow(featureMat),
                       ncol = length(substrate.list))
  colnames(predMatrix) <- names(substrate.list)
  rownames(predMatrix) <- rownames(featureMat)
  
  tmp.list = lapply(seq(length(substrate.list)), function(i) {
    positive.train <- featureMat[substrate.list[[i]],]
    positive.cls <- rep(1, length(substrate.list[[i]]))
    negative.pool <- featureMat[!(rownames(featureMat) %in%
                                    substrate.list[[i]]), ]
    if (verbose)
      message(paste(i, ".", sep = ""))
    tmp_col = predMatrix[,i]
    for (e in seq_len(ensembleSize)) {
      negativeSize <- length(substrate.list[[i]])
      idx <- sample(seq_len(nrow(negative.pool)),
                    size = negativeSize, replace = FALSE)
      negative.samples <- rownames(negative.pool)[idx]
      negative.train <- featureMat[negative.samples,]
      negative.cls <- rep(2, length(negative.samples))
      train.mat <- rbind(positive.train, negative.train)
      cls <- as.factor(c(positive.cls, negative.cls))
      names(cls) <- rownames(train.mat)
      pred <- multiAdaSampling(train.mat,
                               test.mat = featureMat, label = cls,
                               kernelType = "radial", iter = iter)
      tmp_col <- tmp_col[names(pred[, 1])] + pred[, 1]
    }
    tmp_col
  })
  predMatrix = matrix(unlist(tmp.list), ncol = ncol(predMatrix))
  colnames(predMatrix) = names(substrate.list)
  rownames(predMatrix) = rownames(featureMat)
  
  predMatrix <- predMatrix/ensembleSize
  if (verbose)
    message("done")
  return(predMatrix)
}

substrateList = function(phosScoringMatrices, top, cs, inclusion) {
  kinaseSel <- c()
  substrate.list <- list()
  count <- 0
  for (i in seq_len(ncol(phosScoringMatrices$combinedScoreMatrix))) {
    sel <- names(which(sort(phosScoringMatrices$combinedScoreMatrix[, i],
                            decreasing = TRUE)[seq_len(top)] > cs))
    if (length(sel) >= inclusion) {
      count <- count + 1
      substrate.list[[count]] <- sel
      kinaseSel <- c(kinaseSel,
                     colnames(phosScoringMatrices$combinedScoreMatrix)[i])
    }
  }
  if (length(substrate.list)) { # If substrate.list is not empty
    names(substrate.list) <- kinaseSel
  }
  substrate.list
}


#' @import e1071
multiAdaSampling <- function(train.mat, test.mat,
                             label, kernelType, iter = 5) {
  
  X <- train.mat
  Y <- label
  
  model <- c()
  prob.mat <- c()
  
  for (i in seq_len(iter)) {
    tmp <- X
    rownames(tmp) <- NULL
    model <- e1071::svm(tmp, factor(Y),
                        kernel = kernelType, probability = TRUE)
    prob.mat <- attr(predict(model, train.mat,
                             decision.values = FALSE, probability = TRUE),
                     "probabilities")
    
    X <- c()
    Y <- c()
    tmp = lapply(seq(ncol(prob.mat)), function(j) {
      voteClass <- prob.mat[label == colnames(prob.mat)[j], ]
      idx <- c()
      idx <- sample(seq_len(nrow(voteClass)),
                    size = nrow(voteClass), replace = TRUE,
                    prob = voteClass[, j])
      X <- rbind(X, train.mat[rownames(voteClass)[idx],])
      Y <- c(Y, label[rownames(voteClass)[idx]])
      cbind(X,Y)
    })
    tmp = do.call(rbind, tmp)
    Y = tmp[,tail(seq(ncol(tmp)), 1)]
    X = tmp[,-tail(seq(ncol(tmp)),1)]
    
  }
  
  pred <- attr(predict(model, newdata = test.mat,
                       probability = TRUE), "prob")
  return(pred)
}
