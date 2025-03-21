scgsva_custom <- function (obj, annot = NULL, assay = NULL, slot = "counts", batch = 1000, maxRank_UCell = 1500,
          method = "ssgsea", kcdf = "Poisson", abs.ranking = FALSE, 
          min.sz = 1, max.sz = Inf, mx.diff = TRUE, ssgsea.norm = TRUE, 
          useTerm = TRUE, BPPARAM = SnowParam(), cores = 4, verbose = TRUE, 
          sc.keep = TRUE, ...) {
  tau = switch(method, gsva = 1, ssgsea = 0.25, NA)
  if (is.null(annot)) {
    stop("Please provide annotation object or data.frame")
  }
  else {
    if (isTRUE(useTerm)) {
      annotation <- split(annot[, 1], annot[, 3])
    }
    else {
      annotation <- split(annot[, 1], annot[, 2])
    }
  }
  if (inherits(x = obj, what = "Seurat")) {
    if (is.null(assay)) 
      assay <- "RNA"
    input <- GetAssayData(obj, assay = assay, layer = slot)
    input <- input[tabulate(summary(input)$i) != 0, , drop = FALSE]
    input <- as.matrix(input)
  }
  else if (inherits(x = obj, what = "SingleCellExperiment")) {
    input <- counts(obj)
    if (!"logcounts" %in% names(assays(obj))) {
      libsizes <- colSums(assay(obj, "counts"))
      size.factors <- libsizes/mean(libsizes)
      logcounts(obj) <- as.matrix(log(t(t(input)/size.factors) + 
                                        1))
    }
    else {
      logcounts(obj) <- as.matrix(logcounts(obj))
    }
    input <- input[tabulate(summary(input)$i) != 0, , drop = FALSE]
    input <- as.matrix(input)
    obj <- as.Seurat(obj)
  }
  else {
    input <- obj
  }
  input <- input[rowSums(input > 0) != 0, ]
  if (method == "UCell") {
    out <- suppressWarnings(UCell::ScoreSignatures_UCell(input, maxRank = maxRank_UCell,
                                                         features = annotation, chunk.size = batch, ncores = cores, 
                                                         BPPARAM = SerialParam(progressbar = verbose), ...))
    colnames(out) <- gsub(" ", "\\.", sub("_UCell", "", colnames(out)))
    out <- as.data.frame(out)
  }
  else {
    if (ncol(obj) > batch) {
      split.data <- split.data.matrix(matrix = input, chunk.size = batch)
      out <- lapply(split.data, function(x) .sgsva(input = x, 
                                                   annotation = annotation, method = method, kcdf = kcdf, 
                                                   abs.ranking = abs.ranking, min.sz = min.sz, max.sz = max.sz, 
                                                   cores = cores, tau = tau, ssgsea.norm = FALSE, 
                                                   verbose = verbose))
      out <- do.call(rbind, out)
      if (isTRUE(ssgsea.norm)) {
        rng <- range(out)
        if (any(is.na(rng) | !is.finite(rng))) 
          rng <- range(out, na.rm = TRUE)
        out <- out[1:nrow(out), , drop = FALSE]/(rng[2] - 
                                                   rng[1])
      }
    }
    else {
      out <- .sgsva(input = input, annotation = annotation, 
                    method = method, kcdf = kcdf, abs.ranking = abs.ranking, 
                    min.sz = min.sz, max.sz = max.sz, cores = cores, 
                    tau = tau, ssgsea.norm = ssgsea.norm, verbose = verbose)
    }
  }
  annot <- annot[annot[, 1] %in% rownames(input), ]
  if (isTRUE(useTerm)) {
    annot <- annot[order(annot[, 3]), ]
  }
  else {
    annot <- annot[order(annot[, 2]), ]
  }
  if (!isTRUE(sc.keep)) {
    if (is.null(assay)) 
      assay <- "RNA"
    empty_counts <- Matrix::Matrix(0, nrow = 0, ncol = 0, 
                                   sparse = TRUE)
    empty_counts <- as(as(empty_counts, "generalMatrix"), 
                       "CsparseMatrix")
    obj <- SetAssayData(object = obj, assay = "RNA", layer = "counts", 
                        new.data = empty_counts)
  }
  res <- new("GSVA", obj = obj, gsva = out, annot = annot)
  return(res)
}
environment(scgsva_custom) <- environment(scgsva)
