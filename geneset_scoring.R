library(scGSVA)
source("~/scgsva_custom.R")
# source("~/Downloads/load_msigdb_datasets.R")
#' Score Gene Sets from Single-Cell Data Stored in an H5AD File
#'
#' @description
#' This function reads a single-cell dataset from an H5AD file and computes
#' gene set enrichment scores via the [scGSVA](https://cran.r-project.org/package=scGSVA)
#' framework. It supports multiple enrichment methods, including "UCell"
#' (from the [UCell](https://cran.r-project.org/package=UCell) package) and
#' "ssgsea"/"gsva" from GSVA-based algorithms. 
#'
#' After loading the H5AD file using [zellkonverter::readH5AD()],
#' the function configures assay data for subsequent scoring and removes duplicated
#' genes. It then builds an MSigDB annotation object based on user-specified criteria,
#' filtering out gene sets that do not match the expression data. For "UCell" scoring,
#' the function further restricts gene sets by their size (relative to \code{maxRank_UCell})
#' and optionally subsets them if desired. Finally, it uses the customized
#' \code{scgsva_custom()} function from `scGSVA` to generate enrichment scores
#' which are returned as a \code{SingleCellExperiment} object (unless otherwise configured).
#'
#' @param h5ad_path `character(1)`  
#'   Path to the H5AD file containing single-cell data.  
#'   Defaults to \code{NULL}.  
#'
#' @param assay_name `character(1)`  
#'   Name of the assay in the SingleCellExperiment object on which gene set
#'   scoring should be computed.  
#'   Default: \code{"log1p"}.  
#'
#' @param scoring_method `character(1)`  
#'   The enrichment scoring method to use. Options typically include
#'   \code{"ssgsea"}, \code{"gsva"}, or \code{"UCell"}.  
#'   Default: \code{"UCell"}.  
#'
#' @param anntype `character(1)`  
#'   The type of annotation to use for building the MSigDB object (e.g., \code{"GO"},
#'   \code{"KEGG"}, etc.).  
#'   No default is set inside the function (must be specified if \code{scoring_method} 
#'   requires it).  
#'
#' @param num_workers `numeric(1)`  
#'   Number of parallel workers to use. By default, it is set to \code{snowWorkers()},
#'   or \code{min(num_workers, snowWorkers())} if you specify it.  
#'
#' @param batch `numeric(1)`  
#'   Number of cells processed in each batch during scoring.  
#'   Default: 1000.  
#'
#' @param kcdf `character(1)`  
#'   Kernel to use during non-parametric estimation of the cumulative distribution
#'   function when using GSVA-based methods (e.g. \code{"Poisson"} or \code{"Gaussian"}).  
#'   Default: \code{"Poisson"}.  
#'
#' @param maxRank_UCell `numeric(1)`  
#'   Maximum number of genes to rank per cell for the UCell method. Above this rank,
#'   genes are considered not expressed.  
#'   Default: 1500.  
#'
#' @param abs.ranking `logical(1)`  
#'   Flag used only when \code{mx.diff=TRUE}; see \code{\link[GSVA]{gsva}} for details.  
#'   Default: \code{FALSE}.  
#'
#' @param min.sz `numeric(1)`  
#'   Minimum size of the resulting gene sets to consider for enrichment.  
#'   Default: 1.  
#'
#' @param max.sz `numeric(1)`  
#'   Maximum size of the resulting gene sets to consider for enrichment.  
#'   Default: \code{Inf}.  
#'
#' @param mx.diff `logical(1)`  
#'   Switch determining how the enrichment statistic (ES) is calculated for the
#'   GSVA-based methods.  
#'   Default: \code{TRUE}.  
#'
#' @param ssgsea.norm `logical(1)`  
#'   Whether or not to normalize the resulting ssgsea scores; applies only
#'   if \code{method="ssgsea"}.  
#'   Default: \code{TRUE}.  
#'
#' @param useTerm `logical(1)`  
#'   If \code{TRUE}, the function uses descriptive Terms from the MSigDB annotation
#'   instead of numeric IDs.  
#'   Default: \code{TRUE}.  
#'
#' @param verbose `logical(1)`  
#'   If \code{TRUE}, prints progress messages and warnings.  
#'   Default: \code{TRUE}.  
#'
#' @param sc.keep `logical(1)`  
#'   If \code{TRUE}, the returned object retains the original single-cell data.  
#'   If \code{FALSE}, only gene set enrichment scores are kept.  
#'   Default: \code{FALSE}.  
#'
#' @details
#' \strong{Data Reading & Cleaning}: The function uses [zellkonverter::readH5AD()] to
#' read the single-cell dataset. It then sets the `counts` and `logcounts` slots using
#' the `"raw"` and `"log1p"` assay data, respectively. Any duplicated gene names are
#' detected, printed, and removed.
#'
#' \strong{Annotation}: The user must specify \code{anntype} (e.g. "GO", "KEGG", etc.).
#' The function builds an MSigDB object for \code{species="human"} and filters it by
#' gene symbols found in \code{rownames(so)}.
#'
#' \strong{UCell Restriction}: If \code{scoring_method="UCell"}, gene sets with sizes
#' larger than \code{maxRank_UCell} are excluded, and for demonstration, the code
#' retains only the first 20 gene sets (remove or adjust this line in production).
#'
#' \strong{scGSVA Call}: The final scoring is done by [scgsva_custom()] (a wrapper
#' from `scGSVA`) with arguments mapped to the user inputs. The function returns
#' a \code{SingleCellExperiment} object containing either just the scores or both
#' the scores and the single-cell data, depending on \code{sc.keep}.
#'
#' For details on the scoring methods, see:
#' - \pkg{scGSVA} [scgsva() documentation](https://cran.r-project.org/package=scGSVA)  
#' - \pkg{UCell} [ScoreSignatures_UCell() documentation](https://cran.r-project.org/package=UCell)  
#'
#' @return
#' A \code{SingleCellExperiment} object with additional gene set enrichment scores stored.
#' If \code{sc.keep=TRUE}, the object contains the full single-cell data along with
#' gene set scores. Otherwise, it is reduced to enrichment-relevant information.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' result_sce <- geneset_score(
#'   h5ad_path     = "path/to/file.h5ad",
#'   assay_name    = "log1p",
#'   scoring_method= "UCell",
#'   anntype       = "GO",
#'   num_workers   = 4,
#'   batch         = 1000,
#'   kcdf          = "Poisson",
#'   maxRank_UCell = 1500,
#'   abs.ranking   = FALSE,
#'   min.sz        = 1,
#'   max.sz        = Inf,
#'   mx.diff       = TRUE,
#'   ssgsea.norm   = TRUE,
#'   useTerm       = TRUE,
#'   verbose       = TRUE,
#'   sc.keep       = TRUE
#' )
#'
#' # The returned object can be further processed or accessed as needed:
#' result_scores <- SingleCellExperiment::assay(result_sce, "scgsva_score")
#' head(result_scores[, 1:5])
#' }
#'
#' @import zellkonverter
#' @import scGSVA
#' @importFrom BiocParallel SnowParam
#' @importFrom methods slot
#' @export


geneset_score <- function(h5ad_path = NULL,
                          assay_name = "log1p",
                          scoring_method = "UCell",
                          anntype = anntype,
                          num_workers = snowWorkers(),
                          batch = 1000,
                          kcdf = "Poisson",
                          maxRank_UCell = 1500, # Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed. Note: this parameter is ignored if precalc.ranks are specified
                          abs.ranking = FALSE,
                          min.sz = 1,
                          max.sz = Inf,
                          mx.diff = TRUE,
                          ssgsea.norm = TRUE,
                          useTerm = TRUE,
                          verbose = TRUE,
                          sc.keep=FALSE
                          
                          ){
  so <- zellkonverter::readH5AD(h5ad_path)
  counts(so) <- assays(so)$raw
  logcounts(so) <- assays(so)$log1p
  
  gene_names <- rownames(assays(so)$raw)
  gene_names <- gsub(pattern = "_", replacement = "-", gene_names)
  duplicated_gene_indices <- which(duplicated(gene_names))
  ## print duplicated gene names
  print(gene_names[duplicated_gene_indices])
  ## remove duplicated gene names
  so <- so[-duplicated_gene_indices,]
  

  num_workers <- min(num_workers, snowWorkers())
  parallel_obj <- SnowParam(workers = num_workers)
  annot_obj <- buildMSIGDB(species="human",keytype="SYMBOL",anntype=anntype)
  annot_obj@annot <- annot_obj@annot[annot_obj@annot$gene_symbol %in% rownames(so),]
  
  ###subset annot_obj to only include genesets of size less than maxRank_UCell
  if(scoring_method == "UCell"){
    annot_df <- ii@annot
    annot_list <- split(annot_df, annot_df$Annot)
    annot_gs_size <- sapply(annot_list, nrow)
    annot_gs_size_df <- data.frame(Annot = names(annot_gs_size), gs_size = annot_gs_size)
    annot_df <- annot_df %>% dplyr::left_join(annot_gs_size_df, by = "Annot")
    annot_df <- annot_df[annot_df$gs_size <= maxRank_UCell,]
    annot_df$gs_size <- NULL
    annot_df <- annot_df[annot_df$Annot %in% unique(annot_df$Annot)[1:20],]
    annot_obj@annot <- annot_df
  }
  
  
  print(head(annot_obj))
  
  enr_obj <- scgsva_custom(
    so,
    annot = annot_obj,
    assay = assay_name,
    batch = batch,
    method = scoring_method,
    kcdf = kcdf,
    abs.ranking = abs.ranking,
    min.sz = min.sz,
    max.sz = max.sz,
    mx.diff = mx.diff,
    ssgsea.norm = ssgsea.norm,
    useTerm = useTerm,
    BPPARAM = parallel_obj,
    cores = num_workers,
    verbose = verbose,
    sc.keep=sc.keep,
  )
  return(enr_obj)
}

## test eg:
# kk <- geneset_score(h5ad_path = "/Users/u0119129/WORK/COST/b_macro_subset.h5ad",
#                     assay_name = "log1p",
#                     scoring_method = "UCell",
#                     anntype = "GO",
#                     num_workers = 4,
#                     batch = 1000,
#                     kcdf = "Poisson",
#                     maxRank_UCell = 1500, # Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed. Note: this parameter is ignored if precalc.ranks are specified
#                     abs.ranking = FALSE,
#                     min.sz = 1,
#                     max.sz = Inf,
#                     mx.diff = TRUE,
#                     ssgsea.norm = TRUE,
#                     useTerm = TRUE,
#                     verbose = TRUE,
#                     sc.keep=TRUE
# )