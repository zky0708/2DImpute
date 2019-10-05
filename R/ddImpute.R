#' Run ddImpute on an expression matrix
#'
#' ddImpute is an imputation algorithm designed for
#' recovering dropouts in single-cell RNA-seq data
#' by using information from two-dimensional relationships 
#' (cells and genes), as described in Zhu et al. 2019.
#' 
#' @param exprs A log-transformed gene expression data matrix,
#' where rows and columns correspond to genes and cells, 
#' respectively.
#' 
#' @param t Threshold (between 0 and 1) for determining 
#' similar cells in the step of dropout identification. 
#' Default is 0.2.
#' 
#' @param genes A vector containing gene symbols that will 
#' get imputed. Default is NULL, in which case all available 
#' genes in the matrix exprs will be imputed.
#' 
#' @param mi_threshold MI threshold (between 0 to 1) to 
#' determine the gene membership of signatures. Genes with 
#' greater than or equal to the threshold are considered to be 
#' associated with the signature. Default is 0.4.
#' 
#' @param k Number of neighbors to be used in the k-nearest 
#' neighbor regression in the step of \code{imputeByCells}.
#' Default is 10.
#' 
#' @param verbose Whether to show the progress of imputation.
#' Default is TRUE.
#'
#'
#' @return The imputed version of expression data matrix.
#' 
#' 
#' @examples 
#' data(ge_10x_sample)
#' imputed <- ddImpute(exprs = ge_10x_sample)
#' 
#'
#' @export


ddImpute <- function(exprs, t = 0.2, genes = NULL, mi_threshold = 0.4, k = 10, verbose = TRUE){
  
  results1 <- findDropouts(exprs = exprs, t = t, genes = genes, verbose = verbose)
    
  attractors <- attractorFinding(exprs = exprs, dropout_ind = results1$dropout_ind, 
                                 mi_threshold = mi_threshold, verbose = verbose)
  
  results2 <- imputeByAttractors(exprs = exprs, dropout_ind = results1$dropout_ind, 
                                 attractor_list = attractors, mi_threshold = mi_threshold,
                                 verbose = verbose)
  
  imputed <- imputeByCells(exprs = results2$imputed, dropout_ind = results2$dropout_ind, 
                           k = k, verbose = verbose)
  
  return(imputed)
  
}
