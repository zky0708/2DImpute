#' Run 2DImpute on an expression matrix
#'
#' 2DImpute is an imputation algorithm designed for
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
#' @param k Number of neighbors to be used in the k-nearest 
#' neighbor regression in the step of \code{\link{imputeByCells}}.
#' Default is 10.
#' 
#' @param ncores Number of cores to be used. Default is 1.
#' 
#' @param log_input Whether the input data matrix is log-transformed.
#' Default is TRUE. 
#' 
#' @param return_J Whether to return the calculated pairwise
#' Jaccard matrix between cells in the step of \code{\link{findDropouts}}. 
#' Default is FALSE. 
#' 
#' @param return_attractors Whether to return the co-expressed
#' gene attractors found by the function of \code{\link{attractorFinding}}.
#' Default is FALSE.
#' 
#' @param verbose Whether to show the progress of imputation.
#' Default is TRUE.
#'
#'
#' @return A list the following components:
#' \item{imputed}{The imputed version of expression data matrix}
#' \item{attractors}{(optional) Identified co-expressed gene attractor signatures}
#' 
#' 
#' @examples 
#' data(ge_10x_sample)
#' 
#' # Fast demonstration
#' res <- run_2DImpute(exprs = ge_10x_sample, genes = c('XIST', 'CD3D'), ncores = 2)
#' imputed_exprs <- res$imputed
#' 
#' # Return identified co-expressed attractor signatures
#' res <- run_2DImpute(exprs = ge_10x_sample, genes = c('XIST', 'CD3D'), ncores = 2, return_attractors = TRUE)
#' imputed_exprs <- res$imputed
#' attractors <- res$attractors
#' 
#' # Full imputation
#' res <- run_2DImpute(exprs = ge_10x_sample, ncores = 2)
#' imputed_exprs <- res$imputed
#'
#' @export


run_2DImpute <- function(exprs, t = 0.2, genes = NULL, k = 10, ncores = 1, log_input = TRUE,
                         return_J = FALSE, return_attractors = FALSE, verbose = TRUE){
  
  if(!log_input)
    exprs <- log2(exprs + 1)
  
  
  results1 <- findDropouts(exprs = exprs, t = t, genes = genes, ncores = ncores, 
                           return_J = return_J, verbose = verbose)
  
  dropout1 <- results1$dropout_ind
  
  attractors <- attractorFinding(exprs = exprs, genes_involved = unique(dropout1[,1]), verbose = verbose)
  
  results2 <- imputeByAttractors(exprs = exprs, dropout_ind = dropout1, 
                                 attractor_list = attractors, verbose = verbose)
  
  imputed_exprs <- imputeByCells(exprs = results2$imputed, dropout_ind = results2$dropout_ind, 
                           k = k, ncores = ncores, verbose = verbose)
  
  if(!log_input)
    imputed_exprs <- round(2^imputed_exprs - 1, digits = 2)
  
  if(return_attractors)
    return(list(imputed = imputed_exprs, attractors = attractors))
  else
    return(list(imputed = imputed_exprs))
  
}
