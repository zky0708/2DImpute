#' Imputation based on cell-to-cell relationships
#' 
#' This function imputes dropout-suspected zeros by 
#' borrowing information from similar cells using
#' k-nearest neighbors (kNN) regression. 
#' 
#' @param exprs A log-transformed expression matrix, where
#' rows and columns corresponds to genes and cells, respectively.
#' 
#' @param dropout_ind A matrix where each row contains the 
#' indices of zero entries (first column for gene name, 
#' second column for cell ID) that are suspected as dropouts.
#' 
#' @param k Number of neighbors to be used in the kNN regression.
#' Default is 10.
#' 
#' @param ncores Number of cores to use. Default is 1.
#' 
#' @param verbose Whether to show the progress of imputation.
#' Default is TRUE.
#' 
#' 
#' @return An imputed version of gene expression matrix. 
#' 
#' 
#' @export


imputeByCells <- function(exprs, dropout_ind, k = 10, ncores = 1, verbose = TRUE){
  
  if(verbose)
    cat(nrow(dropout_ind), 'entries remained to impute. \n')
  
  imputed = exprs
  
  # only impute zeros in genes that have not been processed in the previous 'imputeByAttractors'
  genes_all <- rownames(exprs)
  cells_all <- colnames(exprs)
  
  cells_to_impute <- names(sort(table(dropout_ind[,2]), decreasing = F))
  cells_nonzero <- apply(exprs, 1, function(x) names(which(x > 0)))
  
  count = length(cells_to_impute)
  for (c in cells_to_impute) {
    
    genes_to_impute <- dropout_ind[which(dropout_ind[,2] == c),1]
    
    # use genes excluding those to be imputed to calculate the pairwise cell-to-cell correlation
    tmp = setdiff(genes_all, genes_to_impute)

    P <- unlist(parallel::mclapply(cells_all, 
                                   function(x) cor(imputed[tmp, x], imputed[tmp, c]), mc.cores = ncores))
    names(P) <- cells_all
    
    nearest_neighbors <- names(sort(P, decreasing = T)[2:(k+1)])
    tmp = imputed[genes_to_impute, nearest_neighbors]; 
    tmp[tmp == 0] = NA; 
    if(length(genes_to_impute) > 1)
      imputed[genes_to_impute, c] = rowMeans(tmp, na.rm = T)
    else
      imputed[genes_to_impute, c] = mean(tmp, na.rm = T)
    imputed[is.na(imputed)] <- 0
    
    count = count - 1
    if(verbose && count %% 100 == 0)
      cat(count, 'cells remained to impute. \n')
    
    
  }
  
  
  return(imputed)
  
}