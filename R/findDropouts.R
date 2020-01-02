#' Identify dropout-suspected zeros from all zero entries
#'
#' This function distinguishes zeros that are likely to be 
#' dropouts from those biologically true zeros based on 
#' pairwise Jaccard Index between cells.
#'
#' @param exprs An input log-transformed data matrix. The rows and columns
#' correspond to genes and cells, respectively.
#'
#' @param t Threshold (between 0 and 1) for determining similar cells. 
#' Default is 0.2.
#' 
#' @param genes A vector containing gene symbols that will 
#' get imputed. Default is NULL, in which case all available 
#' genes in the matrix exprs will be imputed.
#' 
#' @param ncores Number of cores to use. Default is 1.
#' 
#' @param return_J Whether to return the calculated pairwise Jaccard matrix between cells. 
#' Default is FALSE.
#' 
#' @param verbose Whether to show the progress of dropout identification.
#' Default is TRUE.
#'
#' @return A list with the following components:
#' \item{dropout_ind}{A matrix in which each row contains the indices of a dropout}
#' \item{J}{(optional) Calculated pairwise Jaccard Index between cells}
#'
#' @examples 
#' data(ge_10x_sample)
#' results1 <- findDropouts(ge_10x_sample, ncores = 2, return_J = TRUE)
#' # pairwise Jaccard Index between cells
#' J <- results1$J
#' # indices of identified dropouts
#' dropout_ind <- results1$dropout_ind
#' 
#'
#' @export


# source('~/Downloads/imputation_project_code/jaccard_index.R')

findDropouts <- function(exprs, t = 0.2, genes = NULL, ncores = 1, return_J = FALSE, verbose = TRUE){

  if(verbose)
    cat('Computing pairwise cell-cell Jaccard Index... \n')

  J = jaccard_index(Matrix = exprs, ncores = ncores)

  zero_ind <- which(exprs == 0, arr.ind = T)
  zero_ind[,1] <- rownames(exprs)[zero_ind[,1]]
  zero_ind[,2] <- colnames(exprs)[as.integer(zero_ind[,2])]
  
  if(!is.null(genes))
    zero_ind <- zero_ind[which(zero_ind[,1] %in% genes),]
  
  cells_all <- unique(zero_ind[,2])

  if(verbose)
    cat('Identifying dropout-suspected zeros... \n')

  similar_cells <- parallel::mclapply(cells_all, 
                                      function(c) order(J[,c], decreasing = TRUE)[2:max(11, sum(J[,c] >= 0.5))], 
                                      mc.cores = ncores)
  names(similar_cells) <- cells_all

  prob_list <- parallel::mclapply(cells_all, 
                                  function(c) sapply(zero_ind[which(zero_ind[,2] == c),1], 
                                                     function(x) mean(exprs[x, similar_cells[[c]]] > 0)), 
                                  mc.cores = ncores)

  prob <- unlist(prob_list, use.names = F)
  dropout <- prob > t
  dropout_ind <- zero_ind[which(dropout == 1),]
  
  if(verbose)
    cat(paste(nrow(dropout_ind), 'out of', nrow(zero_ind), 'zeros (',
              round(nrow(dropout_ind)/nrow(zero_ind)*100, digits = 2), 
              '% ) are suspected as dropouts. \n'))
  
  

  if(return_J)
    results = list(dropout_ind = dropout_ind, J = J)
  else
    results = list(dropout_ind = dropout_ind)
  
  return(results)

}

