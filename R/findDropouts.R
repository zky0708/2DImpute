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
#' @param verbose Whether to show the progress of dropout identification.
#' Default is TRUE.
#'
#' @return A list with the following components:
#' \item{dropout_ind}{A matrix in which each row contains the indices of a dropout}
#' \item{J}{Calculated pairwise Jaccard Index between cells}
#'
#' @examples 
#' data(ge_10x_sample)
#' results1 <- findDropouts(ge_10x_sample)
#' # pairwise Jaccard Index between cells
#' J <- results1$J
#' # indices of identified dropouts
#' dropout_ind <- results1$dropout_ind
#' 
#'
#' @export


# source('~/Downloads/imputation_project_code/jaccard_index.R')

findDropouts <- function(exprs, t = 0.2, genes = NULL, verbose = TRUE){

  if(verbose)
    cat('Computing pairwise cell-cell Jaccard Index... \n')

  J = jaccard_index(Matrix = exprs)

  zero_ind <- which(exprs == 0, arr.ind = T)
  zero_ind[,1] <- rownames(exprs)[zero_ind[,1]]
  zero_ind[,2] <- colnames(exprs)[as.integer(zero_ind[,2])]
  cells_all <- unique(zero_ind[,2])

  if(verbose)
    cat('Identifying dropout-suspected zeros... \n')

  prob <- rep(0, length = nrow(zero_ind))
  count = 0
  for(c in cells_all){

    cell_id <- which(zero_ind[,2] == c)

    # identify similar cells according to Jaccard index
    similar_cells <- setdiff(colnames(exprs)[J[,c] >= 0.5], c)

    if(length(similar_cells) < 10)
      # next
      similar_cells <- colnames(exprs)[order(J[, c], decreasing = T)[2:11]]
    
    # estimate dropout likelihood for each zero-expressed gene in cell c
    genes_zero <- zero_ind[cell_id, 1]
    L = length(similar_cells)
    prob[cell_id] <- sapply(genes_zero, function(x) sum(exprs[x, similar_cells] > 0)/L)

    count = count + 1
    if(verbose && count %% 100 == 0)
      cat(count, 'cells completed. \n')

  }

  dropout <- prob > t
  dropout_ind <- zero_ind[which(dropout == 1),]
  
  if(!is.null(genes))
    dropout_ind <- dropout_ind[which(dropout_ind[,1] %in% genes),]

  results = list(dropout_ind = dropout_ind, J = J)
  return(results)

}

