#' Calculate pairwise Jaccard Index between cells
#'
#' This function calculates the Jaccard Index between each pair of cells
#' based on binarized expression matrix.
#'
#' @param Matrix An input matrix, where rows and columns represent
#' genes and cells, respectively.
#'
#' @param threshold The threshold for generating binary matrix.
#' Values greater than it are set to 1, otherwise 0. Default is 0.
#'
#' @param ncores Number of cores to use. Default is 1.
#'
#' @return A symmetric matrix containing pairwise Jaccard Index values
#' between cells. 
#' 
#' 
#' @export



jaccard_index <- function(Matrix, threshold = 0, ncores = 1){

  Matrix[Matrix > threshold] = 1
  Matrix[Matrix <= threshold] = 0

  n = ncol(Matrix)
  J = matrix(NA, ncol = n, nrow = n)
  rownames(J) <- colnames(J) <- colnames(Matrix)

  if(ncores > 1){
    J_list <- parallel::mclapply(1:(n-2),
                               function(i) apply(Matrix[, (i+1):n], 2, function(y) jaccard::jaccard(Matrix[,i], y)),
                               mc.cores = ncores)

    for (i in 1:(n-2)) {
      J[i, (i+1):n] <- J[(i+1):n, i] <- J_list[[i]]
    }

  }else{
    
    for(i in 1:(n-2)){
      
      x = Matrix[,i]
      tmp = Matrix[, (i+1):n]
      J[i, (i+1):n] <- J[(i+1):n, i] <- apply(tmp, 2, function(y) jaccard::jaccard(x,y))
      
    }
    
    
  }

  J[(n-1), n] <- J[n, n-1] <- jaccard::jaccard(Matrix[, (n-1)], Matrix[, n])
  
  diag(J) <- 1

  return(J)

}
