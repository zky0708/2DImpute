#' Find all co-expressed signatures in the expression matrix
#' 
#' This function uses heuristic exhaustive search to find
#' existing co-expressed signatures in the expression matrix,
#' with attractor finding algorithm implemented in the
#' R package \pkg{cafr}.
#' 
#' @param exprs A log-transformed expression data matrix of which 
#' rows and columns represent genes and cells, respectively.
#' 
#' @param dropout_ind A matrix in which each row contains the indices
#' of zero entries that are suspected as dropouts.
#' 
#' @param a Exponent of mutual information (MI) used to create weighted
#' vector for metagenes, one of parameters in function 
#' \code{findAttractors} in package \pkg{cafr}. Default is 5.
#' 
#' @param epsilon Convergence threshold, one of parameters in function
#' \code{findAttractors} in package \pkg{cafr}. Default is 1E-4.
#' 
#' @param mi_threshold MI threshold (between 0 to 1) to determine the 
#' gene membership of signatures. Genes with greater than or equal to 
#' the threshold are considered to be associated with the signature. 
#' Default is 0.4.
#' 
#' @param verbose Whether to show the finding progress. Default is TRUE. 
#' 
#' 
#' @return An attractor list in which each element contains the 
#' vector of MIs between the genes and the converged metagene 
#' for each attractor signature. 
#' 
#' @examples 
#' data(ge_10x_sample)
#' data(results1)
#' attractors <- attractorFinding(ge = ge_10x_sample, dropout_ind = results1$dropout_ind)
#' 
#' 
#' @export


attractorFinding <- function(exprs, dropout_ind, a = 5, epsilon = 1E-4, mi_threshold = 0.4, verbose = TRUE){
  
  genes_involved = unique(dropout_ind[,1])
  
  if(length(genes_involved) == nrow(exprs)){
    
    # Order genes according to standard deviation (SD), starting with the one with highest SD
    genes_zero <- unique(dropout_ind[,1])
    sd_genes <- apply(exprs[genes_zero,], 1, sd)
    genes_zero <- genes_zero[order(sd_genes, decreasing = T)]
    
    seeds <- genes_zero
    
  }else{
    
    seeds <- genes_involved
  
  }
  
  attractors <- NULL
  topgenes <- NULL
  # seeds <- genes_zero
  
  if(verbose)
    cat('Identifying co-expressed attractor signatures...\n')
  
  while(length(seeds) > 0){
    
    i = seeds[1]
    
    mi <- cafr::getAllMIWz(exprs, exprs[i,], negateMI = T, sorting = T)
    
    # Skip the gene if the highest MI with other genes < 0.2
    if(mi[2] < 0.2){
      
      seeds <- seeds[-1]
      next
      
    }else{
      
      # metagene <- cafr::findAttractor(ge, ge[i,], a = a, epsilon = epsilon, verbose = F)
      metagene <- findAttractor_new(exprs, exprs[i,], i, a=a, threshold = mi_threshold, epsilon=epsilon)
      if(is.null(metagene)){
        seeds <- seeds[-1]
        next
      }
      
      topgene <- names(metagene[1])
      
      # Check if the attractor has already been found before
      if(topgene %in% topgenes){
        
        seeds <- seeds[-1]
        
      }else{
        
        # Remove genes from the seed pool with scores greater than that of the seed 
        # for computational efficiency.
        seeds <- setdiff(seeds[-1], names(metagene)[which(metagene >= metagene[i])])

        # Only retain attractors in which the 5th ranked gene has score >= mi_threshold (default: 0.4)
        # Only retain attractors in which the score difference between the 1st and 2nd genes < 0.2
        if(metagene[5] >= mi_threshold && (metagene[1]-metagene[2]) < 0.2){
          
          topgenes <- c(topgenes, topgene)
          attractors[[topgene]] <- metagene[which(metagene > 0)]
          if(verbose)
            cat(length(attractors), 'attractors found and', length(seeds), 'genes left to scan. \n')
          
        }
        
      }
      
      
    }
    
    
  }
  
  return(attractors)
  
}

findAttractor_new <- function(data, vec, seed, threshold = 0.5, a=5, maxIter = 100, 
                              epsilon=1E-14, bin = 6, so = 3,rankBased = FALSE,  
                              negateMI = TRUE, verbose = FALSE){
  
  # modified from findAttractor() function in cafr R package (https://github.com/weiyi-bitw/cafr)
  
  dataIn <- data
  if(rankBased){
    vec <- rank(vec)
    for(i in 1:nrow(data)){
      dataIn[i,] <- rank(data[i,])
    } 
  }
  mi <- cafr::getAllMIWz(dataIn, vec, bin=bin, so=so, negateMI = negateMI)
  premi <- mi
  #   w <- abs(mi)^a / sum(abs(mi)^a)
  #   w[mi < 0] <- 0
  mi[mi < 0] <- 0
  w <- abs(mi)^a / sum(abs(mi)^a)
  metagene <- w %*% data
  if(rankBased) metagene <- rank(metagene)
  c <- 0
  while(c < maxIter){
    mi <- cafr::getAllMIWz(dataIn, metagene, bin=bin, so=so, negateMI = negateMI)
    
    if(mi[seed] < threshold)  return(NULL)
    
    delta <- sum((mi - premi)^2)
    if(verbose){
      cat("Iteration ", (c+1), "\tDelta = ", delta, "\n", sep="")
      print(mi[order(mi, decreasing=T)[1:20]])
      flush.console()
    }
    if(delta < epsilon){
      break
    }
    premi <- mi
    mi[mi < 0] <- 0
    w <- abs(mi)^a / sum(abs(mi)^a)
    # w[mi < 0] <- 0
    metagene <- w %*% data
    if(rankBased) metagene <- rank(metagene)
    c <- c + 1
  }
  if(c >= maxIter) return (NULL)
  return (sort(mi, decreasing=T))
}


