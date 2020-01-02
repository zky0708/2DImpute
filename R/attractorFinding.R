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
#' @param genes_involved A vector containing genes that will be used
#' as seeds in the attractor finding process. Default is NULL, in which
#' case all genes in the expression matrix will be used.
#' 
#' @param a Exponent of mutual information (MI) used to create weighted
#' vector for metagenes, one of parameters in function 
#' \code{\link[cafr]{findAttractor}} in package \pkg{cafr}. Default is 5.
#' 
#' @param epsilon Convergence threshold, one of parameters in function
#' \code{\link[cafr]{findAttractor}} in package \pkg{cafr}. Default is 1E-4.
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
#' attractors <- attractorFinding(exprs = ge_10x_sample, genes_involved = c('CD3D', 'XIST'))
#' 
#' 
#' @export


attractorFinding <- function(exprs, genes_involved = NULL, a = 5, epsilon = 1E-4, mi_threshold = 0.4, verbose = TRUE){
  
  if(is.null(genes_involved))
    genes_involved = rownames(exprs)
  
  # Order genes according to standard deviation (SD), starting with the one with highest SD
  sd_genes <- apply(exprs[genes_involved,], 1, sd)
  # seeds <- genes_involved[order(sd_genes, decreasing = T)[1:min(2000, length(genes_involved))]]
  seeds <- genes_involved[order(sd_genes, decreasing = T)]
  
  attractors <- NULL
  topgenes <- NULL

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

        # Only retain attractors in which 
        # (1) the 5th ranked gene has score >= mi_threshold (default: 0.4) AND
        # (2) the score difference between the 1st and 2nd genes < 0.2
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