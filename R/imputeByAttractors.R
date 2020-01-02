#' Imputation based on co-expressed gene attractor signatures
#' 
#' This function imputes dropout-suspected zeros by borrowing 
#' information from co-expressed genes in the identified 
#' co-expressed attractor signatures.
#' 
#' For each gene of interest, we use the top 20 genes involved 
#' in the same signature as the predictive variables to fit a 
#' linear regression model using cells where all those genes are
#' expressed. 
#' 
#' 
#' @param exprs A log-transformed expression matrix, where
#' rows and columns represent genes and cells, respectively.
#' 
#' @param dropout_ind A matrix where each row contains the 
#' indices of zero entries (first column for gene name, 
#' second column for cell ID) that are suspected as dropouts.
#' 
#' @param attractor_list A list where each element contains
#' the identified co-expressed signature. 
#' 
#' @param mi_threshold MI threshold (between 0 to 1) to determine the 
#' gene membership of signatures. Genes with greater than or equal to 
#' the threshold are considered to be associated with the signature. 
#' Default is 0.2.
#' 
#' 
#' @return A list consisting of the following components:
#' \item{imputed}{The imputed version of the expression matrix}
#' \item{dropout_ind}{An updated indices matrix of dropout-suspected
#' zero entries, where in each contains the indices of the remaining
#' dropouts after this step of imputation.}
#' 
#' 
#' 
#' @export


imputeByAttractors <- function(exprs, dropout_ind, attractor_list, mi_threshold = 0.4, verbose = TRUE){
  
  imputed <- exprs
  
  cells_all = colnames(exprs)

  N = length(attractor_list)
  
  if(N > 1){
    
    attractor_genes <- sapply(1:N, function(x) names(attractor_list[[x]][which(attractor_list[[x]] > mi_threshold)]))
    names(attractor_genes) <- names(attractor_list)
    
  }else{
    
    attractor_genes <- attractor_list
    attractor_genes[[1]] <- names(which(attractor_genes[[1]] > mi_threshold))
    
  }
  
  genes_to_impute <- unlist(attractor_genes)
  
  if(verbose)
    cat(length(unique(genes_to_impute)), 'genes are involved in attractors. \n')
  
  ## Assign gene to the most associated attractor if it is involved in multiple attractors
  genes_dup = names(which(table(genes_to_impute) > 1))
  
  if(length(genes_dup) > 0){
    
    for(g in genes_dup){
      
      indices = grep(g, attractor_genes)
      MI = sapply(indices, function(x) attractor_list[[x]][g])
      index = indices[which.max(MI)]
      index_to_remove = setdiff(indices, index)
      
      for(i in index_to_remove){
        j = grep(g, attractor_genes[[i]])
        attractor_genes[[i]] = attractor_genes[[i]][-j]        
      }
      
    }
    
  }
  
  # Filter out attractors in which there are < 5 genes left
  L = sapply(attractor_genes, length)
  attractor_genes <- attractor_genes[which(L >= 5)]
  N = length(attractor_genes)

  
  ## Impute genes involved in the attractor using the top 20 genes
  genes_imputed <- NULL
  
  for(a in 1:N){
    
    metagene = attractor_genes[[a]]
    metagene_members = metagene[1:min(20, length(metagene))]
    
    for(g in metagene){
      
      cells_zero = dropout_ind[which(dropout_ind[,1] == g),2]
      
      # remove genes that are not expressed in any cell of 'cells_zero'
      if(length(cells_zero) > 1){
        metagene_members_pos <- c(g, metagene_members[rowSums(imputed[metagene_members, cells_zero] == 0) < length(cells_zero)])
      }else if(length(cells_zero) == 1){
        metagene_members_pos <- c(g, metagene_members[imputed[metagene_members, cells_zero] > 0])
      }else{
        next
      }

      df <- as.data.frame(t(imputed[metagene_members_pos, ]))
      
      # use only cells in which all the genes in the attractor are expressed
      pos_id <- sapply(metagene_members_pos, function(x) which(imputed[x,] > 0), simplify = F)
      cells_subset <- cells_all[purrr::reduce(pos_id, intersect)]
      
      # if the eligible cells are too few, we reduce the number of metagene members and repeat the process above
      n = length(metagene_members_pos)-1
      while(length(cells_subset) < (n*10) && length(cells_subset) < 50 && n > 1){
        
        metagene_members_pos <- metagene_members_pos[1:n]
        df <- as.data.frame(t(imputed[metagene_members_pos, ]))
        
        pos_id <- sapply(metagene_members_pos, function(x) which(imputed[x,] > 0), simplify = F)
        cells_subset <- cells_all[purrr::reduce(pos_id, intersect)]
        n = length(metagene_members_pos)-1
        
      }
      
      
      if(n >= 1){
        
        df1 <- df[c(cells_zero, cells_subset),]
        for(x in metagene_members_pos){
          df1[cells_subset, x] <- capping(df1[cells_subset, x])
        }
        
        xnam <- paste0('x', 1:(length(metagene_members_pos)-1))
        colnames(df1) <- c('y', xnam)
        fo <- as.formula(paste('y ~ 0 + ', paste(xnam, collapse = '+')))
        fit <- lm(fo, df1[-(1:length(cells_zero)),])
        
        predicted <- predict(fit, df1[cells_zero,])
        predicted[predicted < 0] <- 0
        imputed[g, cells_zero] <- predicted
        
        genes_imputed <- c(genes_imputed, g)
        
      }
      
    }
    
  }
  
  dropout_ind1 <- dropout_ind[which(imputed[dropout_ind] == 0),]
  results = list(imputed = imputed, dropout_ind = dropout_ind1)
  return(results)
  
}