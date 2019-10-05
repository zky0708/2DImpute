#' Example single-cell RNA-seq dataset
#' 
#' The dataset contains 3,774 genes and 338 cells, 
#' a subset of the 50\%:50\% Jurkat:293T cell mixture 
#' data downloaded from the 10X Genomics website.
#'
#' @docType data
#' 
#' @usage data(ge_10x_sample)
#'  
#' @references Zheng, G.X. et al. Massively parallel 
#' digital transcriptional profiling of single cells. 
#' \emph{Nat Commun} 8, 14049 (2017). (\href{https://www.ncbi.nlm.nih.gov/pubmed/28091601}{PubMed})
#'
#' @source \href{https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/jurkat:293t_50:50}{10X Genomics}.

"ge_10x_sample"

#' Example output of function \code{\link{findDropouts}}
#' 
#' The results of performing dropout identification on data 
#' \code{\link{ge_10x_sample}}. A list consisting of
#' \describe{
#'     \item{dropout_ind}{A matrix in which each row contains the indices of a dropout}
#'     \item{J}{Calculated pairwise Jaccard Index between cells}
#' }
#' 
#' @usage data(results1)

"results1"


#' Example output of function \code{\link{imputeByAttractors}}
#' 
#' The results of performing imputation using identified 
#' co-expressed attractors for data \code{\link{ge_10x_sample}}. 
#' A list consisting of 
#' \describe{
#'     \item{imputed}{The imputed version of the expression matrix}
#'     \item{dropout_ind}{An updated indices matrix of dropout-suspected}
#' }
#' 
#' @usage data(results2)

"results2"


#' Example output of function \code{\link{attractorFinding}}
#' 
#' An attractor list containing the two co-expressed gene signatures 
#' identified in \code{\link{ge_10x_sample}} by running function 
#' \code{\link{attractorFinding}}. 
#' 
#' Each signature is a vector consisting of the normalized mutual 
#' information between the corresponding gene and the converged metagene. 
#' 
#' @usage data(attractors)

"attractors"


#' Example output of function \code{\link{ddImpute}} 
#' 
#' The imputed version of the example scRNA-seq data
#' \code{\link{ge_10x_example}} by ddImpute algorithm.
#' 
#' @usage data(imputed)

"imputed"






