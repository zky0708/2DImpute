# 2DImpute
An imputation method that recovers dropout events in single cell RNA sequencing (scRNA-seq) data from relationships in two dimensions.

## Installation
You can install 2DImpute from GitHub with:
```R
# install.packages("devtools")
require(devtools)
install_github("weiyi-bitw/cafr")   # The 'cafr' R package is required for running 2DImpute.
install_github("zky0708/2DImpute")
```

## Quick start
`2DImpute` takes a log-transformed normalized scRNA-seq data matrix (genes in rows and cells in columns) as input, and outputs the imputed version with the same dimension.
1. Load the R package ("R2DImpute")
```R
require(R2DImpute)
```

2. The imputation task can be done with one single function `run_2DImpute`:
```R
res <- run_2DImpute(
  exprs,                      # The input expression matrix.
  t = 0.2,                    # Threshold set for dropout identification.
  genes = NULL,               # A vector of genes of which the zeros will get imputed. If it is NULL (default), all genes will be considered.
  k = 10,                     # Number of nearest neighbors.
  ncores = 1,                 # Number of cores to be used.
  return_J = FALSE,           # Whether to return the calculated pairwise Jaccard matrix between cells.
  return_attractors = FALSE,  # Whether to return the identified co-expressed gene attractor signatures.
  verbose = TRUE              # Whether to show the progress of imputation.
)
imputed_exprs <- res$imputed

```
For example:
```R
# An example 10x scRNA-seq dataset containing 3,775 genes and 338 cells
data(ge_10x_sample)          

# Fast demonstration (only imputes zeros in genes 'XIST' and 'CD3D')
res <- run_2DImpute(
  exprs = ge_10x_sample, 
  genes = c('XIST', 'CD3D'), 
  ncores = 2
  )
imputed_exprs <- res$imputed

# Return identified co-expressed attractor signatures
res <- run_2DImpute(
  exprs = ge_10x_sample, 
  genes = c('XIST', 'CD3D'), 
  ncores = 2, 
  return_attractors = TRUE
  )
imputed_exprs <- res$imputed
attractors <- res$attractors

# Full imputation
res <- run_2DImpute(
  exprs = ge_10x_sample, 
  ncores = 2
  )
imputed_exprs <- res$imputed
```

