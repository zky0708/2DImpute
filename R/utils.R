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
