#' Extreme values treatment
#' 
#' This function treats extreme values by winsorization. 
#' 
#' For values that lie outside the 1.5*IQR limits, we cap 
#' them by replacing those observations outside the lower
#' limit with the value of 5th %ile, and those that lie 
#' above the upper limit with the value of 95th %ile.
#' 
#' @param x A vector containing numeric numbers.
#' 
#' @return A vector with the same length as the input,
#' in which extreme values have been treated.
#' 
#' 
#' @export


capping <- function(x){

  qnt <- quantile(x, probs = c(.25, .75), na.rm = T)
  caps <- quantile(x, probs = c(.05, .95), na.rm = T)
  H <- 1.5*IQR(x, na.rm = T)
  x[x < (qnt[1] - H)] <- caps[1]
  x[x > (qnt[2] + H)] <- caps[2]

  return(x)
}
