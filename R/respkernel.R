#' A function kernelizes the response data
#' @param Y Vector of resposne data
#' @param kernel Desired kernel for the response data
#' @param nobs number of observations in the dataset
#' @keywords response
#' @export
#' @examples
#' respkernel(Y, kernel, nobs)




respkernel <- function(Y, kernel, nobs){
  ####missing data
  Y[is.na(Y)] <- mean(Y, na.rm=T)
  
  if (kernel=="linear"){
    L <- tcrossprod(Y)
  } else if (kernel=="delta"){
    #nobs <- dim(Y)[[1]]
    yf <- as.factor(Y)
    L<-matrix(0,nobs,nobs)
    for (i in levels(yf)){
      tmp<-yf==i
      L[tmp,tmp]<-1
    } 
  } else {
    stop("Please select a valid kernel, linear kernel or delta kernel")
  } 
  return(L)
}
