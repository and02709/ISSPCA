#' This function decomposes a matrix using singular value decomposition.  
#' The resulting eigenvector corresponding to the loadings in then shrunk
#' using an L1 penalty.
#' @param X The training feature space
#' @param Y The response vector
#' @param xts The testing feature space
#' @param sumabsv sparseness penalty
#' @param niter
#' @param npc number of principal components
#' @param npc.iter number of pls iterations
#' @param orth
#' @param trace
#' @param v
#' @param center
#' @param cnames
#' @param vpos
#' @param vneg
#' @param compute.pve
#' @param Ykernel
#' @param decomposition.type
#' @keywords parent function
#' @export
#' @examples SMD(x, sumabsu, sumabsv, niter=20,trace=TRUE, v, upos, uneg, vpos, vneg)

iSSPCA <- function(X, Y, xts, sumabsv=4, niter=20, npc=1, npc.iter=1, orth=TRUE, 
                   trace=FALSE, v=NULL, center=TRUE, cnames=NULL, 
                   vpos=FALSE, vneg=FALSE, compute.pve=TRUE, 
                   Ykernel=c("linear", "delta"), 
                   decomposition.type=c("default","RSpectra")){
  X <- as.matrix(X)
  xts <- as.matrix(xts)
  Y <- as.matrix(Y)
  
  nobs <- nrow(X)
  p <- ncol(X)
  
  L <- respkernel(Y, Ykernel, nobs)
  H <- diag(1, nobs) - 1/nobs*rep(1, nobs)%*%t(rep(1, nobs))
  Ii <- diag(p)
  
  if(center){
    xmeans <- colMeans(X)
    X <- t(apply(X, 1, function(x) x-xmeans))
    xts <- t(apply(xts, 1, function(x) x-xmeans))
  }
  
  i <- 1
  
  if(sumabsv==0){
    while(i < (npc.iter+1)){
      if(i==1){
        out <- SPCA(X=X, Y=Y, npc=npc, decomposition.type=decomposition.type, nobs=nobs, p=p, L=L, H=H)
        Ut <- as.matrix(out$v)
        W <- Ut
      }
      else{
        ti <- X%*%W
        Pi <- t(X)%*%ti%*%solve(t(ti)%*%ti)
        Qi <- Ii - Pi%*%solve(t(Pi)%*%Pi)%*%t(Pi)
        Xi <- X%*%Qi
        Psi <- t(Delta)%*%H%*%Xi
        out <- SPCA(X=Xi, Y=Y, npc=npc, decomposition.type=decomposition.type, nobs=nobs, p=p, L=L, H=H)
        Ut <- Re(as.matrix(out$v))
        W <- as.matrix(cbind(W,Ut))
      }
      
      i <- i+1
    }
  }
  else{
    while(i < (npc.iter+1)){
      if(i==1){
        Eigendecomp <- eigen(L)
        U <- Eigendecomp$vectors
        EV <- Eigendecomp$values
        Sigmat <- diag(sqrt(zapsmall(EV)))
        Delta <- zapsmall(U%*%Sigmat%*%t(U))
        Psi <- t(Delta)%*%H%*%X
        out <- PMDL1(Psi,sumabsu=sqrt(nrow(X)), sumabsv=sumabsv, niter=niter,K=npc,orth=orth,trace=trace,v=v,center=FALSE,cnames=cnames, upos=FALSE, uneg=FALSE, vpos=vpos, vneg=vneg)
        Ut <- as.matrix(out$v)
        W <- Ut
      }
      else{
        ti <- X%*%W
        Pi <- t(X)%*%ti%*%solve(t(ti)%*%ti)
        Qi <- Ii - Pi%*%solve(t(Pi)%*%Pi)%*%t(Pi)
        Xi <- X%*%Qi
        Psi <- t(Delta)%*%H%*%Xi
        out <- PMDL1(Psi,sumabsu=sqrt(nrow(X)), sumabsv=sumabsv, niter=niter,K=npc,orth=orth,trace=trace,v=v,center=FALSE,cnames=cnames, upos=FALSE, uneg=FALSE, vpos=vpos, vneg=vneg)
        Ut <- Re(as.matrix(out$v))
        W <- as.matrix(cbind(W,Ut))
      }
      
      i <- i+1
    }
  }
  
  
  Z <- X%*%W
  z <- xts%*%W
  
  return(list(V=W, Z=Z, z=z))
}