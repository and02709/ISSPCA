#' This function decomposes a matrix using singular value decomposition.  
#' The resulting eigenvector corresponding to the loadings in then shrunk
#' using an L1 penalty.
#' @param X This contains the training feature space
#' @param Y This contains the response vector
#' @param npc This is the number principal components desired
#' @param decomposition.type Determines which decomposition method is used.
#' @param nobs This contains the number of observations
#' @param p This contains the number of predictors
#' @param L This is the response kernel matrix
#' @param H This is the centering matrix
#' @keywords supervised PCA
#' @export
#' @examples SMD(x, sumabsu, sumabsv, niter=20,trace=TRUE, v, upos, uneg, vpos, vneg)
#' 

SPCA <- function(X, Y, npc, decomposition.type=c("default","RSpectra"), nobs, p, L, H){
  #nobs <- nrow(X)
  #p <- ncol(X)
  
  #L <- respkernel(Y, Ykernel, nobs)
  
  #H <- diag(1, nobs) - 1/nobs*rep(1, nobs)%*%t(rep(1, nobs))
  
  if(p > nobs){
    Eigendecomp <- eigen(L)
    # Need to retain eigenvectors to reconstruct Delta
    U <- Eigendecomp$vectors
    # Eigenvalues of kernel response matrix L
    EV <- Eigendecomp$values
    # Generation of diagonal matrix of square root eigenvalues of L
    Sigmat <- diag(zapsmall(sqrt(zapsmall(EV))))
    # Generation of Delta such that t(Delta)%*%Delta=L
    Delta <- U%*%Sigmat%*%t(U)
    # generate identity matrix that has dimensions equal to the training set
    Psi <- t(Delta)%*%H%*%X
    # generate Dual Supervised form
    M <- Psi%*%t(Psi)
    if(decomposition.type=="RSpectra"){
      decomp.obj <- RSpectra::eigs_sym(M, k=npc, which="LA")
    }
    else{
      decomp.obj <- eigen(M, only.values = F)
    }
    
    if(npc==1){
      Sigg <- sqrt(decomp.obj$values[1:npc])
      v <- t(Psi)%*%decomp.obj$vectors[,1:npc]%*%t(solve(Sigg))
    }
    else{
      Sigg <- diag(sqrt(decomp.obj$values[1:npc]))
      v <- t(Psi)%*%decomp.obj$vectors[,1:npc]%*%t(solve(Sigg))
    }
    return(list(e.values=decomp.obj$values[1:npc], v=v))
  }
  else{
    M <- t(X)%*%H%*%L%*%H%*%X
    
    if(decomposition.type=="RSpectra"){
      decomp.obj <- RSpectra::eigs_sym(M, k=npc, which="LA")
    }
    else{
      decomp.obj <- eigen(M, only.values = F)
    }
    return(list(e.values= decomp.obj$values[1:npc], v=decomp.obj$vectors[,1:npc]))
  }
}