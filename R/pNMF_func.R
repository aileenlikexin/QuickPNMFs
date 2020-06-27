#' Fast Projective Nonnegative Matrix Factorization Realizatiton based on euclidean distance / KL Divergence / Discriminant pNMF.
#'
#' @param X Input data matrix, where rows represent genes, columns represent cells.
#' @param rank Specification of the factorization rank.
#' @param tol A threshold below which would be considered converged.
#' @param maxIter Number of max iteration times.
#' @param verboseN A boolean value indicating whether to print number of iterations.
#' @param zerotol A threshold on basis loadings below which would be considered zero.
#' @param method A character string indicating which method to be used. One of "EucDist", "KL", or "DPNMF".
#' @param label A character vector indicating the cluster type for each cell. Required only when method = "DPNMF".
#' @param mu A numerical value.
#' @param lambda A numerical value.
#' @param seed random seed
#' 
#' @return A list of fitted NMF model, our pNMF basis, and mapped scores.
#' 
#' @importFrom irlba irlba
#' @export
#' 
#' @examples 
#' k <- matrix(rpois(600, lambda=5), nrow=20, ncol=30)
#' estim <- PNMFfun(k, rank=3, method="EucDist")
#'
#'
#'
PNMFfun <- function(X, rank=10, tol=1e-3, maxIter=500, verboseN=FALSE, zerotol=1e-10, method="KL", label=NULL, mu=1, lambda=0.01, seed=123) {
  #nmfmod <- NMF::nmf(X, rank)
  set.seed(seed)
  Init <- irlba(X, nv = rank)
  Winit <- Init$u
  Winit[Winit < 0] <- -Winit[Winit < 0]
  
  if (method == "EucDist") {
    W <- PNMF_EucDistC(X, Winit, tol, maxIter, verboseN, zerotol) 
    ld <- t(X) %*% W
  }
  else if (method == "KL") {
    W <- PNMF_KLC(X, Winit, tol, maxIter, verboseN, zerotol) 
    ld <- t(X) %*% W
  }
  else if (method == "DPNMF") {
    if (length(label) != dim(X)[2]) {
      stop("Cluster labels must have same length as number of cells.")
    }
    cluvec <- as.factor(label)
    cluvec.num <- as.numeric(cluvec)
    cluvec.ord <- order(cluvec.num)
    Xord <- X[, cluvec.ord]
    clunum <- as.integer(table(cluvec))
    
    W <- DPNMFC(X, Winit, tol, maxIter, verboseN, zerotol, Xord, clunum, mu, lambda) 
    ld <- t(X) %*% W
  }
  
  return(list(basis=W, coef=ld))
}







