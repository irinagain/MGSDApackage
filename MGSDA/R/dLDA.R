dLDA <- function(xtrain, ytrain, lambda, Vinit=NULL, eps=1e-6, maxiter=1000, rho=1){ 
    if (any(is.na(xtrain))|any(is.na(ytrain))){
      stop("Missing values are not allowed")
    }
    
    n <- length(ytrain)
    if (nrow(xtrain) != n){
        stop("Dimensions of xtrain and ytrain don't match!")
    } 
  
    if (any(apply(xtrain,2,sd) < 1e-13)){
      stop(paste("Some features have standard deviation less than 1e-13!", sep = ""))
    }
    G <- max(ytrain)
    if (!is.null(Vinit)){
      if ((nrow(Vinit) != ncol(xtrain))|(ncol(Vinit) != G-1)){
          stop("Supplied initial value for Vinit has wrong dimensions!")
        }
    }
  
    #center and scale X
    Xadj <- scale(xtrain)
    coef <- attr(Xadj, which = "scaled:scale")
    
    Ytilde <- .createY(ytrain)  
    if (rho != 1){
        Xadj <- Xadj + (sqrt(rho) - 1)*(Ytilde %*% crossprod(Ytilde,Xadj))/nrow(xtrain)
        Ytilde <- Ytilde/sqrt(rho)
    }
    if (G == 2){
        V <- .solveMyLasso_c(Xadj, Ytilde, lambda=lambda, eps=eps, maxiter=maxiter, binit=Vinit)
    }else{
        V <- .solveMyLassoF_c(Xadj, Ytilde, lambda=lambda, eps=eps, maxiter=maxiter, binit=Vinit) 
    }
  
  V[abs(V) < eps] <- 0
  diag(1/coef) %*% V
}
