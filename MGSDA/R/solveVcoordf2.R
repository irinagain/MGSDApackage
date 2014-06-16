.solveVcoordf2 <-
function(W,D,lambda,eps=10^(-6),V=NULL,maxiter=10000){

  p=nrow(W)
  r=ncol(D)
  
  if (r==1){
    #########################################
    #V is a vector, this is just lasso
    if (is.null(V)) V=rep(0,p)
    
    if (p==1){
      V=sign(D)*max(0,abs(D)-lambda)
      V
    }    
  } else{
    ##################################################
    #V is a matrix, group lasso
    if (is.null(V)) V=matrix(0,p,r)

    if (p==1){
      V=max(1-lambda/sqrt(D^2),0)*D
      V
    }
  }
  outc=.solveVcoordf2C(D,W,V,lambda,p,r,eps,maxiter)
  if (outc[[9]]>=maxiter) warning(paste("Convergence was not achieved for l=",lambda,sep=""))
  V=matrix(outc[[3]],p,r)
  V[abs(V)<eps]=0
  V
}
