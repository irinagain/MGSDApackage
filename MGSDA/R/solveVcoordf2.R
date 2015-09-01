.solveVcoordf2 <-
function(W,D,lambda,eps=10^(-6),V=NULL,maxiter=10000){

  p=nrow(W)
  r=ncol(D)
  
  if (r==1){
    #V is a vector
    if (is.null(V)) V=rep(0,p)
    
    if (p==1){
      V=sign(D)*max(0,abs(D)-lambda)
      V
    }    
  } else{
    #V is a matrix
    if (is.null(V)) V=matrix(0,p,r)

    if (p==1){
      V=max(1-lambda/sqrt(D^2),0)*D
      V
    }
  }
  niter=0;
  outc=.C("solveVcoordf2", as.double(as.vector(D)),as.double(as.vector(W)),as.double(as.vector(V)),as.double(lambda),as.integer(p),as.integer(r),as.double(eps),as.integer(maxiter),as.integer(niter))
  if (outc[[9]]>=maxiter) warning(paste("Convergence was not achieved for l=",lambda,sep=""))
  V=matrix(outc[[3]],p,r)
  V[abs(V)<eps]=0
  V
}
