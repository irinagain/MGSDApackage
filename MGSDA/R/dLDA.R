dLDA <-function(xtrain,ytrain,lambda,Vinit=NULL,eps=1e-6){ 
  if (nrow(xtrain)!=length(ytrain)){
    stop("Dimensions of xtrain and ytrain don't match!")
  } 
  
  if (any(is.na(xtrain))|any(is.na(ytrain))){
    stop("Missing values are not allowed")
  }
  
  fsd=apply(xtrain,2,sd)
  if (any(fsd<1e-13)){
      stop(paste("Some features have standard deviation less than 1e-13!",sep=""))
  }
  #center and scale X
  Xadj=scale(xtrain)
  coef=attr(Xadj,which="scaled:scale")
  
  G=max(ytrain)
  
  if (G==2){
    n=length(ytrain)
    Z=matrix(0,n,2)
    for (g in 1:2){
        Z[ytrain==g,g]=1
    }
    n1=sum(ytrain==1)
    n2=n-n1
    Ytilde=sqrt(n1*n2/n)*Z%*%c(1/n1,-1/n2)
    V=solveMyLasso_c(Xadj,Ytilde,lambda=lambda)
    
  } else {   
    D=.constructD(Xadj,ytrain)
  
    if (is.null(Vinit)) {
      V=.solveVcoordf2(W=crossprod(Xadj)/(length(ytrain)-1),D=D,lambda=lambda,eps=eps)
    } else {
      if ((nrow(Vinit)!=ncol(xtrain))|(ncol(Vinit)!=G-1)){
        stop("Supplied initial value for Vinit has wrong dimensions!")
      }
      V=.solveVcoordf2(W=crossprod(Xadj)/(length(ytrain)-1),D=D,lambda=lambda,V=Vinit,eps=eps)
    }
  }
  
  diag(1/coef)%*%V
}
