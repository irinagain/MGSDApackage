dLDA <-function(xtrain,ytrain,lambda,Vinit=NULL){ 
  if (nrow(xtrain)!=length(ytrain)){
    stop("Dimensions of xtrain and ytrain don't match!")
  } 
  
  if (any(is.na(xtrain))|any(is.na(ytrain))){
    stop("Missing values are not allowed")
  }
  
  fsd=apply(xtrain,2,sd)
  if (any(fsd)<1e-13){
      stop(paste("Some features have standard deviation less than 1e-13!",sep=""))
  }
  #center and scale X
  Xadj=scale(xtrain)
  coef=attr(Xadj,which="scaled:scale")
  
  G=max(ytrain)
  
  if (G==2){
    n1=sum(ytrain==1)
    n2=sum(ytrain==2)
    Ynew=ytrain
    Ynew[ytrain==1]=-(n1+n2)/n1
    Ynew[ytrain==2]=(n1+n2)/n2
    
    out=glmnet(Xadj,Ynew,family="gaussian",alpha=1,standardize=F,lambda=lambda)
    V=as.matrix(out$beta)
    
  } else {   
    D=.constructD(Xadj,ytrain)
  
    if (is.null(Vinit)) {
      V=.solveVcoordf2(W=crossprod(Xadj)/(length(ytrain)-1),D=D,lambda=lambda)
    } else {
      if ((nrow(Vinit)!=ncol(xtrain))|(ncol(Vinit)!=G-1)){
        stop("Supplied initial value for Vinit has wrong dimensions!")
      }
      V=.solveVcoordf2(W=crossprod(Xadj)/(length(ytrain)-1),D=D,lambda=lambda,V=Vinit)
    }
  }
  
  diag(1/coef)%*%V
}
