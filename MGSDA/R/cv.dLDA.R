cv.dLDA<-function(Xtrain,Ytrain,lambdaval=NULL,nl=100,msep=5,eps=1e-6,l_min_ratio=0.01,myseed=NULL,prior=TRUE){
  if (any(is.na(Xtrain))|any(is.na(Ytrain))) 
    stop("Missing values are not allowed!")
  
  n=length(Ytrain)
  if (nrow(Xtrain)!=n){
    stop(paste("Number of observations in Ytrain (",n,") doesn't match the number of rows in Xtrain, (",nrow(Xtrain),")",sep=""))
  }
  G=max(Ytrain)
  p=ncol(Xtrain)
  
  if (G==2){
    n1=sum(Ytrain==1)
    n2=sum(Ytrain==2)
    Ynew=Ytrain
    Ynew[Ytrain==1]=-n/n1
    Ynew[Ytrain==2]=n/n2

    #calculate lambda path
    if (is.null(lambdaval)){
      out=glmnet(Xtrain,Ynew,family="gaussian",lambda=lambdaval,alpha=1,standardize=T)
      lambdaval=out$lambda
    }
  
    nl=length(lambdaval)
    error=matrix(0,msep,nl)
    features=matrix(p,msep,nl)
    
    if (!is.null(myseed)){set.seed(myseed)}
    id=1:n
    for (i in 1:G){
      id[Ytrain==i] <- sample(rep(seq_len(msep),length.out=sum(Ytrain==i)))
    }
    
    cat("Fold")
    for (i in 1:msep){   
      cat(i)
      xtrain=Xtrain[id!=i,]
      ytrain=Ynew[id!=i]
      xtest=Xtrain[id==i,]
      ytest=Ytrain[id==i]
      
      out=glmnet(xtrain,ytrain,family='gaussian',alpha=1,lambda=lambdaval,standardize=T)
      features[i,]=out$df
      
      for (j in 1:nl){  
       ypred=classifyV(xtrain,Ytrain[id!=i],xtest,out$beta[,j],prior=prior)
       error[i,j]=sum(ypred!=ytest) #how many features are selected
      }
    }
    errormean=colMeans(error)
    j=which.min(errormean)
    obj<-list(lambda=lambdaval[j],error=colMeans(error),f=round(colMeans(features)),lambdaval=lambdaval)
    return(obj) 
  } else{
    #multiple group case
    D=.constructD(scale(Xtrain),Ytrain)
    l_max=max(sqrt(rowSums(D^2)))
    rm(D)
    
    if (!is.null(lambdaval)){
      lambdaval=lambdaval[lambdaval<=l_max]
      nl=length(lambdaval)
      if (nl<2) stop("There should be at least two lambdas")
    }
    else {
      lambdaval=10^seq(log10(l_min_ratio*l_max),log10(l_max),length.out=nl)
    }
    lambdaval=sort(lambdaval,decreasing=T)
    
    error=matrix(0,msep,nl) 
    features=matrix(p,msep,nl)
  
    if (!is.null(myseed)){set.seed(myseed)}
    id=1:n
    for (i in 1:G){
      id[Ytrain==i] <- sample(rep(seq_len(msep),length.out=sum(Ytrain==i)))
    }
  
    cat("Fold")
    for (i in 1:msep){   
      cat(i)
      xtrain=Xtrain[id!=i,]
      Xadj=scale(xtrain)
      coef=attr(Xadj,which="scaled:scale")
      mtrain=attr(Xadj,which="scaled:center")
      ytrain=Ytrain[id!=i]
      xtest=Xtrain[id==i,]
      xtest=scale(xtest,center=mtrain,scale=coef) 
      ytest=Ytrain[id==i]
      
      D=.constructD(Xadj,ytrain)
      Tot=crossprod(Xadj)/(length(ytrain)-1)
      V=matrix(0,p,G-1)
      error[i,1:nl]=length(ytest)
    
      for (j in 1:nl){  
        V=.solveVcoordf2(Tot,D,lambdaval[j],eps=eps,V=V)
        features[i,j]=sum(rowSums(V)!=0) 
      
        if (features[i,j]>p-1){
          ytestpred=classifyV(Xadj,ytrain,xtest,V,prior=prior)
          error[i,j:nl]=sum(ytestpred!=ytest)
          break
        }else if (features[i,j]>0){
          ytestpred=classifyV(Xadj,ytrain,xtest,V,prior=prior)
          error[i,j]=sum(ytestpred!=ytest)
       }
    }
  }
  errormean=colMeans(error)
  j=which.min(errormean)
  obj<-list(lambda=lambdaval[j],error=colMeans(error),f=round(colMeans(features)),lambdaval=lambdaval)
  return(obj) 
  }
}
