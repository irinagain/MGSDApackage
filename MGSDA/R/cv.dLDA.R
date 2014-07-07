cv.dLDA<-function(Xtrain,Ytrain,lambdaval=NULL,nl=100,msep=5,eps=1e-6,l_min_ratio=0.01){
  if (any(is.na(Xtrain))|any(is.na(Ytrain))) 
    stop("Missing values are not allowed!")
  
  n=length(Ytrain)
  if (nrow(Xtrain)!=n){
    stop(paste("Number of observations in Ytrain (",n,") doesn't match the number of rows in Xtrain, (",nrow(Xtrain),")",sep=""))
  }
  G=max(Ytrain)
  p=ncol(Xtrain)
  
  if (G==2){
    # two-group case, use glmnet
    n1=sum(Ytrain==1)
    n2=sum(Ytrain==2)
    Ynew=Ytrain
    Ynew[Ytrain==1]=-n/n1
    Ynew[Ytrain==2]=n/n2
    
    #require(glmnet)
    
    #calculate lambda path
    if (is.null(lambdaval)){
      out=glmnet(Xtrain,Ynew,family="gaussian",lambda=lambdaval,alpha=1,standardize=T)
      lambdaval=out$lambda
    }
  
    nl=length(lambdaval)
    error=matrix(0,msep,nl) #msep by lambda[i]
    features=matrix(p,msep,nl)
    
    #split the dataset into msep parts for each type
    id=1:n
    for (i in 1:G){
      id[Ytrain==i] <- sample(rep(seq_len(msep),length.out=sum(Ytrain==i)))
    }
    
    cat("Fold")
    #for each split of the dataset
    for (i in 1:msep){   
      cat(i)
      # determine test set and train set
      xtrain=Xtrain[id!=i,]
      ytrain=Ynew[id!=i]
      xtest=Xtrain[id==i,]
      ytest=Ytrain[id==i]
      
      #call glmnet once for all those lambda
      out=glmnet(xtrain,ytrain,family='gaussian',alpha=1,lambda=lambdaval,standardize=T)
      features[i,]=out$df
      
      #calculate the error rate
      for (j in 1:nl){  
       ypred=classifyV(xtrain,Ytrain[id!=i],xtest,out$beta[,j])
       error[i,j]=sum(ypred!=ytest) #how many features are selected
      }
    }
    #calculate mean error for all lambda
    errormean=colMeans(error)
    j=which.min(errormean)
    obj<-list(lambda=lambdaval[j],error=colMeans(error),f=round(colMeans(features)),lambdaval=lambdaval)
    return(obj) 
  } else{
    #multiple group case
    #calculate l_max
    D=.constructD(scale(Xtrain),Ytrain)
    l_max=max(sqrt(rowSums(D^2)))
    rm(D)
    
    #calculate lambda path
    if (!is.null(lambdaval)){
      lambdaval=lambdaval[lambdaval<=l_max]
      nl=length(lambdaval)
      if (nl<2) stop("There should be at least two lambdas")
    }
    else {
      lambdaval=10^seq(log10(l_min_ratio*l_max),log10(l_max),length.out=nl)
    }
    lambdaval=sort(lambdaval,decreasing=T)
    
    #keep errors and features
    error=matrix(0,msep,nl) #msep by lambda[i]
    features=matrix(p,msep,nl)
  
    #split the dataset into msep parts for each type
    id=1:n
    for (i in 1:G){
      id[Ytrain==i] <- sample(rep(seq_len(msep),length.out=sum(Ytrain==i)))
    }
  
    cat("Fold")
    #for each split of the dataset
    for (i in 1:msep){   
      cat(i)
      # determine test set and train set
      xtrain=Xtrain[id!=i,]
      Xadj=scale(xtrain)
      coef=attr(Xadj,which="scaled:scale")
    
      ytrain=Ytrain[id!=i]
      xtest=Xtrain[id==i,]
      ytest=Ytrain[id==i]
      
      #calculate D only once for all those lambdaval
      D=.constructD(Xadj,ytrain)
      Tot=crossprod(Xadj)/(length(ytrain)-1)
      V=matrix(0,p,G-1)
      error[i,1:nl]=length(ytest)
    
      #make discrimination with lambda[j],calculate the error rate
      for (j in 1:nl){  
        # no scaling for initial value
        V=.solveVcoordf2(Tot,D,lambdaval[j],eps=eps,V=V)
   
        features[i,j]=sum(rowSums(V)!=0) #how many features are selected
      
        if (features[i,j]>p-1){ #get to the end of the path when all the features are selected
          ytestpred=classifyV(xtrain,ytrain,xtest,diag(1/coef)%*%V)
          error[i,j:nl]=sum(ytestpred!=ytest)
          break
        }else if (features[i,j]>G-1){
          ytestpred=classifyV(xtrain,ytrain,xtest,diag(1/coef)%*%V)
          error[i,j]=sum(ytestpred!=ytest)
       }
    }
  }
  #calculate mean error for all lambda
  errormean=colMeans(error)
  j=which.min(errormean)
  obj<-list(lambda=lambdaval[j],error=colMeans(error),f=round(colMeans(features)),lambdaval=lambdaval)
  return(obj) 
  }
}
