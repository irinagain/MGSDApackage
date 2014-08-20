#Irina Gaynanova
#Nov 25th, 2013

# Classification for the test set based on V
# Xtrain is n1 by p
# Xtest is n2 by p
# Ytrain is n1 by 1
# V is p by g-1

#V should be standardized so that V^tWV=I

classifyV<-function(Xtrain,Ytrain,Xtest,V,prior=T,tol1=1e-10){
  p=ncol(Xtrain)
  if (ncol(Xtest)!=p){
    stop("Dimensions of Xtrain and Xtest don't match!")
  }
  G=max(Ytrain)
  if (length(V)/(G-1)!=p){
    stop("Dimensions of Xtrain and V don't match!")
  }
  
  ntrain=nrow(Xtrain)
  if (length(Ytrain)!=ntrain){
    stop("Dimensions of Xtrain and Ytrain don't match!")
  }
  
  ntest=nrow(Xtest)
  Ytest=rep(0,ntest)
  
  V=as.matrix(V)

  if (G==2){
    trainproj=Xtrain%*%V
    testproj=Xtest%*%V
    
    means=matrix(0,2,1)
    for (i in 1:2){
      means[i,]=mean(trainproj[Ytrain==i,])
    }  
    Dis=matrix(testproj^2,ntest,2)-2*tcrossprod(testproj,means)+matrix(t(means^2),ntest,2,byrow=T)
    if (prior) Dis=Dis-matrix(2*log(c(sum(Ytrain==1),sum(Ytrain==2))/ntrain),ntest,2,byrow=T)
    Ytest=apply(Dis,1,which.min)   
    return(Ytest)
  } else{
    ######################################
    ######### G>2 ########################   
    #need this even if use lda since prevents me from null variables (hopefully)
    trainproj=Xtrain%*%V
    testproj=Xtest%*%V
    myg=as.factor(Ytrain)
    group.means=tapply(trainproj,list(rep(myg,ncol(V)),col(trainproj)),mean)
    A1=var(trainproj-group.means[myg,]) 
    tmp=eigen(A1,symmetric=T)
    if (min(tmp$values)>tol1){ V=V%*%tmp$vectors%*%diag(1/sqrt(tmp$values))
    }else { # V is low rank
        if (sum(tmp$values>tol1)>1){V=V%*%tmp$vectors[,tmp$values>tol1]%*%diag(1/sqrt(tmp$values[tmp$values>tol1]))
        }else {V=V%*%tmp$vectors[,tmp$values==max(tmp$values)]/sqrt(tmp$values[tmp$values==max(tmp$values)])}
    }
        
    trainproj=Xtrain%*%V
    testproj=Xtest%*%V
#     means=matrix(0,G,ncol(V))
#     ngroup=rep(0,G)
#     for (i in 1:G){
#       ngroup[i]=sum(Ytrain==i)
#       if ((ngroup[i]>1)&(ncol(V)>1)){ means[i,]=colMeans(trainproj[Ytrain==i,])
#       }else if (ncol(V)==1){
#           means[i,]=mean(trainproj[Ytrain==i,])
#       }else {means[i,]=trainproj[Ytrain==i,]}
#     }
#     
#     #currently have A as I
#     Dis=matrix(rowSums(testproj^2),ntest,G)-2*tcrossprod(testproj,means)+matrix(rowSums(means^2),ntest,G,byrow=T)
#     if (prior) Dis=Dis-matrix(2*log(ngroup/ntrain),ntest,G,byrow=T)

#use lda classification
    if (prior==T){
        outlda=lda(trainproj,grouping=Ytrain,tol=1e-16)
    } else{
        outlda=lda(trainproj,grouping=Ytrain,prior=rep(1/max(Ytrain),max(Ytrain)),tol=1e-16)
    }
    ypredlda=predict(outlda,testproj)
    return(ypredlda$class)
  }
  
  #Ytest=apply(Dis,1,which.min)   
  #return(Ytest)
}