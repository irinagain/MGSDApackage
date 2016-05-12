cv.dLDA<-function(Xtrain,Ytrain,lambdaval=NULL,nl=100,msep=5,eps=1e-6,l_min_ratio=ifelse(n<p,0.1,0.0001),myseed=NULL,prior=TRUE,rho=0.01){
    
  if (any(is.na(Xtrain))|any(is.na(Ytrain))) 
    stop("Missing values are not allowed!")
  
  n <- length(Ytrain)
  if (nrow(Xtrain)!=n){
    stop(paste("Number of observations in Ytrain (",n,") doesn't match the number of rows in Xtrain, (",nrow(Xtrain),")",sep=""))
  }
  if (any(apply(Xtrain,2,sd)<1e-13)){
      stop(paste("Some features have standard deviation less than 1e-13!",sep=""))
  }
  
  G <- max(Ytrain)
  p <- ncol(Xtrain)
  
  if (!is.null(myseed)){
      set.seed(myseed)
  }
  id <- 1:n
  for (i in 1:G){
      id[Ytrain==i] <- sample(rep(seq_len(msep),length.out=sum(Ytrain==i)))
  }
  
  if (G==2){
      Z <- matrix(0,n,2)
      for (g in 1:2){
          Z[Ytrain==g,g]=1
      }
      n1 <- sum(Ytrain==1)
      n2 <- n-n1
      Ytilde <- sqrt(n1*n2)*Z%*%c(1/n1,-1/n2)
      if (rho!=1){
          X_w <- scale(Xtrain)+(sqrt(rho)-1)*(Ytilde%*%crossprod(Ytilde,scale(Xtrain)))/n
          Y_w <- Ytilde/sqrt(rho)
      }else{
          X_w <- scale(Xtrain)
      }

     #calculate lambda path
    l_max <- max(abs(crossprod(scale(Xtrain),Ytilde)))/n
    if (!is.null(lambdaval)){
        lambdaval <- lambdaval[lambdaval<=l_max]
        nl <- length(lambdaval)
        if (nl<2) stop("There should be at least two lambdas")
    }
    else {
        lambdaval <- 10^seq(log10(l_min_ratio*l_max),log10(l_max),length.out=nl)
    }
    lambdaval <- sort(lambdaval,decreasing=T)
  
    error <- matrix(1,msep,nl)
    features <- matrix(min(n,p),msep,nl)
      
    cat("Fold")
    for (i in 1:msep){   
      cat(i)
      xtrain <- Xtrain[id!=i,]
      Xadj <- scale(xtrain)
      coef <- attr(Xadj,which="scaled:scale")
      mtrain <- attr(Xadj,which="scaled:center")
      ytrain <- Ytrain[id!=i]
      xtest <- Xtrain[id==i,]
      xtest <- scale(xtest,center=mtrain,scale=coef) 
      ytest <- Ytrain[id==i]
      
      nnew <- length(ytrain)
      Znew <- Z[id!=i,]
      n1 <- sum(ytrain==1)
      n2 <- nnew-n1
      Ytilde <- sqrt(n1*n2)*Znew%*%c(1/n1,-1/n2)
      if (rho!=1){
          X_w <- Xadj+(sqrt(rho)-1)*(Ytilde%*%crossprod(Ytilde,Xadj))/nnew
          Y_w <- Ytilde/sqrt(rho)
      }else{
          X_w <- Xadj
      }
      
      V <- matrix(0,p,1)

      for (j in 1:nl){  
        V <- .solveMyLasso_c(X_w,Y_w,lambda=lambdaval[j],binit=V)
        ypred <- classifyV(Xadj,ytrain,xtest,V,prior=prior)
        error[i,j] <- sum(ypred!=ytest)/length(ytest) #how many features are selected
        features[i,j] <- sum(V!=0)
        if (features[i,j]>=min(n,p)){
            break
        }
      }
    }
    errormean <- colMeans(error)
    j <- which.min(errormean)
    obj <- list(lambda=lambdaval[j],error=colMeans(error),f=round(colMeans(features)),lambdaval=lambdaval)
    return(obj) 
  } else{
    #multiple group case
    Ytilde <- .createY(Ytrain)
    if (rho!=1){
        X_w <- scale(Xtrain)+(sqrt(rho)-1)*(Ytilde%*%crossprod(Ytilde,scale(Xtrain)))/n
        Y_w <- Ytilde/sqrt(rho)
    }else{
        X_w <- scale(Xtrain)
    }
    l_max <- max(abs(crossprod(scale(Xtrain),Ytilde)))/n
    if (!is.null(lambdaval)){
        lambdaval <- lambdaval[lambdaval<=l_max]
        nl <- length(lambdaval)
        if (nl<2) stop("There should be at least two lambdas")
    }else{
        lambdaval <- 10^seq(log10(l_min_ratio*l_max),log10(l_max),length.out=nl)
    }
    lambdaval <- sort(lambdaval,decreasing=T)
    
    error <- matrix(1,msep,nl) 
    features <- matrix(min(n,p),msep,nl)
  
    cat("Fold")
    for (i in 1:msep){   
      cat(i)
      xtrain <- Xtrain[id!=i,]
      Xadj <- scale(xtrain)
      coef <- attr(Xadj,which="scaled:scale")
      mtrain <- attr(Xadj,which="scaled:center")
      ytrain <- Ytrain[id!=i]
      xtest <- Xtrain[id==i,]
      xtest <- scale(xtest,center=mtrain,scale=coef) 
      ytest <- Ytrain[id==i]
      
      Ytilde <- .createY(ytrain)
      if (rho!=1){
          X_w <- Xadj+(sqrt(rho)-1)*(Ytilde%*%crossprod(Ytilde,Xadj))/nrow(Xadj)
          Y_w <- Ytilde/sqrt(rho)
      }else{
          X_w <- Xadj
      }
      V <- matrix(0,p,G-1)
    
      for (j in 1:nl){  
        V <- .solveMyLassoF_c(X_w,Y_w,lambda=lambdaval[j],binit=V)
        features[i,j] <- sum(rowSums(V)!=0) 
      
        if (features[i,j]>min(n,p-1)){
          ytestpred <- classifyV(Xadj,ytrain,xtest,V,prior=prior)
          error[i,j] <- sum(ytestpred!=ytest)/length(ytest)
          break
        }else if (features[i,j]>0){
          ytestpred <- classifyV(Xadj,ytrain,xtest,V,prior=prior)
          error[i,j] <- sum(ytestpred!=ytest)/length(ytest)
       }
    }
  }
  errormean <- colMeans(error)
  j <- which.min(errormean)
  obj <- list(lambda=lambdaval[j],error=colMeans(error),f=round(colMeans(features)),lambdaval=lambdaval)
  return(obj) 
  }
}
