cv.dLDA<-function(Xtrain,Ytrain,lambdaval=NULL,nl=100,msep=5,eps=1e-6,l_min_ratio=ifelse(n<p,0.1,0.0001),myseed=NULL,prior=TRUE,rho=1){
    
  if (any(is.na(Xtrain))|any(is.na(Ytrain))){ 
    stop("Missing values are not allowed!")
  }      
  n <- length(Ytrain)
  if (nrow(Xtrain)!=n){
    stop(paste("Number of observations in Ytrain (",n,") doesn't match the number of rows in Xtrain, (",nrow(Xtrain),")",sep=""))
  }
  if (any(apply(Xtrain,2,sd)<1e-13)){
      stop(paste("Some features have standard deviation less than 1e-13!",sep=""))
  }
  
  G <- max(Ytrain)
  p <- ncol(Xtrain)
  
  # Set seed for the assignment to folds
  if (!is.null(myseed)){
      set.seed(myseed)
  }
  
  # Divide samples into folds
  id <- 1:n
  for (i in 1:G){
      id[Ytrain == i] <- sample(rep(seq_len(msep),length.out = sum(Ytrain == i)))
  }
  
  # Form matrix Ytilde and make adjustment if rho!=1
  Ytilde <- .createY(Ytrain)
  if (rho != 1){
      X_w <- scale(Xtrain)+(sqrt(rho)-1)*(Ytilde%*%crossprod(Ytilde,scale(Xtrain)))/n
      Y_w <- Ytilde/sqrt(rho)
  }else{
      X_w <- scale(Xtrain)
      Y_w <- Ytilde
  }
  
  # Calculate lambda path
  l_max <- max(abs(crossprod(scale(Xtrain),Ytilde)))/n
  if (!is.null(lambdaval)){
      lambdaval <- lambdaval[lambdaval <= l_max]
      nl <- length(lambdaval)
      if (nl < 2) stop("There should be at least two lambdas")
  }else{
      lambdaval <- 10^seq(log10(l_min_ratio*l_max),log10(l_max),length.out = nl)
  }
  lambdaval <- sort(lambdaval,decreasing=T)
  
  # Store errors and the number of selected features
  error <- matrix(1, n, nl)
  features <- matrix(min(n,p),msep,nl)
  
    cat("Fold")
    for (i in 1:msep){   
      cat(i)
      xtrain <- Xtrain[id!=i,]
      Xadj <- scale(xtrain)
      coef <- attr(Xadj, which = "scaled:scale")
      mtrain <- attr(Xadj, which = "scaled:center")
      ytrain <- Ytrain[id!=i]
      xtest <- Xtrain[id==i,]
      xtest <- scale(xtest, center = mtrain, scale = coef) 
      ytest <- Ytrain[id==i]
      
      Ytilde <- .createY(ytrain)
      if (rho!=1){
          X_w <- Xadj+(sqrt(rho)-1)*(Ytilde%*%crossprod(Ytilde,Xadj))/nrow(Xadj)
          Y_w <- Ytilde/sqrt(rho)
      }else{
          X_w <- Xadj
          Y_w <- Ytilde
      }
      V <- matrix(0,p,G-1)
    
      for (j in 1:nl){
        if (G != 2){
            V <- .solveMyLassoF_c(X_w,Y_w,lambda=lambdaval[j],binit=V)
            features[i,j] <- sum(rowSums(V) != 0) 
        }else{
            V <- .solveMyLasso_c(X_w,Y_w,lambda=lambdaval[j],binit=V)
            features[i,j] <- sum(V != 0)
        }
        if (features[i,j]>min(n,p-1)){
          ytestpred <- classifyV(Xadj,ytrain,xtest,V,prior=prior)
          error[id == i, j] <- ytestpred != ytest
          break
        }else if (features[i,j]>0){
          ytestpred <- classifyV(Xadj,ytrain,xtest,V,prior=prior)
          error[id == i, j] <- ytestpred != ytest
       }
    } # end loop for lambda
  }# end loop for folds
    error_mean <- colMeans(error)
    error_se <- apply(error,2,sd)/sqrt(n)
    j <- which.min(error_mean)
    obj <- list(lambda_min = lambdaval[j],error_mean = error_mean, error_se = error_se, f = round(colMeans(features)),lambdaval = lambdaval)
    return(obj) 
}
