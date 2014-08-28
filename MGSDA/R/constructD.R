.constructD <-function(X,Y){
  G=max(Y)
  p=ncol(X)  

  ngroup=rep(0,G)
  for (i in (1:G)){
    ngroup[i]=sum(Y==i)
  }
  s=cumsum(ngroup)
  
  D=matrix(0,p,G-1)
  for (i in 1:(G-1)){
    # more than 1 obs in group i+1
    if (ngroup[i+1]>1) D[,i]=sqrt(ngroup[i+1])*(colSums(X[Y<=i,])-s[i]*colMeans(X[Y==i+1,]))/(sqrt(s[i]*s[i+1]))
    # only 1 obs in group i+1
    else D[,i]=sqrt(ngroup[i+1])*(colSums(X[Y<=i,])-s[i]*X[Y==i+1,])/(sqrt(s[i]*s[i+1]))
  }
  
  D/sqrt(s[G]-1)
}
