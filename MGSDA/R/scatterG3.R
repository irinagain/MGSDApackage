scatterG3 <-
function(X,Y,V,k=NULL){
  require(ggplot2)
  #check that the data is in correct format and has right dimensions
  #
  ###############################
  data=data.frame(score=X%*%V,group=as.factor(Y))
  G=length(levels(data$group))
  if (is.null(k)){ k=G-1
  } else k=max(min(k,G-1),1)
  
  if (k>3){
  pairs(data[,1:k],col=Y,pch=as.character(Y), main=paste("Scatterplot of ",k," scores based on ",sum(rowSums(V)!=0) ," selected features",sep=""))
  } else if (k==2){
    plot(data[,1],data[,2],col=Y,pch=as.character(Y), main=paste("Scatterplot of ",k," scores based on ",sum(rowSums(V)!=0) ," selected features",sep=""),xlab="score.1",ylab="score.2")
  } else histG2(X,Y,V[,1])
}
