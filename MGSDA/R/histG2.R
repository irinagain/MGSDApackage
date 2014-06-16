histG2 <-
function(X,Y,V){
  require(ggplot2)
  #check that the data is in correct format and has right dimensions
  #
  G=max(Y)
  if (G!=2){stop("Y specifies more than 2 groups")}
  n=length(Y)
  if (nrow(X)!=length(Y)){stop("Dimensions of X and Y don't match")}
  if (ncol(X)!=length(V)){stop("Dimensions of X and V don't match")}
  ###############################
  V=as.matrix(V)
  data=data.frame(score=X%*%V,group=as.factor(Y))
  colnames(data)=c("score","group")
  
  ggplot(data,aes(x=score,fill=group))+geom_histogram(position="identity",alpha=0.8,binwidth=(max(data$score)-min(data$score))/length(Y))+ggtitle(paste("Histogram of scores based on ",sum(V!=0) ," selected features",sep=""))
}
