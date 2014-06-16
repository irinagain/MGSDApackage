.constructCw <-
function(ytrain){
  
  tmp=sort(ytrain,index=T)
  
  G=max(ytrain)
  ngroup=rep(0,G)
  for (i in 1:G){
    ngroup[i]=sum(ytrain==i)
  }
  s=cumsum(ngroup)
  
  Cw=matrix(0,s[G],s[G])
  
  for (i in 1:G){
    if (i==1) {Cw[1:ngroup[i],1:ngroup[i]]=matrix(1,ngroup[i],ngroup[i])/ngroup[i]
    } else Cw[(s[i-1]+1):s[i],(s[i-1]+1):s[i]]=matrix(1,ngroup[i],ngroup[i])/ngroup[i]
  }
  
  Wtmp=(diag(rep(1,s[G]))-Cw)/(s[G]-1)
  Wtmp[tmp$ix,tmp$ix]
}
