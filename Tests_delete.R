# Assume G=2, using LASSO function to solve MGSDA
G=2
n=200
p=5000
ytrain=rep(1:G,each=n/G)
set.seed(1)
xtrain = matrix(rnorm(p*n),n,p)
# Lasso solver
Xadj=scale(xtrain)
coef=attr(Xadj,which="scaled:scale")
n=length(ytrain)
Z=matrix(0,n,G)
for (g in 1:G){
    Z[ytrain==g,g]=1
}
n1=sum(ytrain==1)
n2=sum(ytrain==2)
Ytilde=sqrt(n1*n2)*Z%*%c(1/n1,-1/n2)
Vlasso1=solveMyLasso_c(Xadj,Ytilde,lambda=0.15)
Vlasso1=diag(1/coef)%*%Vlasso1
# Usual dLDA function, now uses solveMyLasso_c
V=dLDA(xtrain,ytrain,lambda=0.15)
# The two should agree
plot(V,Vlasso1) #and they do

# Does cv function work? yes it does
outcv=cv.dLDA(xtrain,ytrain)

##################
# For G=3, check that my LASSO solver works
G=3
r=G-1
n=200
p=1000
ytrain=rep(1:G,each=n/G)
n=length(ytrain)
set.seed(1)
xtrain = matrix(rnorm(p*n),n,p)
Y=matrix(rnorm(r*n),n,r)
out=solveMyLassoF_c(xtrain,Y,lambda=0.01,maxiter=1000) # works
sum(rowSums(out)!=0)
# Compares solveMyLassoF_c with solveMyLasso
G=2
r=G-1
n=200
p=1000
ytrain=c(rep(1,50),rep(2,130),rep(3,20))
set.seed(1)
xtrain = matrix(rnorm(p*n),n,p)
Y=matrix(rnorm(r*n),n,r)
lambda=0.1
out1=solveMyLasso_c(xtrain,Y,lambda=lambda,eps=1e-6)
sum(out1!=0)
out=solveMyLassoF_c(xtrain,Y,lambda=lambda,maxiter=1000,eps=1e-6) #works
sum(rowSums(out)!=0) 
plot(out1,out) #agree
###################################################
## Need a funnction that creates Y
n=length(ytrain)
# Create matrix Z
Z=matrix(0,n,G)
for (g in 1:G){
    Z[ytrain==g,g]=1
}
# Create cumulative sums
cumsum=rep(sum(ytrain==1),G)
for (i in 2:G){
    cumsum[i]=cumsum[i-1]+sum(ytrain==i)
}

# Create matrix H
H=matrix(0,G,G-1)
for (i in 1:(G-1)){
    # initialize each column
    H[1:i,i]=sqrt(sum(ytrain==i+1)/(cumsum[i]*cumsum[i+1]))
    H[i+1,i]=-sqrt(cumsum[i]/(cumsum[i+1]*sum(ytrain==i+1)))
}

# Create Y
Y=sqrt(n)*(Z%*%H)
D=constructD(xtrain,ytrain)
newD=t(xtrain)%*%Y/n

#Old time versus new time
# For G=3, check that my LASSO solver works
G=3
r=G-1
n=200
p=1000
ytrain=rep(1:G,each=n/G)
n=length(ytrain)
set.seed(1)
xtrain = matrix(rnorm(p*n),n,p)
a=proc.time()
out=dLDA(xtrain,ytrain,lambda=0.05)
b=proc.time()-a
b
library(microbenchmark) # old code, ~512 milliseconds median, min ~424, max ~957
# new code, ~183 milliseconds median, min ~170, max ~455!!!! much faster
microbenchmark(
    out=dLDA(xtrain,ytrain,lambda=0.05)
    )
sum(rowSums(out)!=0)
