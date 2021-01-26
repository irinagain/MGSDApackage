\name{classifyV}
\alias{classifyV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classification for MGSDA
}
\description{
Classify observations in the test set using the supplied matrix of canonical vectors V and the training set.
}
\usage{
classifyV(Xtrain, Ytrain, Xtest, V, prior = T,tol1=1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xtrain}{
A Nxp data matrix; N observations on the rows and p features on the columns.
}
  \item{Ytrain}{
A N vector containing the group labels. Should be coded as 1,2,...,G, where G is the number of groups.
}
  \item{Xtest}{
A Mxp data matrix; M test observations on the rows and p features on the columns.
}
  \item{V}{
A pxr matrix of canonical vectors that is used to classify observations.
}
  \item{prior}{
A logical indicating whether to put larger weights to the groups of larger size; the default value is TRUE.
}
 \item{tol1}{
 Tolerance level for the eigenvalues of \eqn{V^tWV}{V'WV}. If some eigenvalues are less than \code{tol}, the low-rank version of \code{V} is used for classification.
}
}
\details{
 For a new observation with the value x, the classification is performed based on the smallest Mahalanobis distance in the projected space:
 \deqn{\min_{1\le g \le G}(V^tx-Z_g)(V^tWV)^{-1}(V^tx-Z_g)}{
    min_g (V'x-Z_g)(V'WV)^(-1)(V'x-Z_g),}
where \eqn{Z_g}{Z_g} are the group-specific means of the training dataset in the projected space and \eqn{W}{W} is the sample within-group covariance matrix.   

If \code{prior=T}, then the above distance is adjusted by \eqn{-2\log\frac{n_g}{N}}{-2 log n_g/N}, where \eqn{n_g}{n_g} is the size of group g.
}
\value{
Returns a vector of length M with predicted group labels for the test set.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
I.Gaynanova, J.Booth and M.Wells (2016) "Simultaneous Sparse Estimation of Canonical Vectors in the p>>N setting.", JASA, 111(514), 696-706.
}
\author{
Irina Gaynanova
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
### Example 1
# generate training data
n=10
p=100
G=3
ytrain=rep(1:G,each=n)
set.seed(1)
xtrain=matrix(rnorm(p*n*G),n*G,p)
# find V
V=dLDA(xtrain,ytrain,lambda=0.1)
sum(rowSums(V)!=0)
# generate test data
m=20
set.seed(3)
xtest=matrix(rnorm(p*m),m,p)
# perform classification
ytest=classifyV(xtrain,ytrain,xtest,V)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
