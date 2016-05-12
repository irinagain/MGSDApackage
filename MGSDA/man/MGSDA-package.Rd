\name{MGSDA-package}
\alias{MGSDA-package}
\alias{MGSDA}
\docType{package}
\title{
MGSDA: Multi-Group Sparse Discriminant Analysis
}
\description{
Implements Multi-Group Sparse Discriminant Analysis proposal of I.Gaynanova, J.Booth and M.Wells (2015) "Simultaneous sparse estimation of canonical vectors in the p>>n setting", JASA, to appear.
}
\details{
\tabular{ll}{
Package: \tab MGSDA\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2016-05-11\cr
License: \tab GPL(>=2)\cr
}
Three functions for estimation and prediction: \code{cv.dLDA}, \code{dLDA} and \code{classifyV}.

%Two functions for visualization: \code{histG2} and \code{scatterG3}.
}
\author{
Irina Gaynanova

Maintainer: Irina Gaynanova <irinag@stat.tamu.edu>
}
\references{
I.Gaynanova, J.Booth and M.Wells (2015). "Simultaneous sparse estimation of canonical vectors in the p>>N setting", JASA, to appear
}

%\keyword{package}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
### Example 1
# generate training data
n <- 10
p <- 100
G <- 3
ytrain <- rep(1:G,each=n)
set.seed(1)
xtrain <- matrix(rnorm(p*n*G),n*G,p)
# find matrix of canonical vectors V
V <- dLDA(xtrain,ytrain,lambda=0.1)
sum(rowSums(V)!=0)
# generate test data
m <- 20
set.seed(3)
xtest <- matrix(rnorm(p*m),m,p)
# perform classification
ytest <- classifyV(xtrain,ytrain,xtest,V)
}
