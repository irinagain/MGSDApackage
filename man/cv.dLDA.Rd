\name{cv.dLDA}
\alias{cv.dLDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validation for MGSDA
}
\description{
Chooses optimal tuning parameter lambda for function dLDA based on the m-fold cross-validation mean squared error
}
\usage{
cv.dLDA(Xtrain, Ytrain, lambdaval = NULL, nl = 100, msep = 5, eps = 1e-6,
    l_min_ratio = ifelse(n<p,0.1,0.0001),myseed=NULL,prior=TRUE,rho=1)
}
\arguments{
  \item{Xtrain}{
    A Nxp data matrix; N observations on the rows and p features on the columns
}
  \item{Ytrain}{
A N vector containing the group labels. Should be coded as 1,2,...,G, where G is the number of groups
}
  \item{lambdaval}{
Optional user-supplied sequence of tuning parameters; the default value is NULL and \code{cv.dLDA} chooses its own sequence
}
  \item{nl}{
Number of lambda values; the default value is 50
}
  \item{msep}{
Number of cross-validation folds; the default value is 5
}
  \item{eps}{
Tolerance level for the convergence of the optimization algorithm; the default value is 1e-6
}
  \item{l_min_ratio}{
Smallest value for lambda, as a fraction of \code{lambda.max}, the data-derived value for which all coefficients are zero; the default value is 0.1 if the number of samples \code{n} is larger than the number of variables \code{p}, and is 0.001 otherwise.
}
\item{myseed}{
Optional specification of random seed for generating the folds; the default value is NULL.
}
\item{prior}{
A logical indicating whether to put larger weights to the groups of larger size; the default value is TRUE.
}
\item{rho}{A scalar that ensures the objective function is bounded from below; the default value is 1.
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
 \item{lambdaval}{The sequence of tuning parameters used}
 \item{error_mean}{The mean cross-validated number of misclassified observations - a vector of length \code{length(lambdaval)}}
 \item{error_se}{The standard error associated with each value of \code{error_mean}}
 \item{lambda_min}{The value of tuning parameter that has the minimal mean cross-validation error}
 \item{f}{The mean cross-validated number of non-zero features - a vector of length \code{length(lambdaval)}}
}
\references{
I.Gaynanova, J.Booth and M.Wells (2016). "Simultaneous sparse estimation of canonical vectors in the p>>N setting", JASA, 111(514), 696-706.
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
n=10
p=100
G=3
ytrain=rep(1:G,each=n)
set.seed(1)
xtrain=matrix(rnorm(p*n*G),n*G,p)
# find optimal tuning parameter
out.cv=cv.dLDA(xtrain,ytrain)
# find V
V=dLDA(xtrain,ytrain,lambda=out.cv$lambda_min)
# number of non-zero features
sum(rowSums(V)!=0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
