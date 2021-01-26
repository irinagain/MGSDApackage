\name{dLDA}
\alias{dLDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the matrix of discriminant vectors using L_1 penalty on the rows
}
\description{
Solve Multi-Group Sparse Discriminant Anlalysis problem for the supplied value of the tuning parameter lambda.% The problem is formulated as
%\deqn{ \hat V (\lambda)=\arg\min_{V \in \mathbb {R}^{p \times (G-1)}}\frac12\Tr\left(V^tW V\right)+\frac12\|D^tV-I\|^2_F+\lambda\sum_{i=1}^{p}\|v_{i}\|_{2},}{V(\lambda)=\arg\min}
%where \eqn{W} is the within-group sample covariance matrix, \eqn{D} is the matrix of sample orthogonal %mean contrasts between the groups and \eqn{\lambda\sum_{i=1}^{p}\|v_{i}\|_{2}}{\lambda*penalty} is the sparsity penalty.
}
\usage{
dLDA(xtrain, ytrain, lambda, Vinit = NULL,eps=1e-6,maxiter=1000,rho=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xtrain}{
  A Nxp data matrix; N observations on the rows and p features on the columns.
}
  \item{ytrain}{
    A N-vector containing the group labels. Should be coded as 1,2,...,G, where G is the number of groups.
}
  \item{lambda}{
    Tuning parameter.
}
  \item{Vinit}{
    A px(G-1) optional initial value for the optimization algorithm; the default value is NULL.
}
\item{eps}{
    Tolerance level for the convergence of the optimization algorithm; the default value is 1e-6.
}
\item{maxiter}{
    Maximal number of iterations for the optimization algorithm; the default value is 1000.
}
\item{rho}{A scalar that ensures the objective function is bounded from below; the default value is 1.
}
}
\details{
Solves the following optimization problem:
\deqn{\min_V \frac12 Tr(V^tWV+\rho V^tDD^tV)-Tr(D^tV)+\lambda\sum_{i=1}^p\|v_i\|_2}{
min_V \{1/2*Tr(V'WV+rho*V'DD'V-2D'V)+lambda*sum(||v_i||_2)\}.
}
Here W is the within-group sample covariance matrix and D is the matrix of orthogonal contrasts between the group means, both are constructed based on the supplied values of \code{xtrain} and \code{ytrain}.

When \eqn{G=2}{G=2}, the row penalty reduces to vector L_1 penalty.
}
\value{
Returns a px(G-1) matrix of canonical vectors V.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
I.Gaynanova, J.Booth and M.Wells (2016) "Simultaneous Sparse Estimation of Canonical Vectors in the p>>N setting", JASA, 111(514), 696-706.
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
# Example 1
n=10
p=100
G=3
ytrain=rep(1:G,each=n)
set.seed(1)
xtrain=matrix(rnorm(p*n*G),n*G,p)
V=dLDA(xtrain,ytrain,lambda=0.1)
sum(rowSums(V)!=0) # number of non-zero rows
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
