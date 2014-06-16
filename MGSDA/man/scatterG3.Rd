\name{scatterG3}
\alias{scatterG3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scatterplot of projected data for the multi-group case
}
\description{
Creates scatterplot of the data projected onto k dimensions using the matrix of canonical vectors V.
}
\usage{
scatterG3(X, Y, V, k = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{X}{
A Nxp data matrix; N observations on the rows and p features on the columns
}
  \item{Y}{
A N vector containing the group labels. Should be coded as 1,2,...,G, where G is the number of groups
}
  \item{V}{
A px(G-1)-dimensional matrix of canonical vectors
}
  \item{k}{Optional number of vectors to use for projection, must be < = G-1. The default value is G-1.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Creates scatteplots of X projected onto k dimensions using matrix V.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Irina Gaynanova
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
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
# visualize the projected data
scatterG3(xtrain,ytrain,V)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
