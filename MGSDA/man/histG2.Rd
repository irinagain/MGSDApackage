\name{histG2}
\alias{histG2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Histogram of projected data for the two-group case
}
\description{
Creates a histogram of the data projected onto one dimension using canonical vector V.
}
\usage{
histG2(X, Y, V)
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
A p-dimensional canonical vector
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plots a histogram of resulting scores for two groups.
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
G <- 2
ytrain <- rep(1:G,each=n)
set.seed(1)
xtrain <- matrix(rnorm(p*n*G),n*G,p)
# find matrix of canonical vectors V
V <- dLDA(xtrain,ytrain,lambda=0.1)
sum(rowSums(V)!=0)
# visualize the data
histG2(xtrain,ytrain,V)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
