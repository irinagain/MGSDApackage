## Description of R package MGSDA

MGSDA (Multi-Group Sparse Discriminant Analysis) is an R package that implements methods described in 

  * [I.Gaynanova, J.Booth and M. Wells *"Simultaneous sparse estimation of canonical vectors in the p>>n setting."*, JASA, 111(514), 696-706.](https://dx.doi.org/10.1080/01621459.2015.1034318)

The package is available from CRAN.

To install from Github:

``` install
devtools::install_github("irinagain/MGSDApackage")
```

The main functions are `cv.dLDA`(cross-validation), `dLDA`(fitting for specified value of tuning parameter) and `classifyV`(classification). Each function has a documentation  with a simple example which can be accessed using standard ? commands in R (i.e. `?cv.dLDA`).

Please feel free to contact me at irinag [at] stat [dot] tamu [dot] edu if you have any questions or experience problems with the package.

### Simple example

```r
library(MGSDA)

### Example 1
# generate training data
n <- 10
p <- 100
G <- 3
ytrain <- rep(1:G, each = n)
set.seed(1)
xtrain <- matrix(rnorm(p * n * G), n * G, p)

# find matrix of canonical vectors V
V <- dLDA(xtrain, ytrain, lambda = 0.1)
sum(rowSums(V) != 0)

# generate test data
m <- 20
set.seed(3)
xtest <- matrix(rnorm(p * m), m, p)

# perform classification
ytest <- classifyV(xtrain, ytrain, xtest, V)
```