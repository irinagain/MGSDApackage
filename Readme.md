## Description of R package MGSDA

MGSDA (Multi-Group Sparse Discriminant Analysis) is an R package that implements methods described in I.Gaynanova, J.Booth and M. Wells "Simultaneous sparse estimation of canonical vectors in the p>>n setting."[http://arxiv.org/abs/1403.6095] 

The package is still under development, but the most recent beta-version can be installed using MGSDA_0.0.tar.gz.

To install the package, call install.packages(“MGSDA_0.0.tar.gz”,type=“source”) from the working directory in which you save the .tar.gz file and then call library(MGSDA). The main functions are cv.dLDA(cross-validation), dLDA(fitting for specified value of tuning parameter) and classifyV(classification). Each function has a documentation  with a simple example which can be acessed using standard ? commands in R (i.e. ?cv.dLDA).

Please feel free to contact me at ig93 [at] cornell [dot] edu if you have any questions/problems with the package.