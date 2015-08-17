#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// new functions in C for MGSDA package

//function for maximum
double max(double a, double b){
    if (a>=b){
        return a;
    }else{
        return b;
    }
}

//function for inner product of two vectors
double mul(double *a, double *b, int* length){
    int index;
    double result=0;
    
    for (index = 0; index < *length; index++){
        result+=a[index] * b[index];
    }
    return result;
}

//////////////////////////////////////////
// Set of functions for LASSO problem  ///
//////////////////////////////////////////
//function for residual
void res(double *X, double *Y, double *beta, int *n, int *p, double *residual){
    int j,k;
    for (j=0; j<*n;j++){
        residual[j]=Y[j];
        for (k=0; k<*p; ++k){
            residual[j]-=X[j+k*(*n)]*beta[k];
        }
    }
}

//function for calculating normx=colSums(x^2)/n
void colNorm(double *X, int *n, int* p, double* normx){
  int index;
  for (index=0; index<*p; index++){
        normx[index]=mul(&X[index*(*n)],&X[index*(*n)],n)/(*n);
   }
}

//Function for one round of coordinate update
void coordUpdateLasso(double *X, double *Y,double *beta, double *lambda, int *n, int *p, double *errb,double *residual, double *normx)
{
  int k,j;
  double colsum;
  double bold;

  *errb=0;
  for (k=0; k<*p; k++)
    {
      bold=beta[k];
      colsum=mul(&X[k*(*n)],residual,n)/(*n)+normx[k]*bold;
      beta[k]=((colsum>0)-(colsum<0))*max(0,fabs(colsum)-*lambda)/normx[k];
        //update residual
        for (j=0; j<*n;j++){
            residual[j]+=(bold-beta[k])*X[j+k*(*n)];
        }
      *errb=max(*errb,fabs(bold-beta[k]));
    }
}

//Complete lasso solver
void solveMyLasso(double *X, double *Y,double *beta, double *lambda, int *p, int *n, double *eps, int *maxiter,int *niter){

  double errb;
  double residual[*n]; //this is a vector of length n
  double normx[*p]; //this is a vector of length p

    //calculate the current residual between y and Xbeta
    res(X,Y,beta,n,p,&residual[0]);
    //calculate l2 norm of each column of X standardized by n
    colNorm(X,n,p,&normx[0]);

  *niter=0;
    do{
      ++*niter;
      coordUpdateLasso(X,Y,beta,lambda,n,p,&errb,&residual[0],&normx[0]);
    }while ((errb>*eps)&&(*niter<*maxiter));
}
