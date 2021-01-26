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

//function for weighted residual
void resW(double *X, double *Y, double *beta, int *n, int *p, double *residual, double *W){
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
void coordUpdateLasso(double *X,double *beta, double *lambda, int *n, int *p, double *errb,double *residual, double *normx)
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
      coordUpdateLasso(X,beta,lambda,n,p,&errb,&residual[0],&normx[0]);
    }while ((errb>*eps)&&(*niter<*maxiter));
}

//Function for one round of block update, here V is matrix
void blockUpdateLasso(double *X,double *V, double *lambda, int *n, int *p, int *r, double *errV,double *residual, double *normx)
{
  int k,j,l;
  double normt,v_old;
  double T[*r];

  *errV=0;
  for (k=0; k<*p; k++)
    {
      //calculate T
      for (l=0;l<*r;l++){ 
        T[l]=mul(&X[k*(*n)],&residual[l*(*n)],n)/(*n)+normx[k]*V[k+l*(*p)];// be careful with residual indexing
      }
      
      //calculate L_2 norm of T
      normt=mul(&T[0],&T[0],r); //this makes normt=sum(T^2)
      normt=sqrt(normt);
      
      //update kth row of V
      if (normt<= *lambda){
          for (l=0;l<*r;l++){
              v_old=V[k+l*(*p)];
              V[k+l*(*p)]=0; //V[i,j]=V[j*(*Nrow)+i]
              //update residual
              for (j=0; j<*n;j++){
                residual[j+l*(*n)]+=(v_old-V[k+l*(*p)])*X[j+k*(*n)];
              }   
              *errV=max(*errV,fabs(v_old));
          }
      }
      else{
          for (l=0;l<*r;l++){
              v_old=V[k+l*(*p)];
              V[k+l*(*p)]=(1-*lambda/normt)*T[l]/normx[k];
              //update residual
              for (j=0; j<*n;j++){
                residual[j+l*(*n)]+=(v_old-V[k+l*(*p)])*X[j+k*(*n)];
              }  
              *errV=max(*errV,fabs(v_old-V[k+l*(*p)]));
          }
      }
    }
}

//function for residual that is a matrux
void resF(double *X, double *Y, double *beta, int *n, int *p, int *r, double *residual){
    int j,k,l;
    for (j=0; j<*n;j++){
        for (l=0; l<*r; ++l){
            residual[j+l*(*n)]=Y[j+l*(*n)];
            for (k=0; k<*p; ++k){
                residual[j+l*(*n)]-=X[j+k*(*n)]*beta[k+l*(*p)];
            }
        }
    }
}

//Complete lasso solver with Frobenius norm
void solveMyLassoF(double *X, double *Y,double *beta, double *lambda, int *p, int *n, int *r, double *eps, int *maxiter,int *niter){

  double errV;
  double residual[(*n)*(*r)]; //this is a n times r matrix
  double normx[*p]; //this is a vector of length p

    //calculate the current residual between y and Xbeta
    resF(X,Y,beta,n,p,r,&residual[0]);
    //calculate l2 norm of each column of X standardized by n
    colNorm(X,n,p,&normx[0]);

  *niter=0;
    do{
      ++*niter;
      blockUpdateLasso(X,beta,lambda,n,p,r,&errV,&residual[0],&normx[0]);
    }while ((errV>*eps)&&(*niter<*maxiter));
}

//Complete lasso solver for weighted lasso, i.e. X'WX, and X'Y, with Frobenius norm
void solveMyLassoFW(double *X, double *Y,double *W, double *beta, double *lambda, int *p, int *n, int *r, double *eps, int *maxiter,int *niter){
    
    double errV;
    double residual[(*n)*(*r)]; //this is a n times r matrix
    double normx[*p]; //this is a vector of length p
    
    //
    
    //calculate the current residual between y and Xbeta
    resF(X,Y,beta,n,p,r,&residual[0]);
    //calculate l2 norm of each column of X standardized by n
    colNorm(X,n,p,&normx[0]);
    
    *niter=0;
    do{
        ++*niter;
        blockUpdateLasso(X,beta,lambda,n,p,r,&errV,&residual[0],&normx[0]);
    }while ((errV>*eps)&&(*niter<*maxiter));
}