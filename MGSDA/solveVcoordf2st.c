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

////////////////////////
// old MGSDA functions
/////////////////////////
//here B is matrix and we want to multiply each a[i] by B[i,l]
void mulV(double *a, double *B, int* length, int *l, double* runningSum){
  int index;
  for (index=0; index<*length; index++)
    { 
      *runningSum=*runningSum+a[index]*B[index+(*l)*(*length)];
   }
}

//standardized version, W[k,k]=1
void blockUpdate2(double *D, double *W,double *V, double *lambda, int *p, int *r, double *errV)
{
  int k,l;
  double normt;
  double colsum;
  double tmp[*r]; //this is a vector of length r
  double v_old;

  *errV=0;
  for (k=0; k<*p; k++)
    {
      //k=cursample[i]-1; //indexing in C starts from 0 so when I refer to k, need to adjust its number by 1
      
      for (l=0;l<*r;l++){
          colsum=0;
          mulV(&W[k*(*p)],V,p,&l,&colsum);//this multiplies W2[k,] by V[,l]
          tmp[l]=D[k+l*(*p)]-colsum+V[k+l*(*p)];
      }

      normt=0;
      mul(&tmp[0],&tmp[0],r,&normt); //this makes normt=sum(tmp^2)
      normt=sqrt(normt);

      if (normt<= *lambda){
          for (l=0;l<*r;l++){
              v_old=V[k+l*(*p)];
              V[k+l*(*p)]=0; //V[i,j]=V[j*(*Nrow)+i]
              *errV=max(*errV,fabs(v_old-V[k+l*(*p)]));
          }
      }
      else{
          for (l=0;l<*r;l++){
              v_old=V[k+l*(*p)];
              V[k+l*(*p)]=(1-*lambda/normt)*tmp[l];
              *errV=max(*errV,fabs(v_old-V[k+l*(*p)]));
          }
      }
    }
}

void coordUpdate(double *D, double *W,double *V, double *lambda, int *p, double *errV)
{
  int k;
  double colsum;
  double v_old;

  *errV=0;
  for (k=0; k<*p; k++)
    {
      //k=cursample[i]-1; //indexing in C starts from 0 so when I refer to k, need to adjust its number by 1

      colsum=0;
      mul(&W[k*(*p)],V,p,&colsum);//this multiplies W2[k,] by V[,l]
      colsum=D[k]-colsum+V[k];
      
      v_old=V[k];
      V[k]=((colsum>0)-(colsum<0))*max(0,fabs(colsum)-*lambda);
      *errV=max(*errV,fabs(v_old-V[k]));
    }
}

void solveVcoordf2(double *D, double *W,double *V, double *lambda, int *p, int *r,double *eps, int *maxiter,int *niter){

  // int i;
  double errV;
  //int cursample[*p];

  *niter=0;
  if(*r>1){
    do{
      //for now do linear update each time
      // this gives almost the same time as random update in R
      //for (i=0;i<*p;i++){
      //  cursample[i]=i+1;
      // }
      ++*niter;
      blockUpdate2(D,W,V,lambda,p,r,&errV);
    }while ((errV>*eps)&&(*niter<*maxiter));
  }
  else{
    do{
      //for now do linear update each time
      // this gives almost the same time as random update in R
      //for (i=0;i<*p;i++){
      //  cursample[i]=i+1;
      // }
      ++*niter;
      coordUpdate(D,W,V,lambda,p,&errV);
    }while ((errV>*eps)&&(*niter<*maxiter));
  }
}

