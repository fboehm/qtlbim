#include "MatrixUtils.h"

#include <R.h>

#include <math.h>
#include <stdlib.h>


///  Recursive definition of determinate using expansion by minors.
double Determinant(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

/*         m = malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = malloc((n-1)*sizeof(double)); */

   if (n < 1) { 
   Rprintf("\n Major error in C program's determinant function");
   return(det); /* exit(1); not allowed here */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = malloc((n-1)*sizeof(double));
         for(i=0;i<n-1;i++)
            

         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det = det + pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}


///  Find the cofactor matrix of a square matrix copied to b.
void CoFactor(double **a,int n,double **b)
{
   int i,j,ii,jj,i1,j1;
   double det;
   double **c;

   c = malloc((n-1)*sizeof(double *));
   for (i=0;i<n-1;i++)
     c[i] = malloc((n-1)*sizeof(double));

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = Determinant(c,n-1);

         /* Fill in the elements of the cofactor */
         b[i][j] = pow(-1.0,i+j+2.0) * det;
      }
   }
   for (i=0;i<n-1;i++)
      free(c[i]);
   free(c);
}


///  Transpose of a square matrix, do it in place
void Transpose(double **a,int n)
{
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = tmp;
      }
   }
}


/// Inverse of a matrix copied to b
void INVERSE(double **a, int n, double **b)
{
int i,j;
double det;
 det=Determinant(a,n);
 CoFactor(a,n,b);
 Transpose(b,n);
 for(i=0;i<n;i++)
	 for(j=0;j<n;j++)
		 b[i][j]=b[i][j]/det;

 return;
}


/// Copy Matrix "a" to "b" and store it there with an option to add
void MatrixCopy(double **a,double **b,int row,int col,int add)
{
int i,j;
if(add==0)
{
 for(i=0;i<row;i++)
    for(j=0;j<col;j++)
      b[i][j]=a[i][j];
   }
 else if(add==1)
  { 
    for(i=0;i<row;i++)
    for(j=0;j<col;j++)
      b[i][j]=b[i][j]+a[i][j];
   }     
 }


/// Matrix substraction A-B 
void XminusY(double **A, double **B, int row, int col, double **target)
{
int i,j;
for(i=0;i<row;i++)
  for(j=0;j<col;j++)
    target[i][j]=A[i][j]-B[i][j];
}


/// Multiply A'A where A is m x n and A'A = n x n (for m > n) for m < n AA'
void XprimeX(double **A, int row, int col, double **target)
{
 int i,j,j1;
 double sum=0;
if(row >= col)
{
 for(j=0;j<col;j++)
  {
   for(j1=0;j1<col;j1++)
    {
    sum=0;
    for(i=0;i<row;i++)          sum = sum + A[i][j]*A[i][j1];          
     target[j][j1]=sum;   
    }
  }
 }
 else{
 i=row;
 row=col;
 col=i;
 for(j=0;j<col;j++)
    {
   for(j1=0;j1<col;j1++)
      {
      sum=0;
    for(i=0;i<row;i++)          sum = sum + A[j][i]*A[j1][i];          
     target[j][j1]=sum;   
      }
    }
  }
 
 } 


/// Multiply AB where A is m x n and B = n x p target has to be of the dimension Xrow,Ycol
void XprimeY(double **X, double **Y, int Xrow, int Xcol, int Ycol,double **target)
{
 int i,j,j1;
 double sum=0;
 for(j=0;j<Xrow;j++)
  {
   for(j1=0;j1<Ycol;j1++)
    {
    sum=0;
    for(i=0;i<Xcol;i++)          sum = sum + X[j][i]*Y[j1][i];          
     target[j][j1]=sum;   
    }
  }
 
 } 


/// Computes the trace of a square matrix
double Trace(double **a,int n)
{
 double trace=0;
 int i;
 for(i=0;i<n;i++)
       trace = trace + a[i][i];
 return(trace);      
}


/// Computes the Cholesky factor for a positive definite symm square matrix
void Cholesky(double **q,int n,double **lower)
{
double linsum=0,crosssum=0;
int i,j,k;
  for(i=0;i<n;i++)
      for(j=0;j<n;j++)
            lower[i][j]=0;               
  for(i=0;i<n;i++)
   {  linsum=0;
      for(j=0;j<=i;j++)
      {
       if(j==i) lower[i][j] = sqrt(q[i][j]-linsum);
       else {
        for(k=0;k<j;k++) crosssum = crosssum + lower[i][k]*lower[j][k];
        lower[i][j] = (q[j][i] - crosssum)/lower[j][j];
        }
       linsum = linsum + lower[i][j]*lower[i][j];
       crosssum=0;
       }
      }
}   
