#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

double Determinant(double **a,int n);

void CoFactor(double **a,int n,double **b);

void Transpose(double **a,int n);

void INVERSE(double **a, int n, double **b);

void MatrixCopy(double **a,double **b,int row,int col,int add);

void XminusY(double **A, double **B, int row, int col, double **target);

void XprimeX(double **A, int row, int col, double **target);

void XprimeY(double **X, double **Y, int Xrow, int Xcol, int Ycol,double **target);

double Trace(double **a,int n);

void Cholesky(double **q,int n,double **lower);

#endif // MATRIX_UTILS_H
