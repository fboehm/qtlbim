#ifndef STAT_UTILS_H
#define STAT_UTILS_H

double RANDOM();
void ANORMAL(double *R1,double *R2);
void MULTINORMAL(double * PRIOR);
double NormalFunction(double X);
double TrunNormal(double T1,double T2,double B,double V);
void Multivariate_RNORM(double **sigma_cholesky,double *mu,int m,double *sample);

#endif // STAT_UTILS_H
