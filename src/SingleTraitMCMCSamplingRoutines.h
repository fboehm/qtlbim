// QTLBIM - QTL Bayesian Interval Mapping
// Functions relevant to single trait analyses
//********************************************************************

#ifndef SINGLE_TRAIT_H
#define SINGLE_TRAIT_H

void ZeroEffect1(int L);
void ZeroEffect2(int L1,int L2);
void Coefficient(int GENOTYPE);

double Likelihood(double *p,double *G);

void Coefficient0(int I,int L,int QL);

double GenotypeSampling(int I,int L,int II,int QL);

void Mean(double YBAR,double VP);

void ResidualVariance();

void MainEffect(int L,int K);
void MainEffect1(int L);

void EpistaticEffect(int L1,int L2,int K1,int K2);
void EpistaticEffect1(int L1,int L2);

void GBYE_FixedCovariate(int L1,int L2,int K);
void GBYE_FixedCovariate1(int L1,int L2);

void MainVariance(int L,int K,int NU,double TAU);
void MainVariance1(int L,int NU,double TAU);

void EpistaticVariance(int L1,int L2,int K1,int K2,int NU,double TAU);
void EpistaticVariance1(int L1,int L2,int NU,double TAU);

void GBYE_FixedCovariate_Variance(int L1,int L2,int K,int NU,double TAU);
void GBYE_FixedCovariate_Variance1(int L1,int L2,int NU,double TAU);

void FixedCovariate(int L);

void RandomCovariate(int L);

void RanVariance(int L);

void QTLgenotype(int L,int NL,int QL,int I);

void QTLPOSITION(int L,int QLNEW);

int SamplingOnePosition(int L);

void MainEffectIndicator_GROUP0(int L,int K);

void MainEffectIndicator_GROUP1(int L);

void EpistasisIndicator_GROUP0(int L1,int L2,int K1,int K2);

void EpistasisIndicator_GROUP1(int L1,int L2);

void GBYE_FIX_Indicator_GROUP0(int L1,int L2,int K);

void GBYE_FIX_Indicator_GROUP1(int L1,int L2);

#endif // SINGLE_TRAIT_H
