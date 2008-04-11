#include <math.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>

#include "GlobalVars.h"
#include "StatUtils.h"

double RANDOM()
{
	double U=rand()/(RAND_MAX+1.0);
	if(U>=1.0) U=1.0-(1e-10); if(U<=0.0) U=1e-10;
	return(U);
}

//*******************************************************
void ANORMAL(double  *R1, double  *R2)
{
	double U1=RANDOM();
    double U2=RANDOM();
    *R1=sqrt(-2.*log(U1))*cos(6.2830*U2);
    *R2=sqrt(-2.*log(U1))*sin(6.2830*U2);
	return;
}
//********************************************************
void MULTINORMAL(double * PRIOR)
{
	double SUM[NG];int I,J;

    for(I=0;I<NG;I++)
	{
        SUM[I]=0.0;
        for(J=0;J<=I;J++) SUM[I]=SUM[I]+PRIOR[J];
	}

    double R=RANDOM();
    if(R<=SUM[0]) IBD=0;
    else
		for(I=0;I<NG-1;I++)
			if((R>SUM[I])&&(R<=SUM[I+1])) IBD=I+1;
    return;
}

//****************************************************************************************
///Calculate the normal function value
double NormalFunction(double X)
{
	double T,Z,ANS;

	Z=fabs(X)/sqrt(2.0);
	T=1.0/(1.0+0.5*Z);
	ANS=T*exp(-Z*Z-1.26551223+T*(1.00002368+T*(0.37409196+T*(0.09678418+
		T*(-0.18628806+T*(0.27886807+T*(-1.13520398+T*(1.48851587+
		T*(-0.82215223+T*0.17087277)))))))));
	return X>=0.0 ? 0.5*(2-ANS) : 1-0.5*(2-ANS);
}

//*****************************************************************************************
///Sample from a double-truncated normal density
double TrunNormal(double T1,double T2,double B,double V)
{
	double u,p1,p2,p,alpha=0,y=0,yy=0;
	double c0=2.515517,c1=0.802853,c2=0.010328,d1=1.432788,d2=0.189269,d3=0.001308;
//	double T1=CUTPOINT[WW],T2=CUTPOINT[WW+1];

	u=rand()/(RAND_MAX+1.0);

	p1=NormalFunction((T1-B)/sqrt(V));
	p2=NormalFunction((T2-B)/sqrt(V));

	p=p1+u*(p2-p1);
	if(p==0.0) p=1e-10;
	if(p==1.0) p=1-(1e-10);

	if(p>0.0&&p<0.5)
	{
		alpha=p;
		y=sqrt(-2*log(alpha));
	    yy=y-(c0+c1*y+c2*y*y)/(d1*y+d2*y*y+d3*y*y*y+1);
		yy=-yy;
	}
	if(p>0.5&&p<1.0)
	{
		alpha=1-p;
		y=sqrt(-2*log(alpha));
	    yy=y-(c0+c1*y+c2*y*y)/(d1*y+d2*y*y+d3*y*y*y+1);
	}
	if(p==0.5) yy=0;

	return(B+sqrt(V)*yy);
}
