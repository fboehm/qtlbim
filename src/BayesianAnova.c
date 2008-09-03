#include <stdio.h>
#include <time.h>
#include <math.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>

#include "GlobalVars.h"
#include "GlobalVars_SingleTrait.h"
#include "StatUtils.h"
#include "SingleTraitMCMCSamplingRoutines.h"
#include "BayesianAnova.h"

void bayesianAnova()
{

	int I,J,L,K,K1,K2,L0,L1,L2;

	if(SEED==0)	srand(time(NULL));
	else srand(SEED);
	//double R=rand( );
	rand( );


//*************************************************************************
// calculating phenotypic mean and variance, and assign initial values

double YBAR=0,Y2BAR=0,VP=1;
if(CATEGORY==1)
{
	YBAR=0.0,Y2BAR=0.0;
	for(I=0;I<NS;I++)
	{
		YBAR=YBAR+Y[I];
		Y2BAR=Y2BAR+Y[I]*Y[I];
	}
	YBAR=YBAR/NS;
	VP=(1.0/(NS-1))*(Y2BAR-NS*pow(YBAR,2));


	AMU=YBAR;
	for(I=0;I<NS;I++) VE[I]=VP;            //initial values for AMU and VE
	for(I=0;I<NRANCOVA;I++) VRAN[I]=VP;  //initial values for VRAN[I]
}

if(CATEGORY!=1)
{
	AMU=0;
	for(I=0;I<NS;I++) VE[I]=1.0;              //initial values for AMU and VE
	for(I=0;I<NRANCOVA;I++) VRAN[I]=1.0;  //initial values for VRAN[I]
}

//******************************************************************************
//For specify the prior variances of QTL effects and g by e interactions

double CC[NG];
if(CROSS==2)
{
	CC[0]=1.0/2, CC[1]=1.0/4;   //F2
}
else CC[0]=1.0/4;    //BC

double V_FIX[20];   // 20 = the maximum number of fixed effects
for(L=0;L<NFIXCOVA;L++)
	if(GBYE_FIX_INDEX[L]==1)
	{
		double YBAR1=0.0,YBAR2=0.0;
		for(I=0;I<NS;I++) 
		{
			YBAR1=YBAR1+COEF_FIX[I][L];
			YBAR2=YBAR2+pow(COEF_FIX[I][L],2);
		}
		YBAR1=YBAR1/NS;
		V_FIX[L]=(1.0/(NS-1))*(YBAR2-NS*pow(YBAR1,2));
	}

//******************************************************************************
//Give initial threshold values for ordinal traits

if(CATEGORY==2)
{
	CUTPOINT[0]=-1e+10,CUTPOINT[1]=0,CUTPOINT[2]=1e+10;
}

if(CATEGORY==3)
{
	CUTPOINT[0]=-1e+10, CUTPOINT[1]=0;
	for(J=2;J<CN-1;J++) CUTPOINT[J]=(J-1)*1.0/(CN-2);
	CUTPOINT[CN-1]=1, CUTPOINT[CN]=1e+10;
}
//*************************************************************************
// initialize deviance measures

double AMUBAR,GVALUEBAR[NS],VEBAR;

AMUBAR=0.0;
for (I=0; I<NS; I++)
	GVALUEBAR[I]=0.0;
VEBAR=0.0;

//**********************************************************************************
// give initial values for genotypes, effects, variances

for(L=0;L<NQTL;L++)
{
	for(I=0;I<NS;I++)
	{
		if(UPDATEGENO==1)
		{
			double PRR[NG];
			for(K=0;K<NG;K++) PRR[K]=QPROB[QCHR[L]][I][QLOC[L]][K];
			MULTINORMAL(PRR);
			GENO[I][L]=IBD;
			Coefficient(IBD);
		}
		if(UPDATEGENO==0) Coefficient0(I,L,QLOC[L]);
		for(K=0;K<NC;K++) COEF[I][L][K]=X[K];
	}
}


if(GROUP==0)
{
	for(L=0;L<NQTL;L++)
		for(K=0;K<NC;K++)
		{	
			MainVariance(L,K,1,0.1);
			MainEffect(L,K);
		}

	if(EPISTASIS==1)
	{
		for(L1=0;L1<NQTL-1;L1++)
			for(L2=L1+1;L2<NQTL;L2++) 
				for(K1=0;K1<NC;K1++)
					for(K2=0;K2<NC;K2++) 
					{
						EpistaticVariance(L1,L2,K1,K2,1,0.1);
						EpistaticEffect(L1,L2,K1,K2);
					}
	}

	if(GBYE==1)
	{
		for(L1=0;L1<NFIXCOVA;L1++)
			for(L2=0;L2<NQTL;L2++)
				for(K=0;K<NC;K++) 
				{
					GBYE_FixedCovariate_Variance(L1,L2,K,1,0.1);
					GBYE_FixedCovariate(L1,L2,K);
				}
	}
}


if(GROUP==1)
{
	for(L=0;L<NQTL;L++)
	{
		MainVariance1(L,1,0.1);
		for(K=0;K<NC;K++) MainEffect(L,K);
	}

	if(EPISTASIS==1)
	{
		for(L1=0;L1<NQTL-1;L1++)
			for(L2=L1+1;L2<NQTL;L2++)
			{
				EpistaticVariance1(L1,L2,1,0.1);
				for(K1=0;K1<NC;K1++)
					for(K2=0;K2<NC;K2++) 
						EpistaticEffect(L1,L2,K1,K2);
			}
	}

	if(GBYE==1)
	{
		for(L1=0;L1<NFIXCOVA;L1++)
			for(L2=0;L2<NQTL;L2++)
			{
				GBYE_FixedCovariate_Variance1(L1,L2,1,0.1);
				for(K=0;K<NC;K++) 
					GBYE_FixedCovariate(L1,L2,K);
			}
	}
}

//**********************************************************************************

FILE *File1;
File1=fopen(iterfile,"w");

FILE *File2;
File2=fopen(covfile,"w");

FILE *File3;
File3=fopen(mainfile,"w");

FILE *File4;
File4=fopen(pairfile,"w");

FILE *File5;
File5=fopen(gbyefile,"w");

FILE *File6;
File6=fopen(devfile,"w");

// ***********************************************************
// ITERATION STARTS HERE

int ITER,ITER1;
for(ITER=0;ITER<NITER+(int)(1.0*NBURNIN/NTHIN);ITER++)
{ for(ITER1=0;ITER1<NTHIN;ITER1++)
{

//***********************************************************
//UPDATING THE VALUES OF THE LIABILITY

if(CATEGORY!=1)
	for(I=0;I<NS;I++) Y[I]=TrunNormal(CUTPOINT[W[I]],CUTPOINT[W[I]+1],AMU+GVALUE[I],VE[I]);

if(CATEGORY!=1)
{
	YBAR=0.0,Y2BAR=0.0;
	for(I=0;I<NS;I++)
	{
		YBAR=YBAR+Y[I];
		Y2BAR=Y2BAR+Y[I]*Y[I];
	}
	YBAR=YBAR/NS;
	VP=(1.0/(NS-1))*(Y2BAR-NS*pow(YBAR,2));
}

//*********************************************************************************
//UPDATE PARAMETERS

Mean(YBAR,VP);

if(CATEGORY!=2) ResidualVariance();


if(ENV_FACTOR==1)
{
	for(L=0;L<NRANCOVA;L++)
	{
		RandomCovariate(L);
		RanVariance(L);
	}
	for(L=0;L<NFIXCOVA;L++) FixedCovariate(L);
}


int NU=6; double H=0.1,S=2,TAU; //NU=degrees of freedom, H=heritability, TAU=scale

if(GROUP==0)
{
	for(L=0;L<NQTL;L++)
		if(GAMMA_MAIN[L]!=0.0)
			for(K=0;K<NC;K++)
				if(MAIN[L][K]!=0)
				{
					MainEffect(L,K);
					TAU=(NU-2)*H*VP/(NU*CC[K]);	
					MainVariance(L,K,NU,TAU);
				}

	if(EPISTASIS==1)
	{
		for(L1=0;L1<NQTL-1;L1++)
			for(L2=L1+1;L2<NQTL;L2++)
				if(GAMMA_EPISTASIS[L1][L2]!=0.0) 
					for(K1=0;K1<NC;K1++)
						for(K2=0;K2<NC;K2++) 
							if(EPISTATIC[L1][L2][K1][K2]!=0)
							{
								EpistaticEffect(L1,L2,K1,K2);
								TAU=(NU-2)*H*VP/(NU*CC[K1]*CC[K2]);
								EpistaticVariance(L1,L2,K1,K2,NU,TAU);
							}
	}

	if(GBYE==1)
	{
		for(L1=0;L1<NFIXCOVA;L1++)
			if(GBYE_FIX_INDEX[L1]==1)
				for(L2=0;L2<NQTL;L2++)
					if(GAMMA[L2]!=0.0&&GAMMA_GBYE[L1][L2]!=0.0)
						for(K=0;K<NC;K++) 
							if(GBYE_FIX[L1][L2][K]!=0)
							{
								GBYE_FixedCovariate(L1,L2,K);
								TAU=(NU-2)*H*VP/(NU*V_FIX[L1]*CC[K]);
						        GBYE_FixedCovariate_Variance(L1,L2,K,NU,TAU);
							}
	}
}


if(GROUP==1)
{
	for(L=0;L<NQTL;L++)
		if(GAMMA_MAIN[L]!=0.0) 
		{
			for(K=0;K<NC;K++) MainEffect(L,K);
			TAU=S*VP;
			MainVariance1(L,NU,TAU);
		}

	if(EPISTASIS==1)
	{
		for(L1=0;L1<NQTL-1;L1++)
			for(L2=L1+1;L2<NQTL;L2++)
				if(GAMMA_EPISTASIS[L1][L2]!=0.0) 
				{
					for(K1=0;K1<NC;K1++)
						for(K2=0;K2<NC;K2++) 
							EpistaticEffect(L1,L2,K1,K2);
					TAU=S*VP;
					EpistaticVariance1(L1,L2,NU,TAU);
				}
	}

	if(GBYE==1)
	{
		for(L1=0;L1<NFIXCOVA;L1++)
			if(GBYE_FIX_INDEX[L1]==1)
				for(L2=0;L2<NQTL;L2++)
					if(GAMMA[L2]!=0.0&&GAMMA_GBYE[L1][L2]!=0.0) 
					{
						for(K=0;K<NC;K++) 
							GBYE_FixedCovariate(L1,L2,K);
						TAU=S*VP;
						GBYE_FixedCovariate_Variance1(L1,L2,NU,TAU);
					}
	}
}


//************************************************************************************
//UPDATE THE QTL INHERITANCE OF NON-FOUNDERS AND THE GENOTYPIC VALUES

if(UPDATEGENO==1)
{
	for(L=0;L<NQTL;L++)
	if(GAMMA[L]!=0)
	{
		PD1[L]=0.0,PD2[L]=0.0;
		for(I=0;I<NS;I++)
		{
			int K,KK=0;
			for(K=0;K<NG;K++)
			   if(QPROB[QCHR[L]][I][QLOC[L]][K]>0.99)
			   {
				   IBD=K;
				   KK=1;
			   }
			if(KK==0)
			{
				QTLgenotype(L,QCHR[L],QLOC[L],I);
				PD1[L]=PD1[L]+log(PDD1+1e-20);
				PD2[L]=PD2[L]+log(PDD2+1e-20);
			}

			if(IBD!=GENO[I][L])
			{
				GVALUE[I]=GenotypeSampling(I,L,IBD,0);
				GENO[I][L]=IBD;
				Coefficient(IBD);
				for(K=0;K<NC;K++) COEF[I][L][K]=X[K];
			}
		}
	}
}

//**********************************************************************************
//UPDATING QTL POSITIONS

if(UPDATEPOS==1)
{
	for(L=0;L<NQTL;L++)
	if(GAMMA[L]!=0)
	{
		int CAT[4],QLNEW,TEST=0;

		double R=RANDOM();
		CAT[0]=(R>=0&&R<0.25),CAT[1]=(R>=0.25&&R<0.5),CAT[2]=(R>=0.5&&R<0.75),CAT[3]=(R>=0.75&&R<=1.0);
		QLNEW=CAT[0]*(QLOC[L]-2)+CAT[1]*(QLOC[L]-1)+CAT[2]*(QLOC[L]+1)+CAT[3]*(QLOC[L]+2);
		if(QLNEW<0) QLNEW=0;
		if(QLNEW>=NGRID[QCHR[L]]-1) QLNEW=NGRID[QCHR[L]]-1;

		if(QLNEW!=QLOC[L])
		{
			for(L0=0;L0<NQTL;L0++)
				if( QCHR[L0]==QCHR[L]&&L0!=L )
					TEST=TEST+(fabs(GRID[QCHR[L]][QLNEW]-GRID[QCHR[L0]][QLOC[L0]])<=DQQ[QCHR[L]]);

			if(TEST==0) QTLPOSITION(L,QLNEW);
		}
	}
}

//**********************************************************
//update the threshold values

if(CATEGORY==3)
{
	int j; double CUTPOINT0[CN+1],AL0,AL;

	for(j=0;j<=CN;j++) CUTPOINT0[j]=CUTPOINT[j];

	AL0=Likelihood(CUTPOINT,GVALUE);

	for(j=2;j<=CN-2;j++)
	{
		CUTPOINT0[j]=CUTPOINT[j]+0.01*(RANDOM()-0.5);
		if(CUTPOINT0[j]<=CUTPOINT0[j-1]) CUTPOINT0[j]=CUTPOINT0[j-1]+0.01*RANDOM();
		if(CUTPOINT0[j]>CUTPOINT[j+1]) CUTPOINT0[j]=CUTPOINT[j+1]-0.01*RANDOM();
	}

	AL=Likelihood(CUTPOINT0,GVALUE);

	if((AL-AL0)>log(RANDOM()))
		for(j=2;j<=CN-2;j++) CUTPOINT[j]=CUTPOINT0[j];
}

//***************************************************************

}    //ITER1 end here

//***************************************************************

//***************************************************************
//SAVE THE RESULT

if(ITER*ITER1>=NBURNIN)
{

if((ITER!=0)&&(ITER%200==0)&(VERBOSE>0))
{
	Rprintf("%d",ITER);
	Rprintf("\n");
}

if((ITER!=0)&&((ITER%200)==0))
   R_CheckUserInterrupt();

//CALCULATE THE NUMBER OF QTL
int QTL_INCLUDED=0;
for(L=0;L<NQTL;L++) QTL_INCLUDED=QTL_INCLUDED+(GAMMA[L]!=0);

double VE0=0;
for(I=0;I<NS;I++) VE0=VE0+VE[I];
VE0=VE0/NS;

//CALCULATE DEVIANCE
AMUBAR+=AMU;
VEBAR+=VE0;
for (I=0; I<NS; I++) GVALUEBAR[I]+=GVALUE[I];

double DEV=0.0;
for(I=0;I<NS;I++) DEV=DEV+log(2*3.14159265358*VE[I])+pow(Y[I]-AMU-GVALUE[I],2)/VE[I];

//save to "deviance"
fprintf(File6,"%lf\n",DEV);

// save to "iterdiag"
fprintf(File1,"\n");
fprintf(File1,"%d\t",ITER+1);
fprintf(File1,"%d\t",QTL_INCLUDED);
fprintf(File1,"%f\t",AMU);
fprintf(File1,"%f\t",VE0);

// save to "covariates"
if(ENV_FACTOR==1)
{
	fprintf(File2,"\n");
	for(L=0;L<NFIXCOVA;L++) fprintf(File2,"%f\t",FIX[L]);
	for(L=0;L<NRANCOVA;L++) fprintf(File2,"%f\t",VRAN[L]);
}
if(CATEGORY==3)
  for(L=2;L<=CN-2;L++) fprintf(File2,"%f\t",CUTPOINT[L]);



if(GROUP==0)
{

// save to "mainloci"
double VAR1[NG];
for(K=0;K<NC;K++) VAR1[K]=0;
for(L=0;L<NQTL;L++)
	if(GAMMA[L]!=0)
	{
		fprintf(File3,"\n");
		fprintf(File3,"%d\t",ITER+1);
        fprintf(File3,"%d\t",QTL_INCLUDED);
        fprintf(File3,"%d\t",QCHR[L]+1);
//        fprintf(File3,"%f\t",GRID[QCHR[L]][QLOC[L]]);
		fprintf(File3,"%d\t",QLOC[L]);
		for(K=0;K<NC;K++) fprintf(File3,"%f\t",MAIN[L][K]);
		for(K=0;K<NC;K++)
		{
			double SS1=0,SS2=0,SS=0;
			if(MAIN[L][K]!=0)
			{
				for(I=0;I<NS;I++)
				{
					SS1=SS1+pow(COEF[I][L][K],2);
					SS2=SS2+COEF[I][L][K];
				}
				SS=1.0/(NS-1)*(SS1-NS*pow(SS2/NS,2))*pow(MAIN[L][K],2);
			}
			fprintf(File3,"%f\t",SS);
			VAR1[K]=VAR1[K]+SS;
		}
	}
for(K=0;K<NC;K++) fprintf(File1,"%f\t",VAR1[K]);

// save to "pairloci"
double VAR2[NG][NG];
for(K1=0;K1<NC;K1++)
	for(K2=0;K2<NC;K2++) VAR2[K1][K2]=0;
if(EPISTASIS==1)
{
	int N_EPIS=0;
	for(L1=0;L1<NQTL;L1++)
		for(L2=L1+1;L2<NQTL;L2++)
			for(K1=0;K1<NC;K1++)
				for(K2=0;K2<NC;K2++)
					if(EPISTATIC[L1][L2][K1][K2]!=0) N_EPIS=N_EPIS+1;

	for(L1=0;L1<NQTL;L1++)
		for(L2=L1+1;L2<NQTL;L2++)
			if(GAMMA[L1]!=0&&GAMMA[L2]!=0&&GAMMA_EPISTASIS[L1][L2]!=0)
			{
				fprintf(File4,"\n");
				fprintf(File4,"%d\t",ITER+1);
				fprintf(File4,"%d\t",N_EPIS);
				fprintf(File4,"%d\t",QCHR[L1]+1);
//				fprintf(File4,"%f\t",GRID[QCHR[L1]][QLOC[L1]]);
				fprintf(File4,"%d\t",QLOC[L1]);
				fprintf(File4,"%d\t",QCHR[L2]+1);
//				fprintf(File4,"%f\t",GRID[QCHR[L2]][QLOC[L2]]);
				fprintf(File4,"%d\t",QLOC[L2]);
				for(K1=0;K1<NC;K1++)
					for(K2=0;K2<NC;K2++) fprintf(File4,"%f\t",EPISTATIC[L1][L2][K1][K2]);

				for(K1=0;K1<NC;K1++)
				for(K2=0;K2<NC;K2++)
				{
					double SS1=0,SS2=0,SS=0;
					if(EPISTATIC[L1][L2][K1][K2]!=0)
					{
						for(I=0;I<NS;I++)
						{
							SS1=SS1+pow(COEF[I][L1][K1]*COEF[I][L2][K2],2);
							SS2=SS2+COEF[I][L1][K1]*COEF[I][L2][K2];
						}
						SS=1.0/(NS-1)*(SS1-NS*pow(SS2/NS,2))*pow(EPISTATIC[L1][L2][K1][K2],2);
					}
					fprintf(File4,"%f\t",SS);
					VAR2[K1][K2]=VAR2[K1][K2]+SS;
				}
			}
    for(K1=0;K1<NC;K1++)
		for(K2=0;K2<NC;K2++) fprintf(File1,"%f\t",VAR2[K1][K2]);
}


// save to "gbye"
double VAR3[NG];
for(K=0;K<NC;K++) VAR3[K]=0;
if(GBYE==1)
{
	int N_GBYE=0;
	for(L1=0;L1<NFIXCOVA;L1++)
		for(L2=0;L2<NQTL;L2++)
			for(K=0;K<NC;K++)
				if(GBYE_FIX[L1][L2][K]!=0) N_GBYE=N_GBYE+1;

	for(L1=0;L1<NFIXCOVA;L1++)
		 if(GBYE_FIX_INDEX[L1]==1)
			for(L2=0;L2<NQTL;L2++)
				if(GAMMA[L2]!=0.0)
					if(GAMMA_GBYE[L1][L2]!=0.0)
					{
						fprintf(File5,"\n");
						fprintf(File5,"%d\t",ITER+1);
						fprintf(File5,"%d\t",N_GBYE);
						fprintf(File5,"%d\t",L1+1);
						fprintf(File5,"%d\t",QCHR[L2]+1);
					//	fprintf(File5,"%f\t",GRID[QCHR[L2]][QLOC[L2]]);
						fprintf(File5,"%d\t",QLOC[L2]);
						for(K=0;K<NC;K++) fprintf(File5,"%f\t",GBYE_FIX[L1][L2][K]);

						for(K=0;K<NC;K++)
						{
							double SS1=0,SS2=0,SS=0;
							if(GBYE_FIX[L1][L2][K]!=0)
							{
								for(I=0;I<NS;I++)
								{
									SS1=SS1+pow(COEF_FIX[I][L1]*COEF[I][L2][K],2);
									SS2=SS2+COEF_FIX[I][L1]*COEF[I][L2][K];
								}
								SS=1.0/(NS-1)*(SS1-NS*pow(SS2/NS,2))*pow(GBYE_FIX[L1][L2][K],2);
							}
							fprintf(File5,"%f\t",SS);
							VAR3[K]=VAR3[K]+SS;
						}
					}
	for(K=0;K<NC;K++) fprintf(File1,"%f\t",VAR3[K]);
}

double VAR=0;
for(K1=0;K1<NC;K1++)
{
	VAR=VAR+VAR1[K1]+VAR3[K1];
	for(K2=0;K2<NC;K2++) VAR=VAR+VAR2[K1][K2];
}
if(ENV_FACTOR==1) fprintf(File1,"%f\t",VP-VAR-VE0);
fprintf(File1,"%f\t",VAR);

}      //GROUP=0 end


if(GROUP==1)
{

// save to "mainloci"
double VAR1=0;
for(L=0;L<NQTL;L++)
if(GAMMA[L]!=0)
{
	fprintf(File3,"\n");
	fprintf(File3,"%d\t",ITER+1);
    fprintf(File3,"%d\t",QTL_INCLUDED);
    fprintf(File3,"%d\t",QCHR[L]+1);
//    fprintf(File3,"%f\t",GRID[QCHR[L]][QLOC[L]]);
	fprintf(File3,"%d\t",QLOC[L]);
	for(K=0;K<NC;K++) fprintf(File3,"%f\t",MAIN[L][K]);

	double SS1=0,SS2=0,SS=0;
	for(I=0;I<NS;I++)
	{
		double S=0;
		for(K=0;K<NC;K++) S=S+COEF[I][L][K]*MAIN[L][K];
		SS1=SS1+pow(S,2);
		SS2=SS2+S;
	}
	SS=1.0/(NS-1)*(SS1-NS*pow(SS2/NS,2)); //main-effect variance
	fprintf(File3,"%f\t",SS);
	VAR1=VAR1+SS;
}
fprintf(File1,"%f\t",VAR1);


// save to "pairloci"
double VAR2=0;
if(EPISTASIS==1)
{
	int N_EPIS=0;
	for(L1=0;L1<NQTL;L1++)
		for(L2=L1+1;L2<NQTL;L2++)
			if(GAMMA_EPISTASIS[L1][L2]!=0) N_EPIS=N_EPIS+1;

	for(L1=0;L1<NQTL;L1++)
		for(L2=L1+1;L2<NQTL;L2++)
		if(GAMMA[L1]!=0&&GAMMA[L2]!=0&&GAMMA_EPISTASIS[L1][L2]!=0)
		{
			fprintf(File4,"\n");
			fprintf(File4,"%d\t",ITER+1);
			fprintf(File4,"%d\t",N_EPIS);
			fprintf(File4,"%d\t",QCHR[L1]+1);
//			fprintf(File4,"%f\t",GRID[QCHR[L1]][QLOC[L1]]);
			fprintf(File4,"%d\t",QLOC[L1]);
			fprintf(File4,"%d\t",QCHR[L2]+1);
//			fprintf(File4,"%f\t",GRID[QCHR[L2]][QLOC[L2]]);
			fprintf(File4,"%d\t",QLOC[L2]);
			for(K1=0;K1<NC;K1++)
				for(K2=0;K2<NC;K2++) fprintf(File4,"%f\t",EPISTATIC[L1][L2][K1][K2]);

			double SS1=0,SS2=0,SS=0;
			for(I=0;I<NS;I++)
			{
				double S=0;
				for(K1=0;K1<NC;K1++)
					for(K2=0;K2<NC;K2++)
						S=S+COEF[I][L1][K1]*COEF[I][L2][K2]*EPISTATIC[L1][L2][K1][K2];
				SS1=SS1+pow(S,2);
				SS2=SS2+S;
			}
			SS=1.0/(NS-1)*(SS1-NS*pow(SS2/NS,2)); //epistatic variance
			fprintf(File4,"%f\t",SS);
			VAR2=VAR2+SS;
		}
		fprintf(File1,"%f\t",VAR2);
}


// save to "gbye"
double VAR3=0;
if(GBYE==1)
{
	int N_GBYE=0;
	for(L1=0;L1<NFIXCOVA;L1++)
		for(L2=0;L2<NQTL;L2++)
			if(GAMMA_GBYE[L1][L2]!=0.0) N_GBYE=N_GBYE+1;

	for(L1=0;L1<NFIXCOVA;L1++)
		 if(GBYE_FIX_INDEX[L1]==1)
			for(L2=0;L2<NQTL;L2++)
			if(GAMMA[L2]!=0.0)
			if(GAMMA_GBYE[L1][L2]!=0.0)
			{
				fprintf(File5,"\n");
				fprintf(File5,"%d\t",ITER+1);
				fprintf(File5,"%d\t",N_GBYE);
				fprintf(File5,"%d\t",L1+1);
				fprintf(File5,"%d\t",QCHR[L2]+1);
//				fprintf(File5,"%f\t",GRID[QCHR[L2]][QLOC[L2]]);
				fprintf(File5,"%d\t",QLOC[L2]);
				for(K=0;K<NC;K++) fprintf(File5,"%f\t",GBYE_FIX[L1][L2][K]);

				double SS1=0,SS2=0,SS=0;
				for(I=0;I<NS;I++)
				{
					double S=0;
					for(K=0;K<NC;K++) S=S+COEF_FIX[I][L1]*COEF[I][L2][K]*GBYE_FIX[L1][L2][K];
					SS1=SS1+pow(S,2);
					SS2=SS2+S;
				}
				SS=1.0/(NS-1)*(SS1-NS*pow(SS2/NS,2)); //GXE variance
				fprintf(File5,"%f\t",SS);
				VAR3=VAR3+SS;
			}
			fprintf(File1,"%f\t",VAR3);
}

double VAR=VAR1+VAR2+VAR3;  //total genetic variance
if(ENV_FACTOR==1) fprintf(File1,"%f\t",VP-VAR-VE0);   //covariate variance
fprintf(File1,"%f\t",VAR);

} //GROUP=1 end


}

}

//ITER end here

//CALCULATE DHAT
double DHAT=0.0;
AMUBAR/=(double)(NITER);
VEBAR/=(double)(NITER);
for (I=0; I<NS; I++)
{
	GVALUEBAR[I]/=(double)(NITER);
	DHAT=DHAT+log(2*3.14159265358*VEBAR)+pow(Y[I]-AMUBAR-GVALUEBAR[I],2)/VEBAR;
}
//save DHAT
fprintf(File6,"%lf\n",DHAT);

fclose(File1);
fclose(File2);
fclose(File3);
fclose(File4);
fclose(File5);
fclose(File6);
}
