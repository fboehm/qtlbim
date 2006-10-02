

// BAYESIAN MAPPING OF QUANTITATIVE TRAIT LOCI  
//********************************************************************

#include <stdio.h>
#include <time.h>
#include <math.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>

#include "bmq_mcmc.h"

// Defining global variables

int CROSS;         // cross type 0:RILs, 1:BC, 2:F2
int NG;            // number of genotypes
int NS;            // number of individuals 
int NS1;           // number of individuals 
int NLG;           // number of chromosomes
int M;             // number of marker   

//********************************************************************************

int HALDANE;      // 1: Haldane distance; 0: Kosambi distance
double CSTEP;     // length of the grid
double *RECM;     // marker map
double LCHR;      // chr length
int **MK;         // marker genotype
int KCHR;         // number of markers

double ***GenoProb;
double *CHR_GRID;


double **GRID;     // grid points
double ****QPROB;  // genotype probabilities for each individual at each grid

int TNGRID;        // total number of grids
int *NGRID;        // number of grids at each chromosome
int CHL;           // max number of grids 

//*****************************************************************************************
// these parameters are passed from R function bmq.mcmc

double *Y;             // phenotypic data
int CATEGORY;          // 1: normal data; 2: binary data; 3:ordinal data
int CN;                // Categories # for binary or ordinal data

int NITER;             // Number of iterations
int NTHIN;             // Thinning value
int NBURNIN;            // Burnin
int VERBOSE;           // Verbose

int GIBBS;             // 1: Gibbs scaning all effects; 0: Kohn's M-H method for MCMC algorithm
int UPDATEGENO;        // 1: update QTL genotypes; 0: doesn't update QTL genotype

int EPISTASIS;         // 1: epistatic model; 0: non-epistatic model;
int E_NQTL_MAIN;    // expected number of main-effect QTL 
int E_NQTL;         // expected number of all QTL 
int NQTL;              // max QTL #  

double *DQQ;            // distance between flanking two genes 
int *CHR_NQTL;         // max QTL # at each chromosome

int ENV_FACTOR;        // 1:include environmental factors-need to fix this
int GBYE;              // 1: include g by e interactions
int NRANCOVA;          // random effects #
int NFIXCOVA;          // fixed effects # 
int *GBYE_FIX_INDEX;   // indicating which fixed covariates are treated in g-by-e
int *NRAN;             // random effect #
double **COEF_RAN;     // random covariates
double **COEF_FIX;      // fixed covariates

int DEPENDENCE;        // see Chipman's paper
double *C;             // see Chipman's paper 

int SEED;		      // the pseudo-random number generator

//********************************************************************************
 
int GROUP;             // 1: groupedly update all main effects or epistatic effects
int NC;                // the number of main effects at one QTL

int SPH;               // 1: standardized phenotype; 0: original phenotype

//********************************************************************************
// for binary and ordinal traits

int *W;                // ordinal or binary phenotype
double *CUTPOINT;      // threshold values for ordinal traits

//********************************************************************************
// parameters used in prior specification

double W_MAIN;          // prior for main effect indicator
double W_EPISTASIS;     // prior for epistatic effect indicator
double W_GBYE;          // prior for g by e indicator
double  *VMAIN;         // prior variance of main effects
double **VEPISTASIS;    // prior variance of epistatic effects
double **V_GBYE_FIX;    // prior variance of g-by-e effects

//*********************************************************************************
// genetic model parameters

double AMU;             // overall mean
double VE;              // residual variance
double  **MAIN;         // main effects
double ****EPISTATIC;   // epistatic effects

double *GVALUE;         // genotypic valyes
int **GENO;             // QTL genotypes
double ***COEF;         // coefficients of QTL main effects

//*********************************************************************************
// QTL positions, genetic effects indicators

int *QLOC;             // QTL position indicators, position is GRID[QCHR[L]][QLOC[L]] 
int *QCHR;             // chromosomes that QTL locate
int *CHRQTL;           // QTL number at each chromosome 

int *GAMMA;            // QTL indicators
int *GAMMA_MAIN;       // main effects indicators
int **GAMMA_EPISTASIS; // epistatic effects indicators 

//**********************************************************************************
// environmental covariate parameters

double  *FIX;         // effects of fixed covariates 
double  **RAN;        // effects of random covariates
double  *VRAN;        // variances of random covariates
double ***GBYE_FIX;   // interactions of QTL main effects and fixed covariates
double  **GAMMA_GBYE; // g-by-e indicators

//***********************************************************************************

double PDD1, PDD2;
double *PD1, *PD2;

int IBD;

double  *X;

//*************************************************************************************************************

char iterfile[100];
char pairfile[100];
char mainfile[100];
char gbyefile[100];
char covfile[100];


//*******************************************************   
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
double NormalFunction(double X)  //Calculate the normal function value
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
double TrunNormal(int WW,double B,double V) //Sample from a double-truncated normal density
{
	double u,p1,p2,p,alpha=0,y=0,yy=0;
	double c0=2.515517,c1=0.802853,c2=0.010328,d1=1.432788,d2=0.189269,d3=0.001308;
	double T1=CUTPOINT[WW],T2=CUTPOINT[WW+1];
	
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

//*************************************************************************************
double Likelihood(double *p,double *G) //Calculate the likelihood for categorical trait
{
	double CC[CN+1], AL=0.0; int I,J;
	
    for(I=0;I<NS;I++)            
	{
		for(J=0;J<=CN;J++) CC[J]=NormalFunction((p[J]-AMU-G[I])/sqrt(VE));
		double T=0.0;
		for(J=1;J<=CN;J++) T=T+(W[I]==J-1)*(CC[J]-CC[J-1]);
		AL=AL+log(T+1e-20);
	}

	return(AL);
}
//****************************************************************
//Transfering genotypes to regression coefficients

void Coefficient(int GENOTYPE)          //UPDATEGENO=1
{
	int K;

	if(GROUP==1)
	{
		for(K=0;K<NG;K++) X[K]=(GENOTYPE==K)*1.0;
	}

	if(GROUP==0)
	{
		if(CROSS==2)                 //Cockerham model
		{
			X[0]=GENOTYPE-1.0,
			X[1]=GENOTYPE*(2.0-GENOTYPE)-0.5;
		}
		else X[0]=GENOTYPE-0.5;
	}

	return;
}

void Coefficient0(int I,int L,int QL)     //UPDATEGENO=0
{ 
	if(CROSS==2)                //Cockerham model
	{
		X[0]=QPROB[QCHR[L]][I][QL][2]-QPROB[QCHR[L]][I][QL][0],
		X[1]=0.5*(QPROB[QCHR[L]][I][QL][1]-QPROB[QCHR[L]][I][QL][0]-QPROB[QCHR[L]][I][QL][2]);
	}
	else X[0]=0.5*(QPROB[QCHR[L]][I][QL][1]-QPROB[QCHR[L]][I][QL][0]);

	return;
}

//******************************************************************************************
//used in updating QTL genotypes and locations.

double GenotypeSampling(int I,int L,int II,int QL)   //if UPDATEGENO==1(0), dosenot need QL(II) 
{
    int L1,K,K1,K2; double G;
	
	if(UPDATEGENO==1) Coefficient(II);
	if(UPDATEGENO==0) Coefficient0(I,L,QL);
	G=GVALUE[I];
	for(K=0;K<NC;K++) G=G-COEF[I][L][K]*MAIN[L][K]+X[K]*MAIN[L][K];
			
	if(EPISTASIS==1)
		for(L1=0;L1<NQTL;L1++)
		{
			if(L1<L&&GAMMA[L1]!=0&&GAMMA_EPISTASIS[L1][L]!=0)
				for(K1=0;K1<NC;K1++)
					for(K2=0;K2<NC;K2++) 
						G=G-COEF[I][L1][K1]*COEF[I][L][K2]*EPISTATIC[L1][L][K1][K2]+COEF[I][L1][K1]*X[K2]*EPISTATIC[L1][L][K1][K2];

            if(L1>L&&GAMMA[L1]!=0&&GAMMA_EPISTASIS[L][L1]!=0) 
				for(K1=0;K1<NC;K1++)
					for(K2=0;K2<NC;K2++) 
						G=G-COEF[I][L][K1]*COEF[I][L1][K2]*EPISTATIC[L][L1][K1][K2]+X[K1]*COEF[I][L1][K2]*EPISTATIC[L][L1][K1][K2];
		}

	if(GBYE==1)
		for(L1=0;L1<NFIXCOVA;L1++) 
			if(GBYE_FIX_INDEX[L1]==1&&GAMMA_GBYE[L1][L]!=0)
				for(K=0;K<NC;K++) 
					G=G-COEF_FIX[I][L1]*COEF[I][L][K]*GBYE_FIX[L1][L][K]+COEF_FIX[I][L1]*X[K]*GBYE_FIX[L1][L][K];
		
	return(G);
} 

//*************************************************************************************
// update the overall mean

void Mean(double YBAR,double VP)     //prior for AMU is N(YBAR,VP)
{
	int I; double T=0.0,U,U0; 
	
    for(I=0;I<NS;I++) T=T+(Y[I]-GVALUE[I]);
	
	ANORMAL(&U,&U0); 
	AMU=(T+YBAR*VE/VP)/(NS+VE/VP)+U*sqrt(VE/(NS+VE/VP));  
      
	return;
}

//*******************************************************************
//update residual variance

void ResidualVariance(double VP)  //prior for VE is Inv-chisq(NU,TAU)
{	
		int I; double T1=0.0,T2=0.0,U,U0;
		
		int NU=6; double TAU=VP/3;

        for(I=0;I<NS;I++)
		{
            T1=T1+pow(Y[I]-AMU-GVALUE[I],2);
            ANORMAL(&U,&U0);
            T2=T2+U*U; 
		} 
		for(I=0;I<NU;I++)
		{
            ANORMAL(&U,&U0);
            T2=T2+U*U; 
		} 
     	
		VE=(T1+NU*TAU)/T2;       //see Gelman P301
		
		return;
}

//*******************************************************
//update marginal effects

void MainEffect(int L)
{
	int I,K; double G[NS1],T1,T2,T3,U,U0;

	for(K=0;K<NC;K++) 
	{
		if(MAIN[L][K]!=0)
		{
			T1=0,T2=0;
			for(I=0;I<NS;I++)                     
			{
				double Z=COEF[I][L][K];
				G[I]=GVALUE[I]-Z*MAIN[L][K];
				T1=T1+Z*(Y[I]-AMU-G[I]);
				T2=T2+pow(Z,2);
			}	 	
			T3=1/VMAIN[K]+T2/VE;

			ANORMAL(&U,&U0);
			MAIN[L][K]=T1/(VE*T3)+U/sqrt(T3);

			for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF[I][L][K]*MAIN[L][K];
		}
	}

	return;
}

//************************************************************
//update epistatic effects

void EpistaticEffect(int L1,int L2)
{
	int I,K1,K2; double G[NS1],T1,T2,T3,U,U0;

	for(K1=0;K1<NC;K1++)
		for(K2=0;K2<NC;K2++)
		{
			if(EPISTATIC[L1][L2][K1][K2]!=0)
			{
				T1=0,T2=0;
				for(I=0;I<NS;I++)                     
				{
					double Z=COEF[I][L1][K1]*COEF[I][L2][K2];
					G[I]=GVALUE[I]-Z*EPISTATIC[L1][L2][K1][K2];
					T1=T1+Z*(Y[I]-AMU-G[I]);
					T2=T2+pow(Z,2);
				}	 	
				T3=1/VEPISTASIS[K1][K2]+T2/VE;

				ANORMAL(&U,&U0);
				EPISTATIC[L1][L2][K1][K2]=T1/(VE*T3)+U/sqrt(T3);
 
				for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF[I][L1][K1]*COEF[I][L2][K2]*EPISTATIC[L1][L2][K1][K2];
			}
		}

	return;
}

//************************************************************
//update g by e fixed effects

void GBYE_FixedCovariate(int L1,int L2)
{
	int I,K; double G[NS1],T1,T2,T3,U,U0;

	for(K=0;K<NC;K++)
	{
		if(GBYE_FIX[L1][L2][K]!=0)
		{
			T1=0,T2=0;
			for(I=0;I<NS;I++)                     
			{
				double Z=COEF_FIX[I][L1]*COEF[I][L2][K];
				G[I]=GVALUE[I]-Z*GBYE_FIX[L1][L2][K];
				T1=T1+Z*(Y[I]-AMU-G[I]);
				T2=T2+pow(Z,2);
			}	 	
		 	T3=1/V_GBYE_FIX[L1][K]+T2/VE;

			ANORMAL(&U,&U0);
			GBYE_FIX[L1][L2][K]=T1/(VE*T3)+U/sqrt(T3);

			for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF_FIX[I][L1]*COEF[I][L2][K]*GBYE_FIX[L1][L2][K];
		}
	}

	return;
}

//***********************************************************
//updating genetic variances


void MainVariance(double VP)    //prior for VMAIN is Inv-chisq(NU,(NU-2)/NU*TAU),E(VMAIN)=TAU
{
	int L,K,J; double T1,T2,U,U0;

	int NU=6; double TAU=0.2;

	for(K=0;K<NC;K++) 
	{
		int N_MAIN=0;
		for(L=0;L<NQTL;L++)
			if(MAIN[L][K]!=0) N_MAIN=N_MAIN+1; 

		T1=0,T2=0;
		for(L=0;L<NQTL;L++)
			if(MAIN[L][K]!=0) T1=T1+pow(MAIN[L][K],2);

		for(J=0;J<N_MAIN+NU;J++)
		{
			ANORMAL(&U,&U0);
			T2=T2+U*U; 
		} 
     	
		T1=T1/(VMAIN[K]*VP);T2=T2/(VMAIN[K]*VP);
		VMAIN[K]=(T1+NU*(NU-2.0)/NU*TAU)/T2;
	}

	return;
}


void EpistaticVariance(double VP)    //prior for VEPISTASIS is Inv-chisq(NU,(NU-2)/NU*TAU),E(VEPISTASIS)=TAU
{
	int L1,L2,K1,K2,J;double T1,T2,U,U0;

	int NU=6; double TAU=0.1;

	for(K1=0;K1<NC;K1++)
		for(K2=0;K2<NC;K2++)
		{
			int N_EPIS=0;
			for(L1=0;L1<NQTL;L1++)
				for(L2=0;L2<NQTL;L2++)
					if(EPISTATIC[L1][L2][K1][K2]!=0) N_EPIS=N_EPIS+1; 
			
			T1=0,T2=0;
			for(L1=0;L1<NQTL;L1++)
				for(L2=L1+1;L2<NQTL;L2++)
					if(EPISTATIC[L1][L2][K1][K2]!=0) T1=T1+pow(EPISTATIC[L1][L2][K1][K2],2);

			for(J=0;J<N_EPIS+NU;J++)
			{
				ANORMAL(&U,&U0);
				T2=T2+U*U; 
			} 
     	
			T1=T1/(VEPISTASIS[K1][K2]*VP);T2=T2/(VEPISTASIS[K1][K2]*VP);
			VEPISTASIS[K1][K2]=(T1+NU*(NU-2.0)/NU*TAU)/T2;
		}

	return;
}  


//updating genetic variances

void GBYE_FixedCovariate_Variance(int L1,double VP)  //prior for V_GBYE is Inv-chisq(NU,(NU-2)/NU*TAU),E(V_GBYE)=TAU
{
	int L2,K,J; double T1,T2,U,U0;

	int NU=6; double TAU=0.1;

	for(K=0;K<NC;K++) 
	{
		int N_GBYE=0; 
		for(L2=0;L2<NQTL;L2++)
			if(GBYE_FIX[L1][L2][K]!=0) N_GBYE=N_GBYE+1;

		T1=0,T2=0;
		for(L2=0;L2<NQTL;L2++)
			if(GBYE_FIX[L1][L2][K]!=0) T1=T1+pow(GBYE_FIX[L1][L2][K],2);

		for(J=0;J<N_GBYE+NU;J++)
		{
			ANORMAL(&U,&U0);
			T2=T2+U*U; 
		} 
     	
		T1=T1/(V_GBYE_FIX[L1][K]*VP);T2=T2/(V_GBYE_FIX[L1][K]*VP);
		V_GBYE_FIX[L1][K]=(T1+NU*(NU-2.0)/NU*TAU)/T2;
	}

	return;
}
        

//*******************************************************
//update nongenetic effects

void FixedCovariate(int L)
{
	int I; double G[NS1],Z,T1=0.0,T2=0.0,T3,U,U0;

    for(I=0;I<NS;I++)
    {			
		Z=COEF_FIX[I][L];
		G[I]=GVALUE[I]-Z*FIX[L]; 
		T1=T1+Z*(Y[I]-AMU-G[I]);
        T2=T2+Z*Z;
    }
	T3=T2/VE;	               //prior is uniform

	ANORMAL(&U,&U0); 	
	FIX[L]=T1/(VE*T3)+U/sqrt(T3);
        
	for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF_FIX[I][L]*FIX[L];

	return;
}

void RandomCovariate(int L)
{
	int I,K; double G[NS1],T1,T2,T3,U,U0;
		
	for(I=0;I<NS;I++) G[I]=GVALUE[I]-RAN[L][(int)COEF_RAN[I][L]];

	for(K=0;K<NRAN[L];K++)
	{
		T1=0,T2=0;
		for(I=0;I<NS;I++)                    
			if(COEF_RAN[I][L]==K)
			{
				T1=T1+(Y[I]-AMU-G[I]);
				T2=T2+1.0;
			}	 	
		T3=1/VRAN[L]+T2/VE;

		ANORMAL(&U,&U0);
		RAN[L][K]=T1/(VE*T3)+U/sqrt(T3);
	}	
	for(I=0;I<NS;I++) GVALUE[I]=G[I]+RAN[L][(int)COEF_RAN[I][L]];
		
	return;	
}

void RanVariance(int L,double VP)         //prior is Inv-chisq(NU,TAU)
{
	int K; double T1=0,T2=0,U,U0;

	int NU=6; double TAU=VP/3;

	for(K=0;K<NRAN[L];K++) T1=T1+pow(RAN[L][K],2);
	for(K=0;K<NRAN[L]+NU;K++)
	{
        ANORMAL(&U,&U0);
        T2=T2+U*U; 
	} 
     	
	VRAN[L]=(T1+NU*TAU)/T2;      //see Gelman P301
	
	return;
}

//***********************************************************
//update QTL genotypes : generate IBD

void QTLgenotype(int L,int NL,int QL,int I)
{
	int K; double SUMM[NG];              

	for(K=0;K<NG;K++) SUMM[K]=GenotypeSampling(I,L,K,QL); 
              
    double SUM=0.0;
    for(K=0;K<NG;K++)
    {
	   SUMM[K]=exp(-0.5*pow(Y[I]-(AMU+SUMM[K]),2)/VE)*QPROB[NL][I][QL][K];
       SUM=SUM+SUMM[K]; 
    }
	for(K=0;K<NG;K++) SUMM[K]=SUMM[K]/SUM; 

    MULTINORMAL(SUMM);

	PDD1=0.0,PDD2=0.0;
    for(K=0;K<NG;K++)
    {
	   PDD1=PDD1+QPROB[NL][I][QL][K]*(IBD==K);
       PDD2=PDD2+SUMM[K]*(IBD==K);
    }     
	return;
}

//*******************************************************
//UPDATE THE QTL INHERITANCE OF NON-FOUNDERS

void QTLINHERITANCE(int L,int I)
{   
	int K;

	QTLgenotype(L,QCHR[L],QLOC[L],I); 

	GENO[I][L]=IBD;
	Coefficient(IBD);
	for(K=0;K<NC;K++) COEF[I][L][K]=X[K];
    
	PD1[L]=PD1[L]+log(PDD1+1e-20);
    PD2[L]=PD2[L]+log(PDD2+1e-20);		    

	return;
}

//*******************************************************
//UPDATE POSITIONS OF QTLS

void QTLPOSITION(int L,int QLNEW)
{
	int I,GENO1[NS1],K; double G[NS1];
   
	double PROB0=0.0;
    for(I=0;I<NS;I++) PROB0=PROB0-0.5*pow(Y[I]-AMU-GVALUE[I],2)/VE;
               
    double PROB1=0.0,PD10=0.0,PD20=0.0;       
    for(I=0;I<NS;I++)
	{
		if(UPDATEGENO==1)
		{
			int KK=0;
			for(K=0;K<NG;K++)
				if(QPROB[QCHR[L]][I][QLNEW][K]>0.99) 
				{
					GENO1[I]=K;
					KK=1;
				}    
			if(KK==0)
			{
				QTLgenotype(L,QCHR[L],QLNEW,I);                             
           
				PD10=PD10+log(PDD1+1e-20);
				PD20=PD20+log(PDD2+1e-20);

				GENO1[I]=IBD;
			}
			int II=GENO1[I];
			G[I]=GenotypeSampling(I,L,II,QLNEW);  //QLNEW is not used
		}

		if(UPDATEGENO==0) G[I]=GenotypeSampling(I,L,0,QLNEW);  //0 is not used

        PROB1=PROB1-0.5*pow(Y[I]-AMU-G[I],2)/VE;               
	}   

    double S1=(PROB1-PROB0)+(PD10-PD1[L])+(PD2[L]-PD20);
          
    if(S1>log(RANDOM())) 
	{
		QLOC[L]=QLNEW;
		                   
        for(I=0;I<NS;I++)          
		{
			GVALUE[I]=G[I];				
			if(UPDATEGENO==1)
			{
				GENO[I][L]=GENO1[I];
				Coefficient(GENO1[I]);
			}
            if(UPDATEGENO==0) Coefficient0(I,L,QLOC[L]);
			for(K=0;K<NC;K++) COEF[I][L][K]=X[K];
		}		
			
		PD1[L]=PD10,PD2[L]=PD20;
	}
		
	return;
}                              
                     
//*************************************************************
//Sampling a position

int SamplingOnePosition(int L)
{	
	int NL,I,K,CHR0[NLG],T,TT,TTT,GRID00[TNGRID]; double R,GRID0[TNGRID];                                   
	
	for(I=0;I<NLG;I++)
	{
		CHRQTL[I]=0;
		for(K=0;K<NQTL;K++) CHRQTL[I]=CHRQTL[I]+(QCHR[K]==I&&GAMMA[K]!=0);
	}

	T=0;
    for(I=0;I<NLG;I++) 
	{
		if(GRID[I][NGRID[I]-1]!=0&&DQQ[I]!=0)
			CHR0[I]=NGRID[I]-(int)(2*DQQ[I]*CHRQTL[I]*NGRID[I]/GRID[I][NGRID[I]-1]);
		if(DQQ[I]==0) CHR0[I]=NGRID[I]-CHRQTL[I];
		if(GRID[I][NGRID[I]-1]==0) CHR0[I]=0;
		if(CHRQTL[I]>=CHR_NQTL[I]||CHR0[I]<=0) CHR0[I]=0;
		T=T+CHR0[I];
	}

	if(T!=0)
	{
		R=T*RANDOM();
		TT=0;
        for(I=0;I<NLG;I++) 
		{
			if(R>=TT&&R<(TT+CHR0[I])) QCHR[L]=I;
			TT=TT+CHR0[I];
		}		
		NL=QCHR[L];
		if(CHRQTL[NL]+1>CHR_NQTL[NL]||CHR0[NL]==0) T=0;
	}

	if(T!=0)
	{
		if(CHRQTL[NL]==0) QLOC[L]=(int)(RANDOM()*NGRID[NL]);
		if(CHRQTL[NL]>0)
		{
			TT=0;
			for(K=0;K<NQTL;K++)
				if(QCHR[K]==NL&&GAMMA[K]!=0) 
				{
					GRID0[TT]=GRID[NL][QLOC[K]];
					TT=TT+1;
				}

			QLOC[L]=(int)(RANDOM()*NGRID[NL]);
			TT=0;
			for(K=0;K<CHRQTL[NL];K++)
				if(fabs(GRID[NL][QLOC[L]]-GRID0[K])<=DQQ[NL]) TT=TT+1;
	
			if(TT!=0)
			{
				TTT=0;
				for(I=0;I<NGRID[NL];I++)
				{
					GRID00[I]=1;
					for(K=0;K<CHRQTL[NL];K++)
						if(fabs(GRID[NL][I]-GRID0[K])<=DQQ[NL]) GRID00[I]=0;
					TTT=TTT+GRID00[I];
				}
				if(TTT!=0)
				{
					R=RANDOM()*TTT;
					TT=0;
					for(I=0;I<NGRID[NL];I++)
					{
						if(R>=TT&&R<(TT+GRID00[I])) QLOC[L]=I;
						TT=TT+GRID00[I];
					}	
				}
				if(TTT==0) T=0;
			}

			for(K=0;K<CHRQTL[NL];K++)
				if(fabs(GRID[NL][QLOC[L]]-GRID0[K])<=DQQ[NL]) T=0;
		}
	}

	if(T!=0)
	{
		for(I=0;I<NS;I++)         
		{
			if(UPDATEGENO==1)
			{
				double PRR[NG];
				for(K=0;K<NG;K++) PRR[K]=QPROB[NL][I][QLOC[L]][K];
				MULTINORMAL(PRR);                                     
				GENO[I][L]=IBD;
				Coefficient(IBD);
			}
			if(UPDATEGENO==0) Coefficient0(I,L,QLOC[L]);
			for(K=0;K<NC;K++) COEF[I][L][K]=X[K];
		}
	}

	return(T);
}

//***************************************************************
//Deleting QTL with all 0 effects

void ZeroEffect()
{
	int L1,L2; 

    for(L1=0;L1<NQTL;L1++) GAMMA[L1]=0;

    for(L1=0;L1<NQTL;L1++)
	{
		GAMMA[L1]=GAMMA[L1]+(GAMMA_MAIN[L1]!=0);
    
	    if(EPISTASIS==1)
		{
			for(L2=L1+1;L2<NQTL;L2++) GAMMA[L1]=GAMMA[L1]+(GAMMA_EPISTASIS[L1][L2]!=0);   
			
			for(L2=0;L2<=L1-1;L2++) GAMMA[L1]=GAMMA[L1]+(GAMMA_EPISTASIS[L2][L1]!=0);
		}

		if(GBYE==1)
			for(L2=0;L2<NFIXCOVA;L2++) 
				if(GBYE_FIX_INDEX[L2]==1) GAMMA[L1]=GAMMA[L1]+(GAMMA_GBYE[L2][L1]!=0);
	}

	return;
}

void ZeroEffect1(int L)
{
	int L0; 

	GAMMA[L]=(GAMMA_MAIN[L]!=0);

	if(EPISTASIS==1)
	{
		for(L0=0;L0<=L-1;L0++) GAMMA[L]=GAMMA[L]+(GAMMA_EPISTASIS[L0][L]!=0);
		for(L0=L+1;L0<NQTL;L0++) GAMMA[L]=GAMMA[L]+(GAMMA_EPISTASIS[L][L0]!=0);   			
	}

	if(GBYE==1)
	{
		for(L0=0;L0<NFIXCOVA;L0++) 
			if(GBYE_FIX_INDEX[L0]==1) GAMMA[L]=GAMMA[L]+(GAMMA_GBYE[L0][L]!=0);
	}

	return;
}

void ZeroEffect2(int L1, int L2)
{
	int L0; 

	GAMMA[L1]=(GAMMA_MAIN[L1]!=0);
	GAMMA[L2]=(GAMMA_MAIN[L2]!=0);

	if(EPISTASIS==1)
	{
		for(L0=0;L0<=L1-1;L0++) GAMMA[L1]=GAMMA[L1]+(GAMMA_EPISTASIS[L0][L1]!=0);
		for(L0=L1+1;L0<NQTL;L0++) GAMMA[L1]=GAMMA[L1]+(GAMMA_EPISTASIS[L1][L0]!=0);
		
		for(L0=0;L0<=L2-1;L0++) GAMMA[L2]=GAMMA[L2]+(GAMMA_EPISTASIS[L0][L2]!=0);
		for(L0=L2+1;L0<NQTL;L0++) GAMMA[L2]=GAMMA[L2]+(GAMMA_EPISTASIS[L2][L0]!=0); 
	}

	if(GBYE==1)
	{
		for(L0=0;L0<NFIXCOVA;L0++) 
			if(GBYE_FIX_INDEX[L0]==1) 
			{
				GAMMA[L1]=GAMMA[L1]+(GAMMA_GBYE[L0][L1]!=0);
				GAMMA[L2]=GAMMA[L2]+(GAMMA_GBYE[L0][L2]!=0);
			}
	}

	return;
}


//************************************************************************
//Update main effect indicators

void MainEffectIndicator_GROUP0(int L,int K)    
{
	int I,K1;
	double G[NS1],BF_10,GAMMA10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100),T1=0,T2=0,T3,F1=0,F2=0;

	
	for(I=0;I<NS;I++)                     
	{
		double Z=COEF[I][L][K];
		G[I]=GVALUE[I]-Z*MAIN[L][K];
		T1=T1+Z*(Y[I]-AMU-G[I]);
		T2=T2+Z*Z;
	}	 	
	T3=1/VMAIN[K]+T2/VE;	
	
	for(I=0;I<NS;I++)                     
	{
		double Z=COEF[I][L][K];
		F1=F1+pow(Y[I]-AMU-G[I],2)/VE;
		F2=F2+pow(Y[I]-AMU-G[I]-Z*T1/(VE*T3),2)/VE;
	}

	BF_10=-0.5*log(VMAIN[K])-0.5*pow(T1/(VE*T3),2)/VMAIN[K]-0.5*log(T3);
	BF_10=-0.5*F2+0.5*F1+BF_10;      
	
	GAMMA10=log(W_MAIN)-log(1-W_MAIN);
	if(GIBBS==1) GAMMA_1=BF_10+GAMMA10-log(1+exp(BF_10+GAMMA10));
	if(GIBBS==0&&MAIN[L][K]==0) GAMMA_1=BF_10;
	if(GIBBS==0&&MAIN[L][K]!=0) GAMMA_0=-BF_10;  
	      
	double R=log(RANDOM());
	if(R<GAMMA_1) 
	{
		double U,U0;
		ANORMAL(&U,&U0);
		MAIN[L][K]=T1/(VE*T3)+U/sqrt(T3);
		GAMMA_MAIN[L]=1;
		GAMMA[L]=1;
		for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF[I][L][K]*MAIN[L][K];
	}
	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
	{
		MAIN[L][K]=0;
		double SUM=0;
		for(K1=0;K1<NC;K1++) SUM=SUM+fabs(MAIN[L][K1]);
		if(SUM==0) 
		{
			GAMMA_MAIN[L]=0; 
			ZeroEffect1(L);
		}
		for(I=0;I<NS;I++) GVALUE[I]=G[I];
	}

	return;
}


void MainEffectIndicator_GROUP1(int L)   
{	
	int I,K;
	double G[NS1],BF_10,GAMMA10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100),T1[NG],T2[NG],T3[NG],F1=0,F2=0;


	for(I=0;I<NG;I++)
	{
		T1[I]=0, T2[I]=0, T3[I]=0;
	}

	for(I=0;I<NS;I++)
	{
		G[I]=GVALUE[I];
		for(K=0;K<NC;K++) G[I]=G[I]-COEF[I][L][K]*MAIN[L][K];
	}

	for(K=0;K<NC;K++)
	{
		for(I=0;I<NS;I++)                     
		{
			double Z=COEF[I][L][K];
			T1[K]=T1[K]+Z*(Y[I]-AMU-G[I]);
			T2[K]=T2[K]+Z*Z;
		}	 	
		T3[K]=1/VMAIN[K]+T2[K]/VE;
	}	
	
	for(I=0;I<NS;I++)                     
	{
		F1=F1+pow(Y[I]-AMU-G[I],2)/VE;
		double F3=0;
		for(K=0;K<NC;K++) F3=F3+COEF[I][L][K]*T1[K]/(VE*T3[K]);
		F2=F2+pow(Y[I]-AMU-G[I]-F3,2)/VE;
	}

	BF_10=0;
	for(K=0;K<NC;K++) BF_10=BF_10-0.5*log(VMAIN[K])-0.5*pow(T1[K]/(VE*T3[K]),2)/VMAIN[K]-0.5*log(T3[K]);
	BF_10=-0.5*F2+0.5*F1+BF_10;       
	
	GAMMA10=log(W_MAIN)-log(1-W_MAIN);
	if(GIBBS==1) GAMMA_1=BF_10+GAMMA10-log(1+exp(BF_10+GAMMA10));
	if(GIBBS==0&&GAMMA_MAIN[L]==0) GAMMA_1=BF_10;
	if(GIBBS==0&&GAMMA_MAIN[L]!=0) GAMMA_0=-BF_10;  
	      
	double R=log(RANDOM());
	if(R<GAMMA_1) 
	{
		for(K=0;K<NC;K++) 
		{
			double U,U0;
			ANORMAL(&U,&U0);
			MAIN[L][K]=T1[K]/(VE*T3[K])+U/sqrt(T3[K]);
		}
		GAMMA_MAIN[L]=1;
		GAMMA[L]=1;
		for(I=0;I<NS;I++)
		{
			GVALUE[I]=G[I];
			for(K=0;K<NC;K++) GVALUE[I]=GVALUE[I]+COEF[I][L][K]*MAIN[L][K];
		}
	}
	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
	{
		for(K=0;K<NC;K++) MAIN[L][K]=0;
		GAMMA_MAIN[L]=0;
		for(I=0;I<NS;I++) GVALUE[I]=G[I];
		ZeroEffect1(L);
	}

	return;
}

//************************************************************************
//Update epistatic effect indicators


void EpistasisIndicator_GROUP0(int L1,int L2,int K1,int K2)  
{
	int I,K01,K02;
	double G[NS1],BF_10,GAMMA10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100),T1=0,T2=0,T3,F1=0,F2=0;

	if(W_EPISTASIS!=0)
	{
	for(I=0;I<NS;I++)                     
	{
		double Z=COEF[I][L1][K1]*COEF[I][L2][K2];
		G[I]=GVALUE[I]-Z*EPISTATIC[L1][L2][K1][K2];
		T1=T1+Z*(Y[I]-AMU-G[I]);
		T2=T2+Z*Z;
	}	 	
	T3=1/VEPISTASIS[K1][K2]+T2/VE;	
	
	for(I=0;I<NS;I++)                     
	{
		double Z=COEF[I][L1][K1]*COEF[I][L2][K2];
		F1=F1+pow(Y[I]-AMU-G[I],2)/VE;
		F2=F2+pow(Y[I]-AMU-G[I]-Z*T1/(VE*T3),2)/VE;
	}

	BF_10=-0.5*log(VEPISTASIS[K1][K2])-0.5*pow(T1/(VE*T3),2)/VEPISTASIS[K1][K2]-0.5*log(T3);
	BF_10=-0.5*F2+0.5*F1+BF_10;      

	GAMMA10=log(W_EPISTASIS)-log(1-W_EPISTASIS);
	if(GIBBS==1) GAMMA_1=BF_10+GAMMA10-log(1+exp(BF_10+GAMMA10));
	if(GIBBS==0&&EPISTATIC[L1][L2][K1][K2]==0) GAMMA_1=BF_10;
	if(GIBBS==0&&EPISTATIC[L1][L2][K1][K2]!=0) GAMMA_0=-BF_10;

	double R=log(RANDOM());
	if(R<GAMMA_1) 
	{
		double U,U0;
		ANORMAL(&U,&U0);
		EPISTATIC[L1][L2][K1][K2]=T1/(VE*T3)+U/sqrt(T3);
		GAMMA_EPISTASIS[L1][L2]=1;
		GAMMA[L1]=1,GAMMA[L2]=1;
		for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF[I][L1][K1]*COEF[I][L2][K2]*EPISTATIC[L1][L2][K1][K2];
	}
	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
	{
		EPISTATIC[L1][L2][K1][K2]=0;
		double SUM=0;
		for(K01=0;K01<NC;K01++) 
			for(K02=0;K02<NC;K02++) SUM=SUM+fabs(EPISTATIC[L1][L2][K01][K02]);
		if(SUM==0) 
		{
			GAMMA_EPISTASIS[L1][L2]=0;  
			ZeroEffect2(L1,L2);
		}
		for(I=0;I<NS;I++) GVALUE[I]=G[I];
	}

	}

	if(W_EPISTASIS==0&&EPISTATIC[L1][L2][K1][K2]!=0)
	{
		for(I=0;I<NS;I++) GVALUE[I]=GVALUE[I]-COEF[I][L1][K1]*COEF[I][L2][K2]*EPISTATIC[L1][L2][K1][K2];
		EPISTATIC[L1][L2][K1][K2]=0;
		double SUM=0;
		for(K01=0;K01<NC;K01++) 
			for(K02=0;K02<NC;K02++) SUM=SUM+fabs(EPISTATIC[L1][L2][K01][K02]);
		if(SUM==0) 
		{
			GAMMA_EPISTASIS[L1][L2]=0;  
			ZeroEffect2(L1,L2);
		}
	}


	return;	
}


void EpistasisIndicator_GROUP1(int L1,int L2)    
{
	int I,K1,K2;
	double G[NS1],BF_10,GAMMA10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100),T1[NG][NG],T2[NG][NG],T3[NG][NG],F1=0,F2=0;

	if(W_EPISTASIS!=0)
	{

	for(K1=0;K1<NG;K1++)
	   for(K2=0;K2<NG;K2++)
	   {
           T1[K1][K2]=0, T2[K1][K2]=0, T3[K1][K2]=0;
	   }

	for(I=0;I<NS;I++)
	{
		G[I]=GVALUE[I];
		for(K1=0;K1<NC;K1++) 
			for(K2=0;K2<NC;K2++) G[I]=G[I]-COEF[I][L1][K1]*COEF[I][L2][K2]*EPISTATIC[L1][L2][K1][K2];
	}

	for(K1=0;K1<NC;K1++)
		for(K2=0;K2<NC;K2++)
		{
			for(I=0;I<NS;I++)                     
			{
				double Z=COEF[I][L1][K1]*COEF[I][L2][K2];
				T1[K1][K2]=T1[K1][K2]+Z*(Y[I]-AMU-G[I]);
				T2[K1][K2]=T2[K1][K2]+Z*Z;
			}	 	
			T3[K1][K2]=1/VEPISTASIS[K1][K2]+T2[K1][K2]/VE;
		}	
	
	for(I=0;I<NS;I++)                     
	{
		F1=F1+pow(Y[I]-AMU-G[I],2)/VE;
		double F3=0;
		for(K1=0;K1<NC;K1++) 
			for(K2=0;K2<NC;K2++) F3=F3+COEF[I][L1][K1]*COEF[I][L2][K2]*T1[K1][K2]/(VE*T3[K1][K2]);
		F2=F2+pow(Y[I]-AMU-G[I]-F3,2)/VE;
	}

	BF_10=0;
	for(K1=0;K1<NC;K1++) 
		for(K2=0;K2<NC;K2++) 
			BF_10=BF_10-0.5*log(VEPISTASIS[K1][K2])-0.5*pow(T1[K1][K2]/(VE*T3[K1][K2]),2)/VEPISTASIS[K1][K2]-0.5*log(T3[K1][K2]);
	BF_10=-0.5*F2+0.5*F1+BF_10;      

	GAMMA10=log(W_EPISTASIS)-log(1-W_EPISTASIS);
	if(GIBBS==1) GAMMA_1=BF_10+GAMMA10-log(1+exp(BF_10+GAMMA10));
	if(GIBBS==0&&GAMMA_EPISTASIS[L1][L2]==0) GAMMA_1=BF_10;
	if(GIBBS==0&&GAMMA_EPISTASIS[L1][L2]!=0) GAMMA_0=-BF_10;
	      
	double R=log(RANDOM());
	if(R<GAMMA_1) 
	{
		for(K1=0;K1<NC;K1++) 
			for(K2=0;K2<NC;K2++)  
			{
				double U,U0;
				ANORMAL(&U,&U0);
				EPISTATIC[L1][L2][K1][K2]=T1[K1][K2]/(VE*T3[K1][K2])+U/sqrt(T3[K1][K2]);
			}
		GAMMA_EPISTASIS[L1][L2]=1;
		GAMMA[L1]=1,GAMMA[L2]=1;
		for(I=0;I<NS;I++)
		{
			GVALUE[I]=G[I];
			for(K1=0;K1<NC;K1++) 
				for(K2=0;K2<NC;K2++) GVALUE[I]=GVALUE[I]+COEF[I][L1][K1]*COEF[I][L2][K2]*EPISTATIC[L1][L2][K1][K2];
		}
	}
	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
	{
		for(K1=0;K1<NC;K1++) 
			for(K2=0;K2<NC;K2++) EPISTATIC[L1][L2][K1][K2]=0;
		GAMMA_EPISTASIS[L1][L2]=0;
		for(I=0;I<NS;I++) GVALUE[I]=G[I];
		ZeroEffect2(L1,L2);
	}

	}

	if(W_EPISTASIS==0&&GAMMA_EPISTASIS[L1][L2]!=0)
	{
		for(I=0;I<NS;I++)
			for(K1=0;K1<NC;K1++) 
				for(K2=0;K2<NC;K2++) GVALUE[I]=GVALUE[I]-COEF[I][L1][K1]*COEF[I][L2][K2]*EPISTATIC[L1][L2][K1][K2];
		for(K1=0;K1<NC;K1++) 
			for(K2=0;K2<NC;K2++) EPISTATIC[L1][L2][K1][K2]=0;
		GAMMA_EPISTASIS[L1][L2]=0;
		ZeroEffect2(L1,L2);
	}

	return;	
}


//************************************************************************
//Update g by e fixed effect indicators


void GBYE_FIX_Indicator_GROUP0(int L1,int L2,int K)  //use the constraint model
{
	int I,K1;
	double G[NS1],BF_10,GAMMA10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100),T1=0,T2=0,T3,F1=0,F2=0;

	
	for(I=0;I<NS;I++)                     
	{
		double Z=COEF_FIX[I][L1]*COEF[I][L2][K];
		G[I]=GVALUE[I]-Z*GBYE_FIX[L1][L2][K];
		T1=T1+Z*(Y[I]-AMU-G[I]);
		T2=T2+Z*Z;
	}	 	
	T3=1/V_GBYE_FIX[L1][K]+T2/VE;	
	
	for(I=0;I<NS;I++)                     
	{
		double Z=COEF_FIX[I][L1]*COEF[I][L2][K];
		F1=F1+pow(Y[I]-AMU-G[I],2)/VE;
		F2=F2+pow(Y[I]-AMU-G[I]-Z*T1/(VE*T3),2)/VE;
	}

	BF_10=-0.5*log(V_GBYE_FIX[L1][K])-0.5*pow(T1/(VE*T3),2)/V_GBYE_FIX[L1][K]-0.5*log(T3);
	BF_10=-0.5*F2+0.5*F1+BF_10;      

	GAMMA10=log(W_GBYE)-log(1-W_GBYE);
	if(GIBBS==1) GAMMA_1=BF_10+GAMMA10-log(1+exp(BF_10+GAMMA10));
	if(GIBBS==0&&GBYE_FIX[L1][L2][K]==0) GAMMA_1=BF_10;
	if(GIBBS==0&&GBYE_FIX[L1][L2][K]!=0) GAMMA_0=-BF_10;

	double R=log(RANDOM());
	if(R<GAMMA_1) 
	{
		double U,U0;
		ANORMAL(&U,&U0);
		GBYE_FIX[L1][L2][K]=T1/(VE*T3)+U/sqrt(T3);
		GAMMA_GBYE[L1][L2]=1;
		GAMMA[L2]=1;
		for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF_FIX[I][L1]*COEF[I][L2][K]*GBYE_FIX[L1][L2][K];
	}
	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
	{
		GBYE_FIX[L1][L2][K]=0;
		double SUM=0;
		for(K1=0;K1<NC;K1++) SUM=SUM+fabs(GBYE_FIX[L1][L2][K1]);
		if(SUM==0) 
		{
			GAMMA_GBYE[L1][L2]=0;  
			ZeroEffect1(L2);
		}
		for(I=0;I<NS;I++) GVALUE[I]=G[I];
	}

	return;	
}

void GBYE_FIX_Indicator_GROUP1(int L1,int L2)  
{
	int I,K;
	double G[NS1],BF_10,GAMMA10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100),T1[NG],T2[NG],T3[NG],F1=0,F2=0;

	for(I=0;I<NG;I++)
	{
		T1[I]=0, T2[I]=0, T3[I]=0;
	}

	for(I=0;I<NS;I++)
	{
		G[I]=GVALUE[I];
		for(K=0;K<NC;K++) 
			G[I]=G[I]-COEF_FIX[I][L1]*COEF[I][L2][K]*GBYE_FIX[L1][L2][K];
	}

	for(K=0;K<NC;K++)
	{
		for(I=0;I<NS;I++)                     
		{
			double Z=COEF_FIX[I][L1]*COEF[I][L2][K];
			T1[K]=T1[K]+Z*(Y[I]-AMU-G[I]);
			T2[K]=T2[K]+Z*Z;
		}	 	
		T3[K]=1/V_GBYE_FIX[L1][K]+T2[K]/VE;
	}
	
	for(I=0;I<NS;I++)                     
	{
		F1=F1+pow(Y[I]-AMU-G[I],2)/VE;
		double F3=0;
        for(K=0;K<NC;K++) F3=F3+COEF_FIX[I][L1]*COEF[I][L2][K]*T1[K]/(VE*T3[K]);
		F2=F2+pow(Y[I]-AMU-G[I]-F3,2)/VE;
	}

	BF_10=0;
	for(K=0;K<NC;K++) BF_10=BF_10-0.5*log(V_GBYE_FIX[L1][K])-0.5*pow(T1[K]/(VE*T3[K]),2)/V_GBYE_FIX[L1][K]-0.5*log(T3[K]);
	BF_10=-0.5*F2+0.5*F1+BF_10;      

	GAMMA10=log(W_GBYE)-log(1-W_GBYE);
	if(GIBBS==1) GAMMA_1=BF_10+GAMMA10-log(1+exp(BF_10+GAMMA10));
	if(GIBBS==0&&GAMMA_GBYE[L1][L2]==0) GAMMA_1=BF_10;
	if(GIBBS==0&&GAMMA_GBYE[L1][L2]!=0) GAMMA_0=-BF_10;

	double R=log(RANDOM());
	if(R<GAMMA_1) 
	{
		for(K=0;K<NC;K++)
		{
			double U,U0;
			ANORMAL(&U,&U0);
			GBYE_FIX[L1][L2][K]=T1[K]/(VE*T3[K])+U/sqrt(T3[K]);
		}
		GAMMA_GBYE[L1][L2]=1;
		GAMMA[L2]=1;
		for(I=0;I<NS;I++) 
		{
			GVALUE[I]=G[I];
			for(K=0;K<NC;K++) GVALUE[I]=GVALUE[I]+COEF_FIX[I][L1]*COEF[I][L2][K]*GBYE_FIX[L1][L2][K];
		}
	}
	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
	{
		for(K=0;K<NC;K++) GBYE_FIX[L1][L2][K]=0;
		GAMMA_GBYE[L1][L2]=0;  
		for(I=0;I<NS;I++) GVALUE[I]=G[I];
		ZeroEffect1(L2);
	}

	return;	
}
  
//*******************************************************************************    


// ************************************************************

void bmqAnalysis()
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
	 
//	if(SPH==1)
//		for(I=0;I<NS;I++) Y[I]=(Y[I]-YBAR)/sqrt(VP);           //standardizing the phenotype

	AMU=YBAR, VE=VP;            //initial values for AMU and VE
	for(I=0;I<NRANCOVA;I++) VRAN[I]=VP;  //initial values for VRAN[I]
}

if(CATEGORY!=1)
{
	AMU=0, VE=1.0;                        //initial values for AMU and VE
	for(I=0;I<NRANCOVA;I++) VRAN[I]=1.0;  //initial values for VRAN[I]
}

//******************************************************************************
//For specify the prior variances of QTL effects and g by e interactions

double CC[NG];
if(CROSS==2)
{
	CC[0]=1.0/2, CC[1]=1.0/4;
}
else CC[0]=1.0/4;

double V_FIX[20];
for(L=0;L<NFIXCOVA;L++)
	if(GBYE_FIX_INDEX[L]==1) 
	{
		double Y2BAR=0.0;
		for(I=0;I<NS;I++) Y2BAR=Y2BAR+pow(COEF_FIX[I][L],2);
		V_FIX[L]=Y2BAR/NS;
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

// ***********************************************************                         
// ITERATION STARTS HERE
				
int ITER,ITER1;
for(ITER=0;ITER<NITER+(int)(1.0*NBURNIN/NTHIN);ITER++)
{ for(ITER1=0;ITER1<NTHIN;ITER1++)
{

//***********************************************************
//UPDATING THE VALUES OF THE LIABILITY

if(CATEGORY!=1) 
	for(I=0;I<NS;I++) Y[I]=TrunNormal(W[I],AMU+GVALUE[I],VE);      

//********************************************************************************* 
//UPDATE PARAMETERS

Mean(YBAR,VP);

if(CATEGORY!=2) ResidualVariance(VP);


for(L=0;L<NQTL;L++)
	if(GAMMA_MAIN[L]!=0.0) MainEffect(L);	
if(EPISTASIS==1)
	for(L1=0;L1<NQTL-1;L1++)
		for(L2=L1+1;L2<NQTL;L2++)
			if(GAMMA_EPISTASIS[L1][L2]!=0.0) EpistaticEffect(L1,L2);	



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
for(K=0;K<NC;K++) VMAIN[K]=1.0/CC[K];
MainVariance(VP);
if(EPISTASIS==1) 
{
	for(K1=0;K1<NC;K1++)
		for(K2=0;K2<NC;K2++) VEPISTASIS[K1][K2]=1.0/(CC[K1]*CC[K2]);
	EpistaticVariance(VP);
}
if(GBYE==1)
{  
	for(L1=0;L1<NFIXCOVA;L1++)
		if(GBYE_FIX_INDEX[L1]==1)
		{
			for(K=0;K<NC;K++) V_GBYE_FIX[L1][K]=1.0/(V_FIX[L1]*CC[K]);
			GBYE_FixedCovariate_Variance(L1,VP);
		}
}


if(ENV_FACTOR==1)
{ 
	for(L=0;L<NRANCOVA;L++) 
	{
		RandomCovariate(L);		
		RanVariance(L,VP);
	}
	for(L=0;L<NFIXCOVA;L++) FixedCovariate(L);		
}

if(GBYE==1)
{ 
	for(L1=0;L1<NFIXCOVA;L1++)
		if(GBYE_FIX_INDEX[L1]==1)
			for(L2=0;L2<NQTL;L2++)
				if(GAMMA[L2]!=0.0&&GAMMA_GBYE[L1][L2]!=0.0) GBYE_FixedCovariate(L1,L2);	
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
			int L1,K,K1,K2; double G=GVALUE[I];	
			for(K=0;K<NC;K++) G=G-COEF[I][L][K]*MAIN[L][K];
			if(EPISTASIS==1)
				for(L1=0;L1<NQTL;L1++)
				{
					if(L1<L&&GAMMA[L1]!=0&&GAMMA_EPISTASIS[L1][L]!=0)
						for(K1=0;K1<NC;K1++)
							for(K2=0;K2<NC;K2++) 
								G=G-COEF[I][L1][K1]*COEF[I][L][K2]*EPISTATIC[L1][L][K1][K2];

					if(L1>L&&GAMMA[L1]!=0&&GAMMA_EPISTASIS[L][L1]!=0) 
						for(K1=0;K1<NC;K1++)
							for(K2=0;K2<NC;K2++) 
								G=G-COEF[I][L][K1]*COEF[I][L1][K2]*EPISTATIC[L][L1][K1][K2];
				}

			if(GBYE==1)
				for(L1=0;L1<NFIXCOVA;L1++) 
					if(GBYE_FIX_INDEX[L1]==1&&GAMMA_GBYE[L1][L]!=0)
						for(K=0;K<NC;K++) G=G-COEF_FIX[I][L1]*COEF[I][L][K]*GBYE_FIX[L1][L][K];

			
			int KK=0;
			for(K=0;K<NG;K++)
			   if(QPROB[QCHR[L]][I][QLOC[L]][K]>0.99) 
			   {
				   GENO[I][L]=K;
				   Coefficient(K);
				   for(K1=0;K1<NC;K1++) COEF[I][L][K1]=X[K1];
				   KK=1;
			   } 
			if(KK==0) QTLINHERITANCE(L,I);

	
			for(K=0;K<NC;K++) G=G+COEF[I][L][K]*MAIN[L][K];
			if(EPISTASIS==1)
				for(L1=0;L1<NQTL;L1++)
				{
					if(L1<L&&GAMMA[L1]!=0&&GAMMA_EPISTASIS[L1][L]!=0)
						for(K1=0;K1<NC;K1++)
							for(K2=0;K2<NC;K2++) 
								G=G+COEF[I][L1][K1]*COEF[I][L][K2]*EPISTATIC[L1][L][K1][K2];

					if(L1>L&&GAMMA[L1]!=0&&GAMMA_EPISTASIS[L][L1]!=0) 
						for(K1=0;K1<NC;K1++)
							for(K2=0;K2<NC;K2++) 
								G=G+COEF[I][L][K1]*COEF[I][L1][K2]*EPISTATIC[L][L1][K1][K2];
				}

			if(GBYE==1)
				for(L1=0;L1<NFIXCOVA;L1++) 
					if(GBYE_FIX_INDEX[L1]==1&&GAMMA_GBYE[L1][L]!=0)
						for(K=0;K<NC;K++) G=G+COEF_FIX[I][L1]*COEF[I][L][K]*GBYE_FIX[L1][L][K];
				

			GVALUE[I]=G;
		}
	}
}      
//**********************************************************************************
//UPDATING QTL POSITIONS

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
	
//**********************************************************
// UPDATE MARGINAL EFFECT INDICATORS

for(L=0;L<NQTL;L++)
{
	int T;

	if(GROUP==0)
	{
		for(K=0;K<NC;K++)
		{
			T=1;
			if(GIBBS==1) 
			{
				if(GAMMA[L]==0) T=SamplingOnePosition(L);
				if(T!=0) MainEffectIndicator_GROUP0(L,K);
			}
			if(GIBBS==0)
			{
				double R=RANDOM();
				if( (MAIN[L][K]==0&&R<=W_MAIN)||(MAIN[L][K]!=0&&R>W_MAIN) ) 
				{
					if(GAMMA[L]==0) T=SamplingOnePosition(L);
					if(T!=0) MainEffectIndicator_GROUP0(L,K);
				}
			}
		}
	}

	if(GROUP==1)
	{
		T=1;
		if(GIBBS==1) 
		{
			if(GAMMA[L]==0) T=SamplingOnePosition(L);
			if(T!=0) MainEffectIndicator_GROUP1(L);
		}
		if(GIBBS==0)
		{
			double R=RANDOM();
			if( (GAMMA_MAIN[L]==0&&R<=W_MAIN)||(GAMMA_MAIN[L]!=0&&R>W_MAIN) ) 
			{
				if(GAMMA[L]==0) T=SamplingOnePosition(L);
				if(T!=0) MainEffectIndicator_GROUP1(L);
			}
		}
	}

}     

//**********************************************************
// UPDATE TWO ORDER EPISTATIC EFFECT INDICATORS
		
if(EPISTASIS==1)
{
	int T;
	
	for(L1=0;L1<NQTL-1;L1++)
		for(L2=L1+1;L2<NQTL;L2++)
		{
			if(GROUP==0)
			{
				for(K1=0;K1<NC;K1++)
					for(K2=0;K2<NC;K2++)
					{
						if(DEPENDENCE==1)
						{
							if(MAIN[L1][K1]!=0&&MAIN[L2][K2]!=0) W_EPISTASIS=C[0];
							if( (MAIN[L1][K1]!=0&&MAIN[L2][K2]==0)||(MAIN[L1][K1]==0&&MAIN[L2][K2]!=0) ) W_EPISTASIS=C[1];
							if(MAIN[L1][K1]==0&&MAIN[L2][K2]==0) W_EPISTASIS=C[2];
						}

						T=1;
						if(GIBBS==1) 
						{
							if(W_EPISTASIS!=0)
							{
								if(GAMMA[L1]==0&&GAMMA[L2]!=0) T=SamplingOnePosition(L1);
								if(GAMMA[L1]!=0&&GAMMA[L2]==0) T=SamplingOnePosition(L2);
								if(GAMMA[L1]==0&&GAMMA[L2]==0) 
								{
									T=SamplingOnePosition(L1);
									if(T!=0) T=SamplingOnePosition(L2);
									if(CHRQTL[QCHR[L1]]+2>CHR_NQTL[QCHR[L1]]&&QCHR[L1]==QCHR[L2]
										&&(GRID[QCHR[L1]][QLOC[L1]]-GRID[QCHR[L2]][QLOC[L2]])<=DQQ[QCHR[L1]]) T=0;
								}
							}
							if(T!=0) EpistasisIndicator_GROUP0(L1,L2,K1,K2);
						}
						if(GIBBS==0)
						{
							double R=RANDOM();
							if( (EPISTATIC[L1][L2][K1][K2]==0&&R<=W_EPISTASIS)||(EPISTATIC[L1][L2][K1][K2]!=0&&R>W_EPISTASIS) ) 
							{
								if(W_EPISTASIS!=0)
								{
									if(GAMMA[L1]==0&&GAMMA[L2]!=0) T=SamplingOnePosition(L1);
									if(GAMMA[L1]!=0&&GAMMA[L2]==0) T=SamplingOnePosition(L2);
									if(GAMMA[L1]==0&&GAMMA[L2]==0) 
									{
										T=SamplingOnePosition(L1);
										if(T!=0) T=SamplingOnePosition(L2);
										if(CHRQTL[QCHR[L1]]+2>CHR_NQTL[QCHR[L1]]&&QCHR[L1]==QCHR[L2]
											&&(GRID[QCHR[L1]][QLOC[L1]]-GRID[QCHR[L2]][QLOC[L2]])<=DQQ[QCHR[L1]]) T=0;
									}
								}
								if(T!=0) EpistasisIndicator_GROUP0(L1,L2,K1,K2);
							}
						}
					}
			}

			if(GROUP==1)
			{	
				if(DEPENDENCE==1)
				{
					if(GAMMA_MAIN[L1]!=0&&GAMMA_MAIN[L2]!=0) W_EPISTASIS=C[0];
					if( (GAMMA_MAIN[L1]!=0&&GAMMA_MAIN[L2]==0)||(GAMMA_MAIN[L1]==0&&GAMMA_MAIN[L2]!=0) ) W_EPISTASIS=C[1];
					if(GAMMA_MAIN[L1]==0&&GAMMA_MAIN[L2]==0) W_EPISTASIS=C[2];
				}
				
				T=1;
				if(GIBBS==1)
				{
					if(W_EPISTASIS!=0)
					{
						if(GAMMA[L1]==0&&GAMMA[L2]!=0) T=SamplingOnePosition(L1);
						if(GAMMA[L1]!=0&&GAMMA[L2]==0) T=SamplingOnePosition(L2);
						if(GAMMA[L1]==0&&GAMMA[L2]==0)
						{
							T=SamplingOnePosition(L1);
							if(T!=0) T=SamplingOnePosition(L2);
							if(CHRQTL[QCHR[L1]]+2>CHR_NQTL[QCHR[L1]]&&QCHR[L1]==QCHR[L2]
								&&(GRID[QCHR[L1]][QLOC[L1]]-GRID[QCHR[L2]][QLOC[L2]])<=DQQ[QCHR[L1]]) T=0;
						}
					}
					if(T!=0) EpistasisIndicator_GROUP1(L1,L2);
				}
				if(GIBBS==0)
				{
					double R=RANDOM();
					if( (GAMMA_EPISTASIS[L1][L2]==0&&R<=W_EPISTASIS)||(GAMMA_EPISTASIS[L1][L2]!=0&&R>W_EPISTASIS) )
					{
						if(W_EPISTASIS!=0)
						{
							if(GAMMA[L1]==0&&GAMMA[L2]!=0) T=SamplingOnePosition(L1);
							if(GAMMA[L1]!=0&&GAMMA[L2]==0) T=SamplingOnePosition(L2);
							if(GAMMA[L1]==0&&GAMMA[L2]==0) 
							{
								T=SamplingOnePosition(L1);
								if(T!=0) T=SamplingOnePosition(L2);
								if(CHRQTL[QCHR[L1]]+2>CHR_NQTL[QCHR[L1]]&&QCHR[L1]==QCHR[L2]
									&&(GRID[QCHR[L1]][QLOC[L1]]-GRID[QCHR[L2]][QLOC[L2]])<=DQQ[QCHR[L1]]) T=0;
							}
						}
						if(T!=0) EpistasisIndicator_GROUP1(L1,L2);
					}
				}
			}

		}
}  

//************************************************************
//update g by e fixed effects INDICATORS

if(GBYE==1)
{ 
	for(L1=0;L1<NFIXCOVA;L1++)
	if(GBYE_FIX_INDEX[L1]==1)
		for(L2=0;L2<NQTL;L2++)
		if(GAMMA[L2]!=0)
		{	
			if(GROUP==0)
			{
				for(K=0;K<NC;K++)
				{
					if(GIBBS==1) GBYE_FIX_Indicator_GROUP0(L1,L2,K);
					if(GIBBS==0)
					{
						double R=RANDOM();
						if( (GBYE_FIX[L1][L2][K]==0&&R<=W_GBYE)||(GBYE_FIX[L1][L2][K]!=0&&R>W_GBYE) )
							GBYE_FIX_Indicator_GROUP0(L1,L2,K);
					}
				}
			}

			if(GROUP==1)
			{
				if(GIBBS==1) GBYE_FIX_Indicator_GROUP1(L1,L2);

				if(GIBBS==0)
				{
					double R=RANDOM();
					if( (GAMMA_GBYE[L1][L2]==0&&R<=W_GBYE)||(GAMMA_GBYE[L1][L2]!=0&&R>W_GBYE) )
						GBYE_FIX_Indicator_GROUP1(L1,L2);
				}
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
//CALCULATE THE NUMBER OF QTL

int QTL_INCLUDED=0;
for(L=0;L<NQTL;L++) QTL_INCLUDED=QTL_INCLUDED+(GAMMA[L]!=0);

//*************************************************************** 
//SAVE THE RESULT 

if(ITER*ITER1>=NBURNIN)
{

if((ITER!=0)&&(ITER%200==0)&(VERBOSE>0))  
{
	Rprintf("%d",ITER);
	Rprintf("\n");
}

// save to "iterdiag"
fprintf(File1,"\n");
fprintf(File1,"%d\t",ITER+1);
fprintf(File1,"%d\t",QTL_INCLUDED);
fprintf(File1,"%f\t",AMU);
fprintf(File1,"%f\t",VE);

// save to "covariates"
if(ENV_FACTOR==1)
{
	fprintf(File2,"\n");
	for(L=0;L<NFIXCOVA;L++) fprintf(File2,"%f\t",FIX[L]);
	for(L=0;L<NRANCOVA;L++) fprintf(File2,"%f\t",VRAN[L]);
} 
if(CATEGORY==3)
  for(L=2;L<=CN-2;L++) fprintf(File2,"%f\t",CUTPOINT[L]);
      
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
if(ENV_FACTOR==1) fprintf(File1,"%f\t",VP-VAR-VE);
fprintf(File1,"%f\t",VAR);
      
}

}     

//ITER end here

fclose(File1);  
fclose(File2); 
fclose(File3); 
fclose(File4);
fclose(File5);


/*
Rprintf("\n");

Rprintf("simulate %d MCMC steps,", NITER*NTHIN);
Rprintf(" record by %d", NTHIN);
Rprintf("\n");

if(EPISTASIS==0) Rprintf("This is a non-epistatic model");
if(EPISTASIS==1) Rprintf("This is an epistatic model");Rprintf("\n");


Rprintf("Prior number of main-effect QTL: %d", E_NQTL_MAIN);Rprintf("\n");

if(EPISTASIS==1) {Rprintf("Prior number of all QTL: %d", E_NQTL);Rprintf("\n");}

Rprintf("Maximum number of QTL: %d", NQTL);Rprintf("\n");

Rprintf("maximun number of QTL at each chromosome: ");
for(L=0;L<NLG;L++) Rprintf("%d ", CHR_NQTL[L]);Rprintf("\n");

Rprintf("Prior of main effect indicator: %f", W_MAIN);Rprintf("\n");

if(EPISTASIS==1) Rprintf("Prior of epistatic effect indicator: %f", W_EPISTASIS);Rprintf("\n");
 
Rprintf("\n"); */    
      
}

      
