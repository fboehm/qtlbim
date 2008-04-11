// QTLBIM - QTL Bayesian Interval Mapping
// Functions relevant to single trait analyses
//********************************************************************

#include <stdio.h>
#include <math.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>

#include "GlobalVars.h"
#include "StatUtils.h"
#include "SingleTraitMCMCSamplingRoutines.h"




/// Calculate the likelihood for categorical trait
double Likelihood(double *p,double *G)
{
	double CC[CN+1], AL=0.0; int I,J;

    for(I=0;I<NS;I++)
	{
		for(J=0;J<=CN;J++) CC[J]=NormalFunction((p[J]-AMU-G[I])/sqrt(VE[I]));
		double T=0.0;
		for(J=1;J<=CN;J++) T=T+(W[I]==(J-1))*(CC[J]-CC[J-1]);
		AL=AL+log(T+1e-20);
	}

	return(AL);
}



// global variables - Cross,X.
///Transfering genotypes to regression coefficients
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
	int K;

	if(GROUP==1)
	{
		for(K=0;K<NG;K++) X[K]=QPROB[QCHR[L]][I][QL][K];
	}

	if(GROUP==0)
	{
		if(CROSS==2)                //Cockerham model
		{
			X[0]=QPROB[QCHR[L]][I][QL][2]-QPROB[QCHR[L]][I][QL][0],
			X[1]=0.5*(QPROB[QCHR[L]][I][QL][1]-QPROB[QCHR[L]][I][QL][0]-QPROB[QCHR[L]][I][QL][2]);
		}
		else X[0]=0.5*(QPROB[QCHR[L]][I][QL][1]-QPROB[QCHR[L]][I][QL][0]);
	}

	return;
}



//global variables - GAMMA *
///Deleting QTL with all 0 effects
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


//global variables - GAMMA *
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


/// Used for updating QTL genotypes and locations.
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



/// Update the overall mean
void Mean(double YBAR,double VP)     //prior for AMU is N(YBAR,VP)
{
	int I; double T1=0.0,T2=0.0,T3,U,U0;

    for(I=0;I<NS;I++) 
	{
		T1=T1+(Y[I]-GVALUE[I])/VE[I];
		T2=T2+1.0/VE[I];
	}
	T3=T2+1/VP;

	ANORMAL(&U,&U0);
	AMU=(T1+YBAR/VP)/T3+U/sqrt(T3);
	
	return;
}



/// Update the residual variance
void ResidualVariance()  //prior for VE is uniform
{
	int I; double T1=0.0,T2=0.0,U,U0;

    
		for(I=0;I<NS;I++)
		{
            T1=T1+pow(Y[I]-AMU-GVALUE[I],2);
            ANORMAL(&U,&U0);
            T2=T2+U*U;
		}

		for(I=0;I<NS;I++) VE[I]=T1/T2;  
        
		return;
}



/// Update the marginal effects
void MainEffect(int L,int K)
{
	int I; double G[NS1],T1,T2,T3,U,U0;

	T1=0,T2=0;
	for(I=0;I<NS;I++)
	{
		double Z=COEF[I][L][K];
		G[I]=GVALUE[I]-Z*MAIN[L][K];
		T1=T1+Z*(Y[I]-AMU-G[I])/VE[I];
		T2=T2+pow(Z,2)/VE[I];
	}
	double V=VMAIN[L][K];
	if(GROUP==1) V=VMAIN1[L];
	T3=1/V+T2;

	ANORMAL(&U,&U0);
	MAIN[L][K]=T1/T3+U/sqrt(T3);

	for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF[I][L][K]*MAIN[L][K];

	return;
}


/// Update the epistatic effects
void EpistaticEffect(int L1,int L2,int K1,int K2)
{
	int I; double G[NS1],T1,T2,T3,U,U0;

	T1=0,T2=0;
	for(I=0;I<NS;I++)
	{
		double Z=COEF[I][L1][K1]*COEF[I][L2][K2];
		G[I]=GVALUE[I]-Z*EPISTATIC[L1][L2][K1][K2];
		T1=T1+Z*(Y[I]-AMU-G[I])/VE[I];
		T2=T2+pow(Z,2)/VE[I];
	}
	double V=VEPISTASIS[L1][L2][K1][K2];
	if(GROUP==1) V=VEPISTASIS1[L1][L2];
	T3=1/V+T2;

	ANORMAL(&U,&U0);
	EPISTATIC[L1][L2][K1][K2]=T1/T3+U/sqrt(T3);

	for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF[I][L1][K1]*COEF[I][L2][K2]*EPISTATIC[L1][L2][K1][K2];

	return;
}


/// Update g by e fixed effects
void GBYE_FixedCovariate(int L1,int L2,int K)
{
	int I; double G[NS1],T1,T2,T3,U,U0;

	T1=0,T2=0;
	for(I=0;I<NS;I++)
	{
		double Z=COEF_FIX[I][L1]*COEF[I][L2][K];
		G[I]=GVALUE[I]-Z*GBYE_FIX[L1][L2][K];
		T1=T1+Z*(Y[I]-AMU-G[I])/VE[I];
		T2=T2+pow(Z,2)/VE[I];
	}
	double V=V_GBYE_FIX[L1][L2][K];
	if(GROUP==1) V=V_GBYE_FIX1[L1][L2];
	T3=1/V+T2;

	ANORMAL(&U,&U0);
	GBYE_FIX[L1][L2][K]=T1/T3+U/sqrt(T3);

	for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF_FIX[I][L1]*COEF[I][L2][K]*GBYE_FIX[L1][L2][K];

	return;
}


/// Updating genetic variances
void MainVariance(int L,int K,int NU,double TAU)    //prior for VMAIN is Inv-chisq(NU,(NU-2)/NU*TAU),E(VMAIN)=TAU
{
	int J; double T1,T2,U,U0;

	T1=pow(MAIN[L][K],2);
	T2=0;
	for(J=0;J<1+NU;J++)
	{
		ANORMAL(&U,&U0);
		T2=T2+U*U;
	}
	VMAIN[L][K]=(T1+NU*TAU)/T2;

	return;
}

void MainVariance1(int L,int NU,double TAU)    //prior for VMAIN is Inv-chisq(NU,(NU-2)/NU*TAU),E(VMAIN)=TAU
{
	int K,J; double T1,T2,U,U0;

	T1=0,T2=0;
	for(K=0;K<NC;K++) T1=T1+pow(MAIN[L][K],2);
	int N_MAIN=NC;
	for(J=0;J<N_MAIN+NU;J++)
	{
		ANORMAL(&U,&U0);
		T2=T2+U*U;
	}

	VMAIN1[L]=(T1+NU*TAU)/T2;

	return;
}


void EpistaticVariance(int L1,int L2,int K1,int K2,int NU,double TAU)    //prior for VEPISTASIS is Inv-chisq(NU,(NU-2)/NU*TAU),E(VEPISTASIS)=TAU
{
	int J;double T1,T2,U,U0;

	T1=pow(EPISTATIC[L1][L2][K1][K2],2);
	T2=0;
	for(J=0;J<1+NU;J++)
	{
		ANORMAL(&U,&U0);
		T2=T2+U*U;
	}
	VEPISTASIS[L1][L2][K1][K2]=(T1+NU*TAU)/T2;

	return;
}

void EpistaticVariance1(int L1,int L2,int NU,double TAU)    //prior for VEPISTASIS is Inv-chisq(NU,(NU-2)/NU*TAU),E(VEPISTASIS)=TAU
{
	int K1,K2,J;double T1,T2,U,U0;

	T1=0,T2=0;
	for(K1=0;K1<NC;K1++)
		for(K2=0;K2<NC;K2++) T1=T1+pow(EPISTATIC[L1][L2][K1][K2],2);

	int N_EPIS=NC*NC;
	for(J=0;J<N_EPIS+NU;J++)
	{
		ANORMAL(&U,&U0);
		T2=T2+U*U;
	}
	VEPISTASIS1[L1][L2]=(T1+NU*TAU)/T2;

	return;
}



/// Updating genetic variances
void GBYE_FixedCovariate_Variance(int L1,int L2,int K,int NU,double TAU)  //prior for V_GBYE is Inv-chisq(NU,(NU-2)/NU*TAU),E(V_GBYE)=TAU
{
	int J; double T1,T2,U,U0;

	T1=pow(GBYE_FIX[L1][L2][K],2);
	T2=0;
	for(J=0;J<1+NU;J++)
	{
		ANORMAL(&U,&U0);
		T2=T2+U*U;
	}
	V_GBYE_FIX[L1][L2][K]=(T1+NU*TAU)/T2;
	
	return;
}

void GBYE_FixedCovariate_Variance1(int L1,int L2,int NU,double TAU)  //prior for V_GBYE is Inv-chisq(NU,(NU-2)/NU*TAU),E(V_GBYE)=TAU
{
	int K,J; double T1,T2,U,U0;

	T1=0,T2=0;
	for(K=0;K<NC;K++) T1=T1+pow(GBYE_FIX[L1][L2][K],2);

	int N_GBYE=NC;
	for(J=0;J<N_GBYE+NU;J++)
	{
		ANORMAL(&U,&U0);
		T2=T2+U*U;
	}
	V_GBYE_FIX1[L1][L2]=(T1+NU*TAU)/T2;

	return;
}



/// Update nongenetic effects
void FixedCovariate(int L)
{
	int I; double G[NS1],Z,T1=0.0,T2=0.0,T3,U,U0;

    for(I=0;I<NS;I++)
    {
		Z=COEF_FIX[I][L];
		G[I]=GVALUE[I]-Z*FIX[L];
		T1=T1+Z*(Y[I]-AMU-G[I])/VE[I];
        T2=T2+Z*Z/VE[I];
    }
	T3=T2;	               //prior is uniform

	ANORMAL(&U,&U0);
	FIX[L]=T1/T3+U/sqrt(T3);

	for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF_FIX[I][L]*FIX[L];

	return;
}



///
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
				T1=T1+(Y[I]-AMU-G[I])/VE[I];
				T2=T2+1.0/VE[I];
			}
		T3=1/VRAN[L]+T2;

		ANORMAL(&U,&U0);
		RAN[L][K]=T1/T3+U/sqrt(T3);
	}
	for(I=0;I<NS;I++) GVALUE[I]=G[I]+RAN[L][(int)COEF_RAN[I][L]];

	return;
}



///
void RanVariance(int L)         //prior is uniform
{
	int K; double T1=0,T2=0,U,U0;

	for(K=0;K<NRAN[L];K++) T1=T1+pow(RAN[L][K],2);
	for(K=0;K<NRAN[L]-1;K++)
	{
        ANORMAL(&U,&U0);
        T2=T2+U*U;
	}

	VRAN[L]=T1/T2;      //see Gelman P301

	return;
}



/// Update QTL genotypes : generate IBD
void QTLgenotype(int L,int NL,int QL,int I)
{
	int K; double SUMM[NG];

	for(K=0;K<NG;K++) SUMM[K]=GenotypeSampling(I,L,K,QL);

    double SUM=0.0;
    for(K=0;K<NG;K++)
    {
	   SUMM[K]=exp(-0.5*pow(Y[I]-(AMU+SUMM[K]),2)/VE[I])*QPROB[NL][I][QL][K];
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



/// Update the positions of the QTLs
void QTLPOSITION(int L,int QLNEW)
{
	int I,GENO1[NS1],K; double G[NS1];

    double PROB0=0.0,PROB1=0.0,PD10=0.0,PD20=0.0;
    for(I=0;I<NS;I++)
	{
		G[I]=GVALUE[I];
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
			if(GENO[I][L]!=GENO1[I])
			{
				int II=GENO1[I];
				G[I]=GenotypeSampling(I,L,II,QLNEW);  //QLNEW is not used
			}
		}

		if(UPDATEGENO==0) G[I]=GenotypeSampling(I,L,0,QLNEW);  //0 is not used

		if(G[I]!=GVALUE[I])
		{
			PROB0=PROB0-0.5*pow(Y[I]-AMU-GVALUE[I],2)/VE[I];
			PROB1=PROB1-0.5*pow(Y[I]-AMU-G[I],2)/VE[I];
		}
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



/// Sampling a position
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



/// Update main effect indicators
void MainEffectIndicator_GROUP0(int L,int K)
{
	int I,K1;
	double G[NS1],BF_10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100),T1=0,T2=0,T3,F1=0,F2=0;


	for(I=0;I<NS;I++)
	{
		double Z=COEF[I][L][K];
		G[I]=GVALUE[I]-Z*MAIN[L][K];
		T1=T1+Z*(Y[I]-AMU-G[I])/VE[I];
		T2=T2+Z*Z/VE[I];
	}
	T3=1/VMAIN[L][K]+T2;

	for(I=0;I<NS;I++)
	{
		double Z=COEF[I][L][K];
		F1=F1+pow(Y[I]-AMU-G[I],2)/VE[I];
		F2=F2+pow(Y[I]-AMU-G[I]-Z*T1/T3,2)/VE[I];
	}

	BF_10=-0.5*log(VMAIN[L][K])-0.5*pow(T1/T3,2)/VMAIN[L][K]-0.5*log(T3);
	BF_10=-0.5*F2+0.5*F1+BF_10;

	if(MAIN[L][K]==0) GAMMA_1=BF_10;
	if(MAIN[L][K]!=0) GAMMA_0=-BF_10;

	double R=log(RANDOM());
	if(MAIN[L][K]==0&&R<GAMMA_1)
	{
		double U,U0;
		ANORMAL(&U,&U0);
		MAIN[L][K]=T1/T3+U/sqrt(T3);
		GAMMA_MAIN[L]=1;
		GAMMA[L]=1;
		for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF[I][L][K]*MAIN[L][K];
	}
	if(MAIN[L][K]!=0&&R<GAMMA_0)
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



///
void MainEffectIndicator_GROUP1(int L)
{
	int I,K;
	double G[NS1],BF_10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100),T1[NG],T2[NG],T3[NG],F1=0,F2=0;


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
			T1[K]=T1[K]+Z*(Y[I]-AMU-G[I])/VE[I];
			T2[K]=T2[K]+Z*Z/VE[I];
		}
		T3[K]=1/VMAIN1[L]+T2[K];
	}

	for(I=0;I<NS;I++)
	{
		F1=F1+pow(Y[I]-AMU-G[I],2)/VE[I];
		double F3=0;
		for(K=0;K<NC;K++) F3=F3+COEF[I][L][K]*T1[K]/T3[K];
		F2=F2+pow(Y[I]-AMU-G[I]-F3,2)/VE[I];
	}

	BF_10=0;
	for(K=0;K<NC;K++) BF_10=BF_10-0.5*log(VMAIN1[L])-0.5*pow(T1[K]/T3[K],2)/VMAIN1[L]-0.5*log(T3[K]);
	BF_10=-0.5*F2+0.5*F1+BF_10;

	if(GAMMA_MAIN[L]==0) GAMMA_1=BF_10;
	if(GAMMA_MAIN[L]!=0) GAMMA_0=-BF_10;

	double R=log(RANDOM());
	if(GAMMA_MAIN[L]==0&&R<GAMMA_1)
	{
		for(K=0;K<NC;K++)
		{
			double U,U0;
			ANORMAL(&U,&U0);
			MAIN[L][K]=T1[K]/T3[K]+U/sqrt(T3[K]);
		}
		GAMMA_MAIN[L]=1;
		GAMMA[L]=1;
		for(I=0;I<NS;I++)
		{
			GVALUE[I]=G[I];
			for(K=0;K<NC;K++) GVALUE[I]=GVALUE[I]+COEF[I][L][K]*MAIN[L][K];
		}
	}
	if(GAMMA_MAIN[L]!=0&&R<GAMMA_0)
	{
		for(K=0;K<NC;K++) MAIN[L][K]=0;
		GAMMA_MAIN[L]=0;
		for(I=0;I<NS;I++) GVALUE[I]=G[I];
		ZeroEffect1(L);
	}

	return;
}



/// Update epistatic effect indicators
void EpistasisIndicator_GROUP0(int L1,int L2,int K1,int K2)
{
	int I,K01,K02;
	double G[NS1],BF_10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100),T1=0,T2=0,T3,F1=0,F2=0;

	if(W_EPISTASIS!=0)
	{
	for(I=0;I<NS;I++)
	{
		double Z=COEF[I][L1][K1]*COEF[I][L2][K2];
		G[I]=GVALUE[I]-Z*EPISTATIC[L1][L2][K1][K2];
		T1=T1+Z*(Y[I]-AMU-G[I])/VE[I];
		T2=T2+Z*Z/VE[I];
	}
	T3=1/VEPISTASIS[L1][L2][K1][K2]+T2;

	for(I=0;I<NS;I++)
	{
		double Z=COEF[I][L1][K1]*COEF[I][L2][K2];
		F1=F1+pow(Y[I]-AMU-G[I],2)/VE[I];
		F2=F2+pow(Y[I]-AMU-G[I]-Z*T1/T3,2)/VE[I];
	}

	BF_10=-0.5*log(VEPISTASIS[L1][L2][K1][K2])-0.5*pow(T1/T3,2)/VEPISTASIS[L1][L2][K1][K2]-0.5*log(T3);
	BF_10=-0.5*F2+0.5*F1+BF_10;

	if(EPISTATIC[L1][L2][K1][K2]==0) GAMMA_1=BF_10;
	if(EPISTATIC[L1][L2][K1][K2]!=0) GAMMA_0=-BF_10;

	double R=log(RANDOM());
	if(EPISTATIC[L1][L2][K1][K2]==0&&R<GAMMA_1)
	{
		double U,U0;
		ANORMAL(&U,&U0);
		EPISTATIC[L1][L2][K1][K2]=T1/T3+U/sqrt(T3);
		GAMMA_EPISTASIS[L1][L2]=1;
		GAMMA[L1]=1,GAMMA[L2]=1;
		for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF[I][L1][K1]*COEF[I][L2][K2]*EPISTATIC[L1][L2][K1][K2];
	}
	if(EPISTATIC[L1][L2][K1][K2]!=0&&R<GAMMA_0)
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



///
void EpistasisIndicator_GROUP1(int L1,int L2)
{
	int I,K1,K2;
	double G[NS1],BF_10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100),T1[NG][NG],T2[NG][NG],T3[NG][NG],F1=0,F2=0;

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
				T1[K1][K2]=T1[K1][K2]+Z*(Y[I]-AMU-G[I])/VE[I];
				T2[K1][K2]=T2[K1][K2]+Z*Z/VE[I];
			}
			T3[K1][K2]=1/VEPISTASIS1[L1][L2]+T2[K1][K2];
		}

	for(I=0;I<NS;I++)
	{
		F1=F1+pow(Y[I]-AMU-G[I],2)/VE[I];
		double F3=0;
		for(K1=0;K1<NC;K1++)
			for(K2=0;K2<NC;K2++) F3=F3+COEF[I][L1][K1]*COEF[I][L2][K2]*T1[K1][K2]/T3[K1][K2];
		F2=F2+pow(Y[I]-AMU-G[I]-F3,2)/VE[I];
	}

	BF_10=0;
	for(K1=0;K1<NC;K1++)
		for(K2=0;K2<NC;K2++)
			BF_10=BF_10-0.5*log(VEPISTASIS1[L1][L2])-0.5*pow(T1[K1][K2]/T3[K1][K2],2)/VEPISTASIS1[L1][L2]-0.5*log(T3[K1][K2]);
	BF_10=-0.5*F2+0.5*F1+BF_10;

	if(GAMMA_EPISTASIS[L1][L2]==0) GAMMA_1=BF_10;
	if(GAMMA_EPISTASIS[L1][L2]!=0) GAMMA_0=-BF_10;

	double R=log(RANDOM());
	if(GAMMA_EPISTASIS[L1][L2]==0&&R<GAMMA_1)
	{
		for(K1=0;K1<NC;K1++)
			for(K2=0;K2<NC;K2++)
			{
				double U,U0;
				ANORMAL(&U,&U0);
				EPISTATIC[L1][L2][K1][K2]=T1[K1][K2]/T3[K1][K2]+U/sqrt(T3[K1][K2]);
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
	if(GAMMA_EPISTASIS[L1][L2]!=0&&R<GAMMA_0)
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



/// Update g by e fixed effect indicators
void GBYE_FIX_Indicator_GROUP0(int L1,int L2,int K)  //use the constraint model
{
	int I,K1;
	double G[NS1],BF_10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100),T1=0,T2=0,T3,F1=0,F2=0;


	for(I=0;I<NS;I++)
	{
		double Z=COEF_FIX[I][L1]*COEF[I][L2][K];
		G[I]=GVALUE[I]-Z*GBYE_FIX[L1][L2][K];
		T1=T1+Z*(Y[I]-AMU-G[I])/VE[I];
		T2=T2+Z*Z/VE[I];
	}
	T3=1/V_GBYE_FIX[L1][L2][K]+T2;

	for(I=0;I<NS;I++)
	{
		double Z=COEF_FIX[I][L1]*COEF[I][L2][K];
		F1=F1+pow(Y[I]-AMU-G[I],2)/VE[I];
		F2=F2+pow(Y[I]-AMU-G[I]-Z*T1/T3,2)/VE[I];
	}

	BF_10=-0.5*log(V_GBYE_FIX[L1][L2][K])-0.5*pow(T1/T3,2)/V_GBYE_FIX[L1][L2][K]-0.5*log(T3);
	BF_10=-0.5*F2+0.5*F1+BF_10;

	if(GBYE_FIX[L1][L2][K]==0) GAMMA_1=BF_10;
	if(GBYE_FIX[L1][L2][K]!=0) GAMMA_0=-BF_10;

	double R=log(RANDOM());
	if(GBYE_FIX[L1][L2][K]==0&&R<GAMMA_1)
	{
		double U,U0;
		ANORMAL(&U,&U0);
		GBYE_FIX[L1][L2][K]=T1/T3+U/sqrt(T3);
		GAMMA_GBYE[L1][L2]=1;
		GAMMA[L2]=1;
		for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF_FIX[I][L1]*COEF[I][L2][K]*GBYE_FIX[L1][L2][K];
	}
	if(GBYE_FIX[L1][L2][K]!=0&&R<GAMMA_0)
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



///
void GBYE_FIX_Indicator_GROUP1(int L1,int L2)
{
	int I,K;
	double G[NS1],BF_10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100),T1[NG],T2[NG],T3[NG],F1=0,F2=0;

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
			T1[K]=T1[K]+Z*(Y[I]-AMU-G[I])/VE[I];
			T2[K]=T2[K]+Z*Z/VE[I];
		}
		T3[K]=1/V_GBYE_FIX1[L1][L2]+T2[K];
	}

	for(I=0;I<NS;I++)
	{
		F1=F1+pow(Y[I]-AMU-G[I],2)/VE[I];
		double F3=0;
        for(K=0;K<NC;K++) F3=F3+COEF_FIX[I][L1]*COEF[I][L2][K]*T1[K]/T3[K];
		F2=F2+pow(Y[I]-AMU-G[I]-F3,2)/VE[I];
	}

	BF_10=0;
	for(K=0;K<NC;K++) BF_10=BF_10-0.5*log(V_GBYE_FIX1[L1][L2])-0.5*pow(T1[K]/T3[K],2)/V_GBYE_FIX1[L1][L2]-0.5*log(T3[K]);
	BF_10=-0.5*F2+0.5*F1+BF_10;

	if(GAMMA_GBYE[L1][L2]==0) GAMMA_1=BF_10;
	if(GAMMA_GBYE[L1][L2]!=0) GAMMA_0=-BF_10;

	double R=log(RANDOM());
	if(GAMMA_GBYE[L1][L2]==0&&R<GAMMA_1)
	{
		for(K=0;K<NC;K++)
		{
			double U,U0;
			ANORMAL(&U,&U0);
			GBYE_FIX[L1][L2][K]=T1[K]/T3[K]+U/sqrt(T3[K]);
		}
		GAMMA_GBYE[L1][L2]=1;
		GAMMA[L2]=1;
		for(I=0;I<NS;I++)
		{
			GVALUE[I]=G[I];
			for(K=0;K<NC;K++) GVALUE[I]=GVALUE[I]+COEF_FIX[I][L1]*COEF[I][L2][K]*GBYE_FIX[L1][L2][K];
		}
	}
	if(GAMMA_GBYE[L1][L2]!=0&&R<GAMMA_0)
	{
		for(K=0;K<NC;K++) GBYE_FIX[L1][L2][K]=0;
		GAMMA_GBYE[L1][L2]=0;
		for(I=0;I<NS;I++) GVALUE[I]=G[I];
		ZeroEffect1(L2);
	}

	return;
}




void MainEffect1(int L)
{
	int I,K; double G[NS1],T1,T2,T3,U,U0;

	for(K=0;K<NC;K++)
	{
		T1=0,T2=0;
		for(I=0;I<NS;I++)
		{
			double Z=COEF[I][L][K];
			G[I]=GVALUE[I]-Z*MAIN[L][K];
			T1=T1+Z*(Y[I]-AMU-G[I])/VE[I];
			T2=T2+pow(Z,2)/VE[I];
		}
		double V=VMAIN1[L];
		T3=1/V+T2;

		ANORMAL(&U,&U0);
		MAIN[L][K]=T1/T3+U/sqrt(T3);

		for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF[I][L][K]*MAIN[L][K];
	}

	return;
}


/// Update epistatic effects
void EpistaticEffect1(int L1,int L2)
{
	int I,K1,K2; double G[NS1],T1,T2,T3,U,U0;

	for(K1=0;K1<NC;K1++)
		for(K2=0;K2<NC;K2++)
		{
			T1=0,T2=0;
			for(I=0;I<NS;I++)
			{
				double Z=COEF[I][L1][K1]*COEF[I][L2][K2];
				G[I]=GVALUE[I]-Z*EPISTATIC[L1][L2][K1][K2];
				T1=T1+Z*(Y[I]-AMU-G[I])/VE[I];
				T2=T2+pow(Z,2)/VE[I];
			}
			double V=VEPISTASIS1[L1][L2];
			T3=1/V+T2;

			ANORMAL(&U,&U0);
			EPISTATIC[L1][L2][K1][K2]=T1/T3+U/sqrt(T3);

			for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF[I][L1][K1]*COEF[I][L2][K2]*EPISTATIC[L1][L2][K1][K2];
		}

	return;
}

/// Update g by e fixed effects
void GBYE_FixedCovariate1(int L1,int L2)
{
	int I,K; double G[NS1],T1,T2,T3,U,U0;

	for(K=0;K<NC;K++)
	{
		T1=0,T2=0;
		for(I=0;I<NS;I++)
		{
			double Z=COEF_FIX[I][L1]*COEF[I][L2][K];
			G[I]=GVALUE[I]-Z*GBYE_FIX[L1][L2][K];
			T1=T1+Z*(Y[I]-AMU-G[I])/VE[I];
			T2=T2+pow(Z,2)/VE[I];
		}
		double V=V_GBYE_FIX1[L1][L2];
		T3=1/V+T2;

		ANORMAL(&U,&U0);
		GBYE_FIX[L1][L2][K]=T1/T3+U/sqrt(T3);

		for(I=0;I<NS;I++) GVALUE[I]=G[I]+COEF_FIX[I][L1]*COEF[I][L2][K]*GBYE_FIX[L1][L2][K];
	}

	return;
}
