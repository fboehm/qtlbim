

//***************************************************************************

#include "stdio.h"
//#include "bmq_main.h"
#include "bmq_mcmc.h"

//***************************************************************************

void R_AnalysisEngine(int *nind,int *nchr,int *ngen, int *nloci,double *loci,double *prob,
					  double *yvalue,int *traittype,int *ncategory,
					  int *iter,int *thin,int *burnin,
                      int *algorithm,int *genoupdate,
                      int *epis,int *emainqtl,int *eqtl,int *mnqtl,
					  double *interval,int *chrnqtl,
					  int *envi,int *qtl_envi,int *nrancov,int *nfixcov,int *intcov,double *rancoef,double *fixcoef,int *nran,
					  int *depen,double *prop,int *seed,
		      int *verbose)  
{
	int l,i, j, k;

// **************************************************************************
// parameter values passed from R function bmq.mcmc


    NS=*nind;
	NS1=*nind;
	NLG=*nchr;
	NG=*ngen;
	if(NG==3) CROSS=2;
	else CROSS=1;

	NGRID=(int *)S_alloc(NLG, sizeof(int));
    NGRID=nloci;
    TNGRID=0;
	for(i=0;i<NLG;i++) TNGRID=TNGRID+NGRID[i];
	CHL=NGRID[0];
    for(i=1;i<NLG;i++)
		if(NGRID[i]>CHL) CHL=NGRID[i];


	GRID=(double **)S_alloc(NLG, sizeof(double *));
	for(i=0;i<NLG;i++) GRID[i]=(double *)S_alloc(CHL, sizeof(double));
    int pos=0;
	for(i=0;i<NLG;i++)
		for(j=0;j<NGRID[i];j++) 
		{
			GRID[i][j]=loci[pos];
			pos=pos+1;
		}


    QPROB = (double ****)S_alloc(NLG, sizeof(double ***));
    for(i=0; i<NLG; i++) 
	{
		QPROB[i] = (double ***)S_alloc(NS, sizeof(double**));
        for(j=0; j<NS; j++) 
		{
			QPROB[i][j] = (double **)S_alloc(CHL, sizeof(double *));
			for(k=0; k<CHL; k++) 
			{
				QPROB[i][j][k] = (double *)S_alloc(NG, sizeof(double));
			}
		}
	}
	pos=0;
	for(l=0;l<NLG;l++)
		for(k=0;k<NG;k++)
			for(j=0;j<NGRID[l];j++)
				for(i=0;i<NS1;i++)
				{
					QPROB[l][i][j][k]=prob[pos];
					pos=pos+1;
				}



// test
/*
	FILE *File2;
    File2=fopen("C://Temp//test2.dat","a");
	
	for(l=0;l<NLG;l++)
	{
		fprintf(File2,"\n");
	    for(i=0;i<NS1;i++)
		{
			for(j=0;j<NGRID[l];j++)
			{
				for(k=0;k<NG;k++) fprintf(File2,"%f\t", QPROB[l][i][j][k]);
				fprintf(File2,"\n");
			}
		}
	}
    fclose(File2);   */



	Y=(double *)S_alloc(NS1, sizeof(double));
	Y=yvalue;


	CATEGORY=*traittype;
	CN=*ncategory;
	
	
	NITER=*iter;
	NTHIN=*thin;
	NBURNIN=*burnin;
	VERBOSE=*verbose;

	GIBBS=*algorithm;
	UPDATEGENO=*genoupdate;


	EPISTASIS=*epis;
	E_NQTL_MAIN=*emainqtl;
	E_NQTL=*eqtl;
	NQTL=*mnqtl;

	DQQ=(double *)S_alloc(NLG, sizeof(double));
	DQQ=interval;

	CHR_NQTL=(int *)S_alloc(NLG, sizeof(int));
	CHR_NQTL=chrnqtl;


	ENV_FACTOR=*envi; 
	GBYE=*qtl_envi;  
	
	NRANCOVA=*nrancov;  
	NFIXCOVA=*nfixcov;  

	GBYE_FIX_INDEX=(int *)S_alloc(NFIXCOVA, sizeof(int));
	GBYE_FIX_INDEX=intcov; 

	COEF_RAN = (double **)S_alloc(NS1, sizeof(double *));
	for(i=0; i<NS1; i++) COEF_RAN[i] = (double *)S_alloc(NRANCOVA, sizeof(double));
	for(i=0; i<NS1; i++) COEF_RAN[i]=rancoef+i*(*nrancov);

	NRAN=(int *)S_alloc(NRANCOVA,sizeof(int));
	NRAN=nran;   

	COEF_FIX = (double **)S_alloc(NS1, sizeof(double *));
	for(i=0; i<NS1; i++) COEF_FIX[i] = (double *)S_alloc(NFIXCOVA, sizeof(double));
	for(i=0; i<NS1; i++) COEF_FIX[i]=fixcoef+i*(*nfixcov);
 

	DEPENDENCE=*depen;
	C=prop;
	
	SEED=*seed;
     
//***************************************************************************
    
	GROUP=0;
	NC=(GROUP==0)*(NG-1)+(GROUP==1)*NG;

//**************************************************************************
// for binary and ordinal traits

	W=(int *)S_alloc(NS1, sizeof(int));
	CUTPOINT=(double *)S_alloc((CN+1), sizeof(double));

//**************************************************************************
// parameters used in prior specification

	//Calculate the priors of main and epistatic effect indicators
	int NC0=(GROUP==0)*NC+(GROUP==1)*1;
	W_MAIN=0;
	if(NQTL!=0) W_MAIN=1-pow(1-E_NQTL_MAIN*1.0/(NQTL*1.0),1.0/NC0);
	W_EPISTASIS=0;
	if(EPISTASIS==1&&NQTL!=0) W_EPISTASIS=1-pow( (1-E_NQTL*1.0/(NQTL*1.0))/pow(1-W_MAIN,NC0),1.0/(NC0*NC0*(NQTL-1.0)) );

	W_GBYE=0;
	if(NQTL!=0) W_GBYE=(1.0*E_NQTL_MAIN)/NQTL;

	VMAIN=(double *)S_alloc(NG, sizeof(double));

	VEPISTASIS = (double **)S_alloc(NG, sizeof(double *));
	for(i=0; i<NG; i++) VEPISTASIS[i] = (double *)S_alloc(NG, sizeof(double));

	V_GBYE_FIX = (double **)S_alloc(NFIXCOVA, sizeof(double *));
	for(i=0; i<NFIXCOVA; i++) V_GBYE_FIX[i] = (double *)S_alloc(NG, sizeof(double));

//*********************************************************************************
// genetic model parameters

	MAIN = (double **)S_alloc(NQTL, sizeof(double *));
	for(i=0; i<NQTL; i++) MAIN[i] = (double *)S_alloc(NG, sizeof(double));

	EPISTATIC = (double ****)S_alloc(NQTL, sizeof(double ***));
	for(i=0; i<NQTL; i++) 
	{
		EPISTATIC[i] = (double ***)S_alloc(NQTL, sizeof(double**));
		for(j=0; j<NQTL; j++) 
		{
			EPISTATIC[i][j] = (double **)S_alloc(NG, sizeof(double *));
			for(k=0; k<NG; k++) EPISTATIC[i][j][k] = (double *)S_alloc(NG, sizeof(double));
		}
    }

	GVALUE=(double *)S_alloc(NS1, sizeof(double));

	GENO = (int **)S_alloc(NS1, sizeof(int *));
	for(i=0; i<NS1; i++) GENO[i] = (int *)S_alloc(NQTL, sizeof(int));

	COEF = (double ***)S_alloc(NS1, sizeof(double **));
	for(i=0; i<NS1; i++) 
	{
		COEF[i] = (double **)S_alloc(NQTL, sizeof(double *));
		for(j=0;j<NQTL;j++) COEF[i][j]=(double *)S_alloc(NG, sizeof(double));
	}

//*********************************************************************************
// QTL positions, genetic effects indicators

	QLOC=(int *)S_alloc(NQTL, sizeof(int));
	QCHR=(int *)S_alloc(NQTL, sizeof(int));
	CHRQTL=(int *)S_alloc(NLG, sizeof(int));

	GAMMA=(int *)S_alloc(NQTL, sizeof(int));
	GAMMA_MAIN=(int *)S_alloc(NQTL, sizeof(int));

	GAMMA_EPISTASIS = (int **)S_alloc(NQTL, sizeof(int *));
	for(i=0; i<NQTL; i++) GAMMA_EPISTASIS[i] = (int *)S_alloc(NQTL, sizeof(int));


	for(i=0;i<NQTL;i++) 
	{
		QLOC[i]=0;
		QCHR[i]=0;
		GAMMA[i]=0;
		GAMMA_MAIN[i]=0;
		for(j=0;j<NQTL;j++) GAMMA_EPISTASIS[i][j]=0;
	} 	
	for(i=0;i<NLG;i++) CHRQTL[i]=0;  

//**********************************************************************************
// environmental covariate parameters

	FIX=(double *)S_alloc(NFIXCOVA, sizeof(double));

	RAN = (double **)S_alloc(NRANCOVA, sizeof(double *));
	for(i=0; i<NRANCOVA; i++) RAN[i] = (double *)S_alloc(NS1, sizeof(double));

	VRAN=(double *)S_alloc(NRANCOVA, sizeof(double));

	GBYE_FIX = (double ***)S_alloc(NFIXCOVA, sizeof(double **));
	for(i=0; i<NFIXCOVA; i++) 
	{
		GBYE_FIX[i] = (double **)S_alloc(NQTL, sizeof(double *));
		for(j=0;j<NQTL;j++) GBYE_FIX[i][j]=(double *)S_alloc(NG, sizeof(double));
	}

	GAMMA_GBYE = (double **)S_alloc(NFIXCOVA, sizeof(double *));
	for(i=0; i<NFIXCOVA; i++) GAMMA_GBYE[i] = (double *)S_alloc(NQTL, sizeof(double));

//*********************************************************************************

	
	PD1=(double *)S_alloc(NQTL, sizeof(double));
	PD2=(double *)S_alloc(NQTL, sizeof(double));

	X=(double *)S_alloc(NG, sizeof(double));

//**********************************************************************************

	// Assign the binary or ordinal phenotypes

	int I,J,K,NL;

	for(I=0;I<NS1;I++)
		if(CATEGORY!=1) W[I]=(int)Y[I];


	// REMOVE THE INDIVIDUALS WITH MISSING PHENOTYPIC AND COVARIATE VALUES

	int II=-1; 
	for(I=0;I<NS1;I++)
	{ 	
		int MISS=0;			
		for(J=0;J<NFIXCOVA;J++) MISS=MISS+(COEF_FIX[I][J]==999);
		for(J=0;J<NRANCOVA;J++) MISS=MISS+(COEF_RAN[I][J]==999);
	
		if(Y[I]!=999&&MISS==0)
		{
			II=II+1;
			Y[II]=Y[I], W[II]=W[I];
			for(NL=0;NL<NLG;NL++)
				for(J=0;J<NGRID[NL];J++)
					for(K=0;K<NG;K++) QPROB[NL][II][J][K]=QPROB[NL][I][J][K];
			
			for(J=0;J<NFIXCOVA;J++) COEF_FIX[II][J]=COEF_FIX[I][J];
			for(J=0;J<NRANCOVA;J++) COEF_RAN[II][J]=COEF_RAN[I][J];
		}
	}
	NS=II+1;      
	 

//*******************************************
// call mcmc algorithm

	bmqAnalysis();

}
