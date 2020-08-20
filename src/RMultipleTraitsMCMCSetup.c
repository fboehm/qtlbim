#include "GlobalVars.h"
#include "GlobalVars_MultipleTraits.h"
#include "RInterfaces.h"
#include "MultipleTraitsMCMC.h"

#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>

#include <stdio.h>
#include <time.h>
#include <math.h>


void RMultipleTraitsMCMCSetup(int *nind,int *nchr,int *ngen,int *npheno,int *nloci,
double *loci,double *prob,double *yvalue,int *multipletrait,int *traittype,int *ncategory,
int *iter,int *thin,int *burnin,
  int *algorithm,int *genoupdate,
                      int *epis,int *emainqtl,int *eqtl,int *mnqtl,
					  double *interval,int *chrnqtl,
	  int *envi,int *qtl_envi,int *nrancov,int *nfixcov,int *intcov,double *rancoef,
   double *fixcoef,int *nran, int *depen,double *prop,int *seed,int *verbose,
   int *diffloc, int *qtlloc)
              
{
	int l,i,j,k,PH;

// **************************************************************************
// parameter values passed from R function bmq.mcmc
  QTLLOC=*qtlloc;
  NS=*nind;
	NS1=*nind;
	NLG=*nchr;
	NG=*ngen;
	NPHENO=*npheno;
	MULTIPLE=*multipletrait;
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

  
  Y=(double **)S_alloc(NPHENO,sizeof(double *));
	for(i=0;i<NPHENO;i++) 	Y[i]=(double *)S_alloc(NS1, sizeof(double));
  pos=0;
	for(j=0;j<NPHENO;j++)
			{
				for(k=0;k<NS1;k++) 
        {
          Y[j][k]=yvalue[pos];
          pos++;
        }
			}   


	CATEGORY=*traittype;
	CN=*ncategory;
	
	
	NITER=*iter;
	NTHIN=*thin;
	NBURNIN=*burnin;
	VERBOSE=*verbose;
	DiffLocation=*diffloc;

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

	VMAIN=(double ***)S_alloc(NPHENO, sizeof(double **));
  for(i=0;i<NPHENO;i++){
  	VMAIN[i]=(double **)S_alloc(NQTL, sizeof(double *));
  	for(j=0;j<NQTL;j++)   	VMAIN[i][j]=(double *)S_alloc(NG, sizeof(double));
  }


	VEPISTASIS = (double *****)S_alloc(NPHENO, sizeof(double ****));
  	for(l=0;l<NPHENO;l++)
  	{
  	VEPISTASIS[l] = (double ****)S_alloc(NQTL, sizeof(double ***));
    	for(i=0; i<NQTL; i++)
    	{
    		VEPISTASIS[l][i] = (double ***)S_alloc(NQTL, sizeof(double**));
      		for(j=0; j<NQTL; j++)
      		{
      			VEPISTASIS[l][i][j] = (double **)S_alloc(NG, sizeof(double *));
      			for(k=0; k<NG; k++) VEPISTASIS[l][i][j][k] = (double *)S_alloc(NG, sizeof(double));
      		}
        }
    }
    
	V_GBYE_FIX = (double ****)S_alloc(NPHENO, sizeof(double ***));
	for(i=0;i<NPHENO;i++)
	{
  	V_GBYE_FIX[i] = (double ***)S_alloc(NFIXCOVA, sizeof(double **));
  	for(j=0; j<NFIXCOVA; j++) 
  	{
      V_GBYE_FIX[i][j] = (double **)S_alloc(NQTL, sizeof(double *));
      for(k=0;k<NQTL;k++) V_GBYE_FIX[i][j][k] = (double *)S_alloc(NG, sizeof(double));
    }
  }
//*********************************************************************************
// genetic model parameters

  AMU = (double *)S_alloc(NPHENO,sizeof(double));
  SIGMA = (double **)S_alloc(NPHENO,sizeof(double *));
  for(i=0;i<NPHENO;i++) SIGMA[i]=(double *)S_alloc(NPHENO,sizeof(double));

	MAIN = (double ***)S_alloc(NPHENO, sizeof(double **));
	for(i=0; i<NPHENO; i++)
  {
   MAIN[i] = (double **)S_alloc(NQTL, sizeof(double *));
	for(j=0; j<NQTL; j++) MAIN[i][j] = (double *)S_alloc(NG, sizeof(double)); 
  }                                            
  EPISTATIC = (double *****)S_alloc(NPHENO, sizeof(double ****));
  for(PH=0;PH<NPHENO;PH++)
  {  	EPISTATIC[PH] = (double ****)S_alloc(NQTL, sizeof(double ***));
  	for(i=0; i<NQTL; i++) 
  	{
		EPISTATIC[PH][i] = (double ***)S_alloc(NQTL, sizeof(double**));
	 	for(j=0; j<NQTL; j++) 
	   	{
			EPISTATIC[PH][i][j] = (double **)S_alloc(NG, sizeof(double *));
			for(k=0; k<NG; k++) 
          EPISTATIC[PH][i][j][k] = (double *)S_alloc(NG, sizeof(double));
    	}
    }
  }

	GVALUE=(double **)S_alloc(NPHENO, sizeof(double *));
	for(i=0;i<NPHENO;i++) GVALUE[i] = (double *)S_alloc(NS1,sizeof(double));

  GENO = (int ***)S_alloc(NPHENO, sizeof(int **));
  for(k=0;k<NPHENO;k++)
  {
  	GENO[k] = (int **)S_alloc(NS1, sizeof(int *));
  	for(i=0; i<NS1; i++) GENO[k][i] = (int *)S_alloc(NQTL, sizeof(int));
  }	

  COEF = (double ****)S_alloc(NPHENO, sizeof(double ***));
  for(k=0;k<NPHENO;k++)
  {
  COEF[k] = (double ***)S_alloc(NS1, sizeof(double **));
    for(i=0; i<NS1; i++) 
  	{
  		COEF[k][i] = (double **)S_alloc(NQTL, sizeof(double *));
  		for(j=0;j<NQTL;j++) COEF[k][i][j]=(double *)S_alloc(NG, sizeof(double));
  	}
  }
  
//******************************************************************************
// QTL positions, genetic effects indicators
// This is for Bayesian SUR which will be used later

	QLOC=(int **)S_alloc(NPHENO, sizeof(int *));
	for(i=0;i<NPHENO;i++) 	QLOC[i]=(int *)S_alloc(NQTL, sizeof(int));

	QCHR=(int **)S_alloc(NPHENO, sizeof(int *));
	for(i=0;i<NPHENO;i++) 	QCHR[i]=(int *)S_alloc(NQTL, sizeof(int));

	CHRQTL=(int *)S_alloc(NLG, sizeof(int));


	GAMMA=(int **)S_alloc(NPHENO, sizeof(int *));
	for(i=0;i<NPHENO;i++) 	GAMMA[i]=(int *)S_alloc(NQTL, sizeof(int));

	GAMMA_MAIN=(int **)S_alloc(NPHENO, sizeof(int *));
	for(i=0;i<NPHENO;i++) 	GAMMA_MAIN[i]=(int *)S_alloc(NQTL, sizeof(int));

	GAMMA_EPISTASIS = (int ***)S_alloc(NPHENO, sizeof(int **));
	for(i=0;i<NPHENO;i++)
  { 
  GAMMA_EPISTASIS[i] = (int **)S_alloc(NQTL, sizeof(int *));
	for(j=0; j<NQTL; j++) GAMMA_EPISTASIS[i][j] = (int *)S_alloc(NQTL, sizeof(int));
  }


//**********************************************************************************
// environmental covariate parameters
	FIX=(double **)S_alloc(NPHENO, sizeof(double *));
  for(i=0;i<NPHENO;i++)	FIX[i]=(double *)S_alloc(NFIXCOVA, sizeof(double));

	RAN = (double ***)S_alloc(NPHENO, sizeof(double **));
  for(i=0;i<NPHENO;i++)
  {
	RAN[i] = (double **)S_alloc(NRANCOVA, sizeof(double *));
	for(j=0; j<NRANCOVA; j++) RAN[i][j] = (double *)S_alloc(NS1, sizeof(double));
  }

  VRAN=(double **)S_alloc(NPHENO, sizeof(double *));
  for(i=0;i<NPHENO;i++)	VRAN[i]=(double *)S_alloc(NRANCOVA, sizeof(double));

	GBYE_FIX = (double ****)S_alloc(NPHENO, sizeof(double ***));
	for(k=0;k<NPHENO;k++)
	{
	GBYE_FIX[k] = (double ***)S_alloc(NFIXCOVA, sizeof(double **));
	for(i=0; i<NFIXCOVA; i++) 
	 {
		GBYE_FIX[k][i] = (double **)S_alloc(NQTL, sizeof(double *));
		for(j=0;j<NQTL;j++) GBYE_FIX[k][i][j]=(double *)S_alloc(NG, sizeof(double));
	 }
  }
	GAMMA_GBYE = (double ***)S_alloc(NPHENO, sizeof(double **));
  for(i=0;i<NPHENO;i++)
  {
	GAMMA_GBYE[i] = (double **)S_alloc(NFIXCOVA, sizeof(double *));
	for(j=0; j<NFIXCOVA; j++) GAMMA_GBYE[i][j] = (double *)S_alloc(NQTL, sizeof(double));
  }
//*********************************************************************************

	
	PD1=(double *)S_alloc(NQTL, sizeof(double));
	PD2=(double *)S_alloc(NQTL, sizeof(double));

	X=(double *)S_alloc(NG, sizeof(double));

//**********************************************************************************

	// Assign the binary or ordinal phenotypes

	int I,J,K,NL;

	for(I=0;I<NS1;I++)
		if(CATEGORY!=1) W[I]=(int)Y[0][I]; //Don't know if this will work W[I]=(int)Y[I];


	// REMOVE THE INDIVIDUALS WITH MISSING PHENOTYPIC AND COVARIATE VALUES
  
	int II=-1;
	for(I=0;I<NS1;I++)
	{ 	
		int MISS=0;			
    int YMISS=0; 
		for(J=0;J<NFIXCOVA;J++) MISS=MISS+(COEF_FIX[I][J]==999);
		for(J=0;J<NRANCOVA;J++) MISS=MISS+(COEF_RAN[I][J]==999);
  	for(K=0;K<NPHENO;K++) YMISS = YMISS + (Y[K][I]==999);

		if(YMISS==0&&MISS==0)
		{
			II=II+1;
			for(K=0;K<NPHENO;K++) Y[K][II]=Y[K][I];
      W[II]=W[I];
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
multipleTraitsMCMC();
 
}
