#include "GlobalVars.h"
#include "GlobalVars_MultipleTraits.h"
#include "StatUtils.h"
#include "MatrixUtils.h"
#include "MultipleTraitsMCMC.h"

#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>

#include <math.h>
#include <time.h>

// global variables - Cross,X.
///Transfering genotypes to regression coefficients
/*
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
} */


void ResidualVariance_MultipleTraits()
{
 double **newsigma,**G,**chol,*mu,*sample;
 double det,**temp,**omega;
 int I,J,i,j,k;

 newsigma = malloc(NPHENO*sizeof(double *));
 for(i=0;i<NPHENO;i++) newsigma[i]= malloc(NPHENO*sizeof(double));

 omega = malloc(NPHENO*sizeof(double *));
 for(i=0;i<NPHENO;i++) omega[i]= malloc(NPHENO*sizeof(double));

 chol = malloc(NPHENO*sizeof(double *));
 for(i=0;i<NPHENO;i++) chol[i]= malloc(NPHENO*sizeof(double));

 temp = malloc(NPHENO*sizeof(double *));
 for(i=0;i<NPHENO;i++) temp[i]= malloc(NPHENO*sizeof(double));

 G = malloc(NPHENO*sizeof(double *));
 for(i=0;i<NPHENO;i++) G[i]= malloc(NS*sizeof(double));
 
 mu = malloc(NPHENO*sizeof(double));

 sample = malloc(NPHENO*sizeof(double));
   
 for(i=0;i<NPHENO;i++)
 {
  mu[i]=0;
  sample[i]=0;
  for(j=0;j<NPHENO;j++)
    {
     newsigma[i][j]=0;
     temp[i][j]=0;
     omega[i][j]=0;
     chol[i][j]=0;     
    }
    
      for(k=0;k<NS;k++) G[i][k]=0;  

 }

  for(I=0;I<NS;I++)
     for(J=0;J<NPHENO;J++) 
          G[J][I] = Y[J][I] - AMU[J] - GVALUE[J][I];
 
 
 XprimeX(G,NPHENO,NS,omega);

     det = Determinant(omega,NPHENO);
     if(det < 1e-15 ){     
     Rprintf("\n\n Covariance matrix approaching near singularity.. quitting!\n");
     return;
     /* exit(1); */
     }   
     else {
     INVERSE(omega,NPHENO,temp);
     Cholesky(temp,NPHENO,chol);
     }
     

for(i=0;i<NPHENO;i++)
    for(j=0;j<NPHENO;j++)   {
      temp[i][j]=0;
      newsigma[i][j]=0;
      
}     
  for(i=0;i<NS;i++)
  {
     Multivariate_RNORM(chol,mu,NPHENO,sample);  
      for(j=0;j<NPHENO;j++)
        for(k=0;k<NPHENO;k++)
          temp[j][k] = sample[j]*sample[k];
     MatrixCopy(temp,newsigma,NPHENO,NPHENO,1);
     
   } 

//   if(GIBBS==1)   
MatrixCopy(newsigma,SIGMA,NPHENO,NPHENO,0);        
/*
if(GIBBS==0)
{

   det=Determinant(newsigma,NPHENO);
   if(det < 1e-15)  {
    Rprintf("\n Determinant of newly sampled Sigma near 0 det= %f\n",det);
    exit(1); 
   }
   
   double det1,like_old,like_new,alpha,R;
   XprimeY(omega,newsigma,NPHENO,NPHENO,NPHENO,temp);
   like_new = pow(det,NS/2)*exp(-0.5*Trace(temp,NPHENO));
   det1 = Determinant(SIGMA,NPHENO);
   XprimeY(omega,SIGMA,NPHENO,NPHENO,NPHENO,temp);
   like_old = pow(det1,NS/2)*exp(-0.5*Trace(temp,NPHENO));
   alpha = like_new/like_old;

   int accept=0;
   if(alpha < 1)
   {
    R=RANDOM();
    if(R < alpha)
    { MatrixCopy(newsigma,SIGMA,NPHENO,NPHENO,0);        
     accept = 1;
    }
   }
   else {
   MatrixCopy(newsigma,SIGMA,NPHENO,NPHENO,0);        
   accept=1;
   }
 }*/     
    for(i=0;i<NPHENO;i++)
    {
     free(newsigma[i]);
     free(omega[i]);
     free(temp[i]);
     free(chol[i]);
     free(G[i]);
    }
    free(newsigma);
    free(omega);
    free(temp);
    free(chol);
    free(mu);
    free(sample);
    free(G);
 
}


void MT_Coefficient0(int I,int L,int PH,int QL)     //UPDATEGENO=0
{ 
	if(CROSS==2)                //Cockerham model
	{
		X[0]=QPROB[QCHR[PH][L]][I][QL][2]-QPROB[QCHR[PH][L]][I][QL][0],
		X[1]=0.5*(QPROB[QCHR[PH][L]][I][QL][1]-QPROB[QCHR[PH][L]][I][QL][0]-QPROB[QCHR[PH][L]][I][QL][2]);
	}
	else X[0]=0.5*(QPROB[QCHR[PH][L]][I][QL][1]-QPROB[QCHR[PH][L]][I][QL][0]);

	return;
}

//******************************************************************************************
//used in updating QTL genotypes and locations.

double MT_GenotypeSampling(int I,int L,int PH,int II,int QL)   //if UPDATEGENO==1(0), dosenot need QL(II) 
{
    int L1,K,K1,K2; double G;
	
	if(UPDATEGENO==1) Coefficient(II);
	if(UPDATEGENO==0) MT_Coefficient0(I,L,PH,QL);
	G=GVALUE[PH][I];
	for(K=0;K<NC;K++) G=G-COEF[PH][I][L][K]*MAIN[PH][L][K]+X[K]*MAIN[PH][L][K];
			
	if(EPISTASIS==1)
		for(L1=0;L1<NQTL;L1++)
		{
			if(L1<L && GAMMA[PH][L1]!=0 && GAMMA_EPISTASIS[PH][L1][L]!=0)
				for(K1=0;K1<NC;K1++)
					for(K2=0;K2<NC;K2++) 
		G=G-COEF[PH][I][L1][K1]*COEF[PH][I][L][K2]*EPISTATIC[PH][L1][L][K1][K2]+COEF[PH][I][L1][K1]*X[K2]*EPISTATIC[PH][L1][L][K1][K2];

            if(L1>L && GAMMA[PH][L1]!=0 && GAMMA_EPISTASIS[PH][L][L1]!=0) 
				for(K1=0;K1<NC;K1++)
					for(K2=0;K2<NC;K2++) 
						G=G-COEF[PH][I][L][K1]*COEF[PH][I][L1][K2]*EPISTATIC[PH][L][L1][K1][K2]+X[K1]*COEF[PH][I][L1][K2]*EPISTATIC[PH][L][L1][K1][K2];
		}

	if(GBYE==1)
		for(L1=0;L1<NFIXCOVA;L1++) 
			if(GBYE_FIX_INDEX[L1]==1&&GAMMA_GBYE[PH][L1][L]!=0)
				for(K=0;K<NC;K++) 
					G=G-COEF_FIX[I][L1]*COEF[PH][I][L][K]*GBYE_FIX[PH][L1][L][K]
                                    +COEF_FIX[I][L1]*X[K]*GBYE_FIX[PH][L1][L][K];
		
	return(G);
} 
   
//*************************************************************************************
// update the overall mean

void MT_Mean(double YBAR,double VP,int trait)  //prior for AMU is N(YBAR,VP)
{
	int I,PH; double T=0.0,U,U0,G[NPHENO][NS1]; 
 if(MULTIPLE==0)
 {	
    for(I=0;I<NS;I++) T=T+(Y[0][I]-GVALUE[0][I]);
	
	ANORMAL(&U,&U0); 
	AMU[trait]=(T+YBAR*VE/VP)/(NS+VE/VP)+U*sqrt(VE/(NS+VE/VP));  
  }
 else if(MULTIPLE>=1)
 {
    for(I=0;I<NS;I++) 
    {  for(PH=0;PH<NPHENO;PH++)
      {
       if(PH==trait) G[PH][I]=GVALUE[PH][I];
       else G[PH][I]=GVALUE[PH][I] + AMU[PH];         
  
       T=T+(Y[PH][I]-G[PH][I])*SIGMA[trait][PH];
        }   
     }
  ANORMAL(&U,&U0); 
AMU[trait]=(T+YBAR/VP)/(NS*SIGMA[trait][trait]+1/VP)+U/sqrt(NS*SIGMA[trait][trait]+1/VP);
  
     
 } 
	return;
}
//*******************************************************
//update marginal effects
  
void MT_MainEffect(int L,int trait,int forindicators,int effect)
{
	int I,K,ph; double G[NPHENO][NS1],T1,T2,T3,U,U0;
	for(K=0;K<NC;K++) 
	{ if(forindicators==1) K=effect;
		if(MAIN[trait][L][K]!=0 || forindicators==1)
		{
			T1=0,T2=0;
			for(I=0;I<NS;I++)                     
			{
				double Z=COEF[trait][I][L][K];
				if(MULTIPLE==0)
				{
        G[trait][I]=GVALUE[trait][I]-Z*MAIN[trait][L][K];
				T1=T1+Z*(Y[trait][I]-AMU[trait]-G[trait][I]);
				T2=T2+pow(Z,2);
				}
				if(MULTIPLE>=1)
				{ 
				  for(ph=0;ph<NPHENO;ph++)
				  {
            if(ph==trait) G[ph][I]=GVALUE[ph][I]-Z*MAIN[ph][L][K];
            else G[ph][I]=GVALUE[ph][I];
             
	   		  	T1=T1+Z*(Y[ph][I]-AMU[ph]-G[ph][I])*SIGMA[trait][ph];
		  		  }
		  		  T2=T2+Z*Z;            								
				}
			}	
       	
  		ANORMAL(&U,&U0);
      if(MULTIPLE==0)	{      
     		T3=1/VMAIN[trait][L][K]+T2/VE;
	   		MAIN[trait][L][K]=T1/(VE*T3)+U/sqrt(T3);
			}
      else if(MULTIPLE>=1) {
        T3=T2*SIGMA[trait][trait] + 1/VMAIN[trait][L][K];
	   		MAIN[trait][L][K]=T1/T3+U/sqrt(T3);
      }
for(I=0;I<NS;I++) GVALUE[trait][I]=G[trait][I]+COEF[trait][I][L][K]*MAIN[trait][L][K];
		}
		  if(forindicators==1) K=NC-1;
		  
	
	}
	return;
}


//************************************************************
//update epistatic effects

void MT_EpistaticEffect(int L1,int L2,int trait)
{
	int I,K1,K2,ph; double G[NPHENO][NS1],T1,T2,T3,U,U0;

	for(K1=0;K1<NC;K1++)
		for(K2=0;K2<NC;K2++)
		{
			if(EPISTATIC[trait][L1][L2][K1][K2]!=0)
			{
			if(MULTIPLE==0)
			  {
      	T1=0,T2=0;
				for(I=0;I<NS;I++)                     
				{
					double Z=COEF[0][I][L1][K1]*COEF[0][I][L2][K2];
					G[trait][I]=GVALUE[trait][I]-Z*EPISTATIC[0][L1][L2][K1][K2];
					T1=T1+Z*(Y[0][I]-AMU[0]-G[0][I]);
					T2=T2+pow(Z,2);
				}	 	
				T3=1/VEPISTASIS[0][L1][L2][K1][K2]+T2/VE;

				ANORMAL(&U,&U0);
				EPISTATIC[0][L1][L2][K1][K2]=T1/(VE*T3)+U/sqrt(T3);
 
				for(I=0;I<NS;I++) 
          GVALUE[0][I]=G[0][I]+COEF[0][I][L1][K1]*COEF[0][I][L2][K2]*EPISTATIC[0][L1][L2][K1][K2];
				}
			if(MULTIPLE>=1)
			  {
      	T1=0,T2=0;
				for(I=0;I<NS;I++)                     
				{
					double Z=COEF[trait][I][L1][K1]*COEF[trait][I][L2][K2];
				  for(ph=0;ph<NPHENO;ph++)
				  {
        	if(ph==trait)  G[ph][I]=GVALUE[ph][I]-Z*EPISTATIC[ph][L1][L2][K1][K2];
          else G[ph][I]=GVALUE[ph][I];
          
					T1=T1+Z*(Y[ph][I]-AMU[ph]-G[ph][I])*SIGMA[trait][ph];
					}
					T2=T2+Z*Z;

				}	 	
        T3=T2*SIGMA[trait][trait] + 1.0/VEPISTASIS[trait][L1][L2][K1][K2];
	   		
   			ANORMAL(&U,&U0);
				EPISTATIC[trait][L1][L2][K1][K2]=T1/T3+U/sqrt(T3);
 
for(I=0;I<NS;I++) 
GVALUE[trait][I]=G[trait][I]+COEF[trait][I][L1][K1]*COEF[trait][I][L2][K2]*EPISTATIC[trait][L1][L2][K1][K2];
				}
			}
		}

	return;
} 

//************************************************************
//update g by e fixed effects

void MT_GBYE_FixedCovariate(int L1,int L2,int trait)
{
	int I,K,ph; double G[NPHENO][NS1],T1,T2,T3,U,U0;

	for(K=0;K<NC;K++)
	{
		if(GBYE_FIX[trait][L1][L2][K]!=0)
		{
     if(MULTIPLE==0)
     {
			T1=0,T2=0;
			for(I=0;I<NS;I++)                     
			{
				double Z=COEF_FIX[I][L1]*COEF[0][I][L2][K];
				G[0][I]=GVALUE[0][I]-Z*GBYE_FIX[0][L1][L2][K];
				T1=T1+Z*(Y[0][I]-AMU[0]-G[0][I]);
				T2=T2+pow(Z,2);
			}	 	
		 	T3=1/V_GBYE_FIX[0][L1][L2][K]+T2/VE;

			ANORMAL(&U,&U0);
			GBYE_FIX[0][L1][L2][K]=T1/(VE*T3)+U/sqrt(T3);

			for(I=0;I<NS;I++) GVALUE[0][I]=G[0][I]+COEF_FIX[I][L1]*COEF[0][I][L2][K]*GBYE_FIX[0][L1][L2][K];
			}
			if(MULTIPLE>=1)
			{
			T1=0,T2=0;
      for(I=0;I<NS;I++)
      {
			double Z=COEF_FIX[I][L1]*COEF[trait][I][L2][K];
		  for(ph=0;ph<NPHENO;ph++)
				  {
            if(ph==trait) G[ph][I]=GVALUE[ph][I]-Z*GBYE_FIX[ph][L1][L2][K];
            else G[ph][I]=GVALUE[ph][I];

	   		  	T1=T1+Z*(Y[ph][I]-AMU[ph]-G[ph][I])*SIGMA[trait][ph];
	  		  }
		  		  T2=T2+Z*Z;
       }
        T3=T2*SIGMA[trait][trait] + 1.0/V_GBYE_FIX[trait][L1][L2][K];
  			ANORMAL(&U,&U0);
	   		GBYE_FIX[trait][L1][L2][K]=T1/T3+U/sqrt(T3);

			for(I=0;I<NS;I++) 
      GVALUE[trait][I]=G[trait][I]+COEF_FIX[I][L1]*COEF[trait][I][L2][K]*GBYE_FIX[trait][L1][L2][K];

			
			}
			
		}
	}

	return;
}

//***********************************************************
//updating genetic variances


void MT_MainVariance(int L,int K,int NU,double TAU,int trait) //prior for VMAIN is Inv-chisq(NU,(NU-2)/NU*TAU),E(VMAIN)=TAU
{

	int J; double T1,T2,U,U0;

	T1=pow(MAIN[trait][L][K],2);
	T2=0;
	for(J=0;J<1+NU;J++)
	{
		ANORMAL(&U,&U0);
		T2=T2+U*U;
	}
	VMAIN[trait][L][K]=(T1+NU*TAU)/T2;

	return;

}

/*  
void MT_EpistaticVariance(double VP,int trait)    //prior for VEPISTASIS is Inv-chisq(NU,(NU-2)/NU*TAU),E(VEPISTASIS)=TAU
{
	int L1,L2,K1,K2,J;double T1,T2,U,U0;

	int NU=6; double TAU=0.1;

	for(K1=0;K1<NC;K1++)
		for(K2=0;K2<NC;K2++)
		{
			int N_EPIS=0;
			for(L1=0;L1<NQTL;L1++)
				for(L2=0;L2<NQTL;L2++)
					if(EPISTATIC[trait][L1][L2][K1][K2]!=0) N_EPIS=N_EPIS+1; 
			
			T1=0,T2=0;
			for(L1=0;L1<NQTL;L1++)
				for(L2=L1+1;L2<NQTL;L2++)
					if(EPISTATIC[trait][L1][L2][K1][K2]!=0) T1=T1+pow(EPISTATIC[trait][L1][L2][K1][K2],2);

			for(J=0;J<N_EPIS+NU;J++)
			{
				ANORMAL(&U,&U0);
				T2=T2+U*U; 
			} 
     	
			T1=T1/(VEPISTASIS[trait][K1][K2]*VP);T2=T2/(VEPISTASIS[trait][K1][K2]*VP);
			VEPISTASIS[trait][K1][K2]=(T1+NU*(NU-2.0)/NU*TAU)/T2;
		}

	return;
} */ 

void MT_EpistaticVariance(int L1,int L2,int K1,int K2,int NU,double TAU,int trait)    //prior for VEPISTASIS is Inv-chisq(NU,(NU-2)/NU*TAU),E(VEPISTASIS)=TAU
{
	int J;double T1,T2,U,U0;

	T1=pow(EPISTATIC[trait][L1][L2][K1][K2],2);
	T2=0;
	for(J=0;J<1+NU;J++)
	{
		ANORMAL(&U,&U0);
		T2=T2+U*U;
	}
	VEPISTASIS[trait][L1][L2][K1][K2]=(T1+NU*TAU)/T2;

	return;
}

//updating genetic variances
/*
void MT_GBYE_FixedCovariate_Variance(int L1,double VP,int trait)  //prior for V_GBYE is Inv-chisq(NU,(NU-2)/NU*TAU),E(V_GBYE)=TAU
{
	int L2,K,J; double T1,T2,U,U0;

	int NU=6; double TAU=0.1;

	for(K=0;K<NC;K++) 
	{
		int N_GBYE=0; 
		for(L2=0;L2<NQTL;L2++)
			if(GBYE_FIX[trait][L1][L2][K]!=0) N_GBYE=N_GBYE+1;

		T1=0,T2=0;
		for(L2=0;L2<NQTL;L2++)
			if(GBYE_FIX[trait][L1][L2][K]!=0) T1=T1+pow(GBYE_FIX[trait][L1][L2][K],2);

		for(J=0;J<N_GBYE+NU;J++)
		{
			ANORMAL(&U,&U0);
			T2=T2+U*U; 
		} 
     	
		T1=T1/(V_GBYE_FIX[trait][L1][K]*VP);T2=T2/(V_GBYE_FIX[trait][L1][K]*VP);
		V_GBYE_FIX[trait][L1][K]=(T1+NU*(NU-2.0)/NU*TAU)/T2;
	}

	return;
} */

void MT_GBYE_FixedCovariate_Variance(int L1,int L2,int K,int NU,double TAU,int trait)  //prior for V_GBYE is Inv-chisq(NU,(NU-2)/NU*TAU),E(V_GBYE)=TAU
{
	int J; double T1,T2,U,U0;

	T1=pow(GBYE_FIX[trait][L1][L2][K],2);
	T2=0;
	for(J=0;J<1+NU;J++)
	{
		ANORMAL(&U,&U0);
		T2=T2+U*U;
	}
	V_GBYE_FIX[trait][L1][L2][K]=(T1+NU*TAU)/T2;
	
	return;
}

        
 
//*******************************************************
//update nongenetic effects

void MT_FixedCovariate(int L,int trait)
{
	int I,ph; double G[NPHENO][NS1],Z,T1=0.0,T2=0.0,T3,U,U0;

    for(I=0;I<NS;I++)
    {			
		Z=COEF_FIX[I][L];
  		for(ph=0;ph<NPHENO;ph++)
  		{
        if(ph==trait)		G[ph][I]=GVALUE[ph][I]-Z*FIX[ph][L]; 
        else G[ph][I]=GVALUE[ph][I];		
      	T1=T1+Z*(Y[ph][I]-AMU[ph]-G[ph][I])*SIGMA[trait][ph];
  		}
     T2=T2+Z*Z;
    }
	T3=T2*SIGMA[trait][trait];	               //prior is uniform

	ANORMAL(&U,&U0); 	
	FIX[trait][L]=T1/T3+U/sqrt(T3);
        
	for(I=0;I<NS;I++) GVALUE[trait][I]=G[trait][I]+COEF_FIX[I][L]*FIX[trait][L];

	return;
}

void MT_RandomCovariate(int L, int trait)
{
	int I,K,PH; double G[NPHENO][NS1],T1,T2,T3,U,U0;
		
	for(I=0;I<NS;I++) 
	 for(PH=0;PH<NPHENO;PH++)
	   if(PH==trait)
         G[PH][I]=GVALUE[PH][I]-RAN[PH][L][(int)COEF_RAN[I][L]];
      else G[PH][I]=GVALUE[PH][I];

	for(K=0;K<NRAN[L];K++)
	{
		T1=0,T2=0;
		for(I=0;I<NS;I++)                    
			if(COEF_RAN[I][L]==K)
			{
  		for(PH=0;PH<NPHENO;PH++) T1=T1+(Y[PH][I]-AMU[PH]-G[PH][I])*SIGMA[trait][PH];
				T2=T2+1.0;
			}	 	
		T3=1/VRAN[trait][L]+T2*SIGMA[trait][trait];

		ANORMAL(&U,&U0);
		RAN[trait][L][K]=T1/T3+U/sqrt(T3);
	}	
	for(I=0;I<NS;I++) GVALUE[trait][I]=G[trait][I]+RAN[trait][L][(int)COEF_RAN[I][L]];
		
	return;	
}

void MT_RanVariance(int L,double VP,int trait)         //prior is Inv-chisq(NU,TAU) Now prior is uniform
{
	int K; double T1=0,T2=0,U,U0;

//	int NU=6; double TAU=VP/3;

	for(K=0;K<NRAN[L];K++) T1=T1+pow(RAN[trait][L][K],2);
	for(K=0;K<NRAN[L]-2;K++)
	{
        ANORMAL(&U,&U0);
        T2=T2+U*U; 
	} 
     	
	VRAN[trait][L]=T1/T2;      //see Gelman P301
	
	return;
}
    
//***********************************************************
//update QTL genotypes : generate IBD

void MT_QTLgenotype(int trait,int L,int NL,int QL,int I)    //trait not used for same.loc
{ 
	int K,PH,PH1; double SUMM[NG],temp[NPHENO],SUMM_SL[NPHENO][NG];              
  
  if(DiffLocation==1)
	for(K=0;K<NG;K++) SUMM[K]=MT_GenotypeSampling(I,L,trait,K,QL); 

  if(DiffLocation==0)
  for(PH=0;PH<NPHENO;PH++)
  	for(K=0;K<NG;K++) SUMM_SL[PH][K]=MT_GenotypeSampling(I,L,PH,K,QL); 
              
    double SUM=0.0;
    for(K=0;K<NG;K++)
    {
     if(MULTIPLE==0)
    	  {
         SUMM[K]=exp(-0.5*pow(Y[0][I]-(AMU[0]+SUMM[K]),2)/VE)*QPROB[NL][I][QL][K];
         SUM=SUM+SUMM[K]; 
         }
     if(MULTIPLE>0)
        {
         double sum=0.0;  
         if(DiffLocation==1)
         {
            for(PH=0;PH<NPHENO;PH++)
            {
              temp[PH]=0.0;
              for(PH1=0;PH1<NPHENO;PH1++)
               { 
               if(PH1==trait)  temp[PH] = temp[PH] + (Y[PH1][I]-AMU[PH1]-SUMM[K])*SIGMA[PH][PH1];
               else temp[PH] = temp[PH] + (Y[PH1][I]-AMU[PH1]-GVALUE[PH1][I])*SIGMA[PH][PH1];
                }
               if(PH==trait) sum = sum +temp[PH]*(Y[PH][I]-AMU[PH]-SUMM[K]);
               else sum = sum + temp[PH]*(Y[PH][I]-AMU[PH]-GVALUE[PH][I]);   
             }
           }
         if(DiffLocation==0)
         {
            for(PH=0;PH<NPHENO;PH++)
            {
              temp[PH]=0.0;
              for(PH1=0;PH1<NPHENO;PH1++)
              temp[PH] = temp[PH] + (Y[PH1][I]-AMU[PH1]-SUMM_SL[PH1][K])*SIGMA[PH][PH1];

              sum = sum +temp[PH]*(Y[PH][I]-AMU[PH]-SUMM_SL[PH][K]);
             }
           }
               
         SUMM[K]=exp(-0.5*sum)*QPROB[NL][I][QL][K];
         SUM=SUM+SUMM[K]; 
        }
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

void MT_QTLINHERITANCE(int L,int I,int trait)
{   
	int K;

	MT_QTLgenotype(trait,L,QCHR[trait][L],QLOC[trait][L],I); 

	GENO[trait][I][L]=IBD;
	Coefficient(IBD);
	for(K=0;K<NC;K++) COEF[trait][I][L][K]=X[K];
    
	PD1[L]=PD1[L]+log(PDD1+1e-20);
    PD2[L]=PD2[L]+log(PDD2+1e-20);		    

	return;
}
 
double LogLikelihood(double **yvalue,double **gvalue,double *mu,double **sigma)
{
 double likelihood=0.0,temp[NPHENO],sum=0.0;
 int i,ph,ph1;
 for(i=0;i<NS;i++)
 {
sum=0.0;
    for(ph=0;ph<NPHENO;ph++)
    {
      temp[ph]=0.0;
      for(ph1=0;ph1<NPHENO;ph1++)
      { 
      temp[ph] = temp[ph] + (yvalue[ph1][i]-mu[ph1]-gvalue[ph1][i])*sigma[ph][ph1];
       }
       sum = sum + temp[ph]*(yvalue[ph][i]-mu[ph]-gvalue[ph][i]);
     }  
    likelihood = likelihood + sum;
 }

 likelihood = -0.5*likelihood;
 return(likelihood);

} 
 
//*******************************************************
//UPDATE POSITIONS OF QTLS

void MT_QTLPOSITION(int L,int QLNEW,int trait)
{
	int I,ph,GENO1[NS1],K; 
   double **G;
   G = malloc(NPHENO*sizeof(double *));
   for(ph=0;ph<NPHENO;ph++) G[ph]= malloc(NS1*sizeof(double));

	double PROB0=0.0,PROB1=0.0,PD10=0.0,PD20=0.0;

 if(MULTIPLE==0)  for(I=0;I<NS;I++)  PROB0=PROB0-0.5*pow(Y[0][I]-AMU[0]-GVALUE[0][I],2)/VE;

 if(MULTIPLE>=1)  PROB0=LogLikelihood(Y,GVALUE,AMU,SIGMA);

for(I=0;I<NS;I++)
	{
		if(UPDATEGENO==1)
		{
			int KK=0;
			for(K=0;K<NG;K++)
				if(QPROB[QCHR[trait][L]][I][QLNEW][K]>0.99) 
				{
					GENO1[I]=K;
					KK=1;
				}    
			if(KK==0)
			{
        MT_QTLgenotype(trait,L,QCHR[trait][L],QLNEW,I);                             
				PD10=PD10+log(PDD1+1e-20);
				PD20=PD20+log(PDD2+1e-20);
				GENO1[I]=IBD;
			}
			int II=GENO1[I];
           if(MULTIPLE==0) 
             {
                G[0][I] = MT_GenotypeSampling(I,L,0,II,QLNEW); //0 is not 
                PROB1=PROB1-0.5*pow(Y[0][I]-AMU[0]-G[0][I],2)/VE;               
             }
          if(MULTIPLE==1) 
            {
            for(ph=0;ph<NPHENO;ph++)
        			G[ph][I]=MT_GenotypeSampling(I,L,ph,II,QLNEW);  //QLNEW is not used
        		}	
          if(MULTIPLE==2) 
            {
              for(ph=0;ph<NPHENO;ph++)
              {
                 if(DiffLocation==0)
                 {
                  if(GAMMA[ph][L]!=0) G[ph][I]=MT_GenotypeSampling(I,L,ph,II,QLNEW); 
                  else G[ph][I]=GVALUE[ph][I];
                  }
                 if(DiffLocation==1)
                 { 
                  if(GAMMA[ph][L]!=0 && ph==trait) 
                        G[ph][I]=MT_GenotypeSampling(I,L,ph,II,QLNEW); 
                  else G[ph][I]=GVALUE[ph][I];
                 } 
              }
        		}	
		} 
        
		if(UPDATEGENO==0) 
      {  
       if(MULTIPLE==0) 
            {
            G[0][I] = MT_GenotypeSampling(I,L,0,0,QLNEW); //0 is not 
            PROB1=PROB1-0.5*pow(Y[0][I]-AMU[0]-G[0][I],2)/VE;               
            }
          if(MULTIPLE==1) 
            for(ph=0;ph<NPHENO;ph++)
                    G[ph][I]=MT_GenotypeSampling(I,L,ph,0,QLNEW); //0 is not used

          if(MULTIPLE==2) 
            for(ph=0;ph<NPHENO;ph++)
            {
             if(DiffLocation==0)
             {
              if(GAMMA[ph][L]!=0) G[ph][I]=MT_GenotypeSampling(I,L,ph,0,QLNEW); 
              else G[ph][I]=GVALUE[ph][I];
              }
             if(DiffLocation==1)
              { 
              if(GAMMA[ph][L]!=0 && ph==trait) G[ph][I]=MT_GenotypeSampling(I,L,ph,0,QLNEW); 
              else G[ph][I]=GVALUE[ph][I];
              }
            }
        
      }  
	}   


if(MULTIPLE>=1)  PROB1 = LogLikelihood(Y,G,AMU,SIGMA);

    double S1=(PROB1-PROB0)+(PD10-PD1[L])+(PD2[L]-PD20); //Needed for genoupdate
         
 if(S1>log(RANDOM())) 
	{
		QLOC[trait][L]=QLNEW;
   for(I=0;I<NS;I++)          
		{
			if(UPDATEGENO==1)
			{
				GENO[trait][I][L]=GENO1[I];
				Coefficient(GENO1[I]);
			}  

      if(UPDATEGENO==0) MT_Coefficient0(I,L,trait,QLOC[trait][L]);

			for(K=0;K<NC;K++) COEF[trait][I][L][K]=X[K];

      for(ph=0;ph<NPHENO;ph++)  GVALUE[ph][I]=G[ph][I];				
  			
		}		
	 		
		PD1[L]=PD10,PD2[L]=PD20;
	}
		   for(ph=0;ph<NPHENO;ph++)  free(G[ph]);
		   free(G);
	return;
}                              




void QTLPOSITION_SameLocation(int L,int QLNEW)
{
	int I,ph,K;
// int GENO1[NS1];   
   double **G;
   G = malloc(NPHENO*sizeof(double *));
   for(ph=0;ph<NPHENO;ph++) G[ph]= malloc(NS1*sizeof(double));

	double PROB0=0.0;

 if(MULTIPLE==0)  for(I=0;I<NS;I++)  PROB0=PROB0-0.5*pow(Y[0][I]-AMU[0]-GVALUE[0][I],2)/VE;

 if(MULTIPLE>=1)  PROB0=LogLikelihood(Y,GVALUE,AMU,SIGMA);

            
    double PROB1=0.0,PD10=0.0,PD20=0.0;       
for(I=0;I<NS;I++)
	{
/*		if(UPDATEGENO==1)
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
				MT_QTLgenotype(L,QCHR[L],QLNEW,I);                             
           
				PD10=PD10+log(PDD1+1e-20);
				PD20=PD20+log(PDD2+1e-20);

				GENO1[I]=IBD;
			}
			int II=GENO1[I];
			G[I]=MT_GenotypeSampling(I,L,II,QLNEW);  //QLNEW is not used
		} 
        */
		if(UPDATEGENO==0) 
      {  
       if(MULTIPLE==0) 
            {
            G[0][I] = MT_GenotypeSampling(I,L,0,0,QLNEW); //0 is not 
            PROB1=PROB1-0.5*pow(Y[0][I]-AMU[0]-G[0][I],2)/VE;               
            }
          if(MULTIPLE==1) 
            for(ph=0;ph<NPHENO;ph++)
                    G[ph][I]=MT_GenotypeSampling(I,L,ph,0,QLNEW); //0 is not used

          if(MULTIPLE==2) 
            for(ph=0;ph<NPHENO;ph++)
            {
              if(GAMMA[ph][L]!=0) G[ph][I]=MT_GenotypeSampling(I,L,ph,0,QLNEW); 
              else G[ph][I]=GVALUE[ph][I];
            }
        
      }  
	}   


if(MULTIPLE>=1)  PROB1 = LogLikelihood(Y,G,AMU,SIGMA);

    double S1=(PROB1-PROB0)+(PD10-PD1[L])+(PD2[L]-PD20); //Needed for genoupdate
         
 if(S1>log(RANDOM())) 
	{
    QLOC[0][L]=QLNEW;
   for(I=0;I<NS;I++)          
		{
      if(UPDATEGENO==0) MT_Coefficient0(I,L,0,QLOC[0][L]);
			for(K=0;K<NC;K++) COEF[0][I][L][K]=X[K];

    for(ph=0;ph<NPHENO;ph++)  GVALUE[ph][I]=G[ph][I];				
      
/*			if(UPDATEGENO==1)
			{
				GENO[I][L]=GENO1[I];
				Coefficient(GENO1[I]);
			}  */
			
		}		
	 	
		PD1[L]=PD10,PD2[L]=PD20;
	}
		   for(ph=0;ph<NPHENO;ph++)  free(G[ph]);
		   free(G);
	return;
}                              
                     
//*************************************************************
//Sampling a position

int MT_SamplingOnePosition(int L,int trait)
{	
	int NL,PH,I,K,CHR0[NLG],T,TT,TTT,GRID00[TNGRID]; double R,GRID0[TNGRID];                                   
	for(I=0;I<NLG;I++)
	{
		CHRQTL[I]=0;
  if(DiffLocation==1)		
		for(K=0;K<NQTL;K++) CHRQTL[I]=CHRQTL[I]+(QCHR[trait][K]==I&&GAMMA[trait][K]!=0);

	if(DiffLocation==0)	
		for(K=0;K<NQTL;K++) 
  		{ T=0;
        if(QCHR[trait][K]==I)
    		  {
           for(PH=0;PH<NPHENO;PH++) T = T+GAMMA[PH][K];
            CHRQTL[I]=CHRQTL[I]+(T>=1);
          } 
      }    

		
	}
// CHRQTL contains the number of QTLs in each chromosome that are currently in the model
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
//CHR0 stores the number of putative QTL positions that is going to be considered for the new position.
	if(T!=0)
	{
		R=T*RANDOM();
		TT=0;
    for(I=0;I<NLG;I++) 
		{
			if(R>=TT&&R<(TT+CHR0[I])) QCHR[trait][L]=I;
			TT=TT+CHR0[I];
		}		
		NL=QCHR[trait][L];
		if(CHRQTL[NL]+1>CHR_NQTL[NL]||CHR0[NL]==0) T=0;
	}

	if(T!=0)
	{
		if(CHRQTL[NL]==0) QLOC[trait][L]=(int)(RANDOM()*NGRID[NL]);
		if(CHRQTL[NL]>0)
		{
			TT=0;
			if(DiffLocation==1)
    		{	for(K=0;K<NQTL;K++)
      				if(QCHR[trait][K]==NL&&GAMMA[trait][K]!=0) 
      				{
      					GRID0[TT]=GRID[NL][QLOC[trait][K]];
      					TT=TT+1;
      				}
      	}
      	
    	if(DiffLocation==0)
      	{	for(K=0;K<NQTL;K++)
      			if(QCHR[trait][K]==NL)
      				{
      				 int temp=0;
        		   for(PH=0;PH<NPHENO;PH++) temp = temp+GAMMA[PH][K];
      				 if(temp>=1)
               { 
      					GRID0[TT]=GRID[NL][QLOC[trait][K]];
      					TT=TT+1;
      					}
      				}
         }
			
//GRID0 contains the locations of all QTLs on chr NL (the new chromosome) currently in the model
			QLOC[trait][L]=(int)(RANDOM()*NGRID[NL]);
			TT=0;
			for(K=0;K<CHRQTL[NL];K++)
				if(fabs(GRID[NL][QLOC[trait][L]]-GRID0[K])<=DQQ[NL]) TT=TT+1;
	
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
// GRID00 is the vector of indicators, indicating the putative locations that can be considered for the new locus based on the constraint DQQ (minimum distance between flanking QTLs).
				if(TTT!=0)
				{
					R=RANDOM()*TTT;
					TT=0;
					for(I=0;I<NGRID[NL];I++)
					{
						if(R>=TT&&R<(TT+GRID00[I])) QLOC[trait][L]=I;
						TT=TT+GRID00[I];
					}	
				}
				if(TTT==0) T=0;
			}

			for(K=0;K<CHRQTL[NL];K++)
				if(fabs(GRID[NL][QLOC[trait][L]]-GRID0[K])<=DQQ[NL]) T=0;
		}
	}

	if(T!=0)
	{
		for(I=0;I<NS;I++)         
		{
			if(UPDATEGENO==1)
			{
				double PRR[NG];
				for(K=0;K<NG;K++) PRR[K]=QPROB[NL][I][QLOC[trait][L]][K];
				MULTINORMAL(PRR);                                     
				GENO[trait][I][L]=IBD;
				Coefficient(IBD);
			}
			if(UPDATEGENO==0) MT_Coefficient0(I,L,trait,QLOC[trait][L]);
		
    	for(K=0;K<NC;K++) COEF[trait][I][L][K]=X[K];
		}
	}

	return(T);
}

//***************************************************************
//Deleting QTL with all 0 effects
/*
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
} */

void MT_ZeroEffect1(int L,int trait)
{
	int L0; 

	GAMMA[trait][L]=(GAMMA_MAIN[trait][L]!=0);

	if(EPISTASIS==1)
	{
		for(L0=0;L0<=L-1;L0++) GAMMA[trait][L]=GAMMA[trait][L]+(GAMMA_EPISTASIS[trait][L0][L]!=0);
		for(L0=L+1;L0<NQTL;L0++) GAMMA[trait][L]=GAMMA[trait][L]+(GAMMA_EPISTASIS[trait][L][L0]!=0);   			
	}

	if(GBYE==1)
	{
	for(L0=0;L0<NFIXCOVA;L0++) 
	if(GBYE_FIX_INDEX[L0]==1) GAMMA[trait][L]=GAMMA[trait][L]+(GAMMA_GBYE[trait][L0][L]!=0.0);
	}

	return;
}

void MT_ZeroEffect2(int L1, int L2, int trait)
{
	int L0; 

	GAMMA[trait][L1]=(GAMMA_MAIN[trait][L1]!=0);
	GAMMA[trait][L2]=(GAMMA_MAIN[trait][L2]!=0);

	if(EPISTASIS==1)
	{
		for(L0=0;L0<=L1-1;L0++) GAMMA[trait][L1]=GAMMA[trait][L1]+(GAMMA_EPISTASIS[trait][L0][L1]!=0);
		for(L0=L1+1;L0<NQTL;L0++) GAMMA[trait][L1]=GAMMA[trait][L1]+(GAMMA_EPISTASIS[trait][L1][L0]!=0);
		
		for(L0=0;L0<=L2-1;L0++) GAMMA[trait][L2]=GAMMA[trait][L2]+(GAMMA_EPISTASIS[trait][L0][L2]!=0);
		for(L0=L2+1;L0<NQTL;L0++) GAMMA[trait][L2]=GAMMA[trait][L2]+(GAMMA_EPISTASIS[trait][L2][L0]!=0); 
	}

	if(GBYE==1)
	{
		for(L0=0;L0<NFIXCOVA;L0++) 
			if(GBYE_FIX_INDEX[L0]==1) 
			{
				GAMMA[trait][L1]=GAMMA[trait][L1]+(GAMMA_GBYE[trait][L0][L1]!=0.0);
				GAMMA[trait][L2]=GAMMA[trait][L2]+(GAMMA_GBYE[trait][L0][L2]!=0.0);
			}
	}

	return;
} 

//************************************************************************
//Update main effect for traditional multiple traits jointly for all phenotypes
/*
void MT_MainEffect_Traditional(int L, int K,FILE *File6)
{
	int I,J,K1,ph,ph1;
double **G,T2=0;
  double **sigma_beta_post,**tsigma1,**tsigma2,**chol;
  double *beta_ph,*beta_ph1,*newbeta;
 
         beta_ph = malloc(NPHENO*sizeof(double));
         beta_ph1 = malloc(NPHENO*sizeof(double));
         newbeta = malloc(NPHENO*sizeof(double));
         G = malloc(NPHENO*sizeof(double *));
         for (ph=0;ph<NPHENO;ph++) G[ph] = malloc(NS*sizeof(double));

  tsigma1 = malloc(NPHENO*sizeof(double *));
  tsigma2 = malloc(NPHENO*sizeof(double *));
  chol = malloc(NPHENO*sizeof(double *));
  sigma_beta_post = malloc(NPHENO*sizeof(double *));
  for (ph=0;ph<NPHENO;ph++)
  {   tsigma1[ph] = malloc(NPHENO*sizeof(double));
      tsigma2[ph] = malloc(NPHENO*sizeof(double));
      chol[ph] = malloc(NPHENO*sizeof(double));
      sigma_beta_post[ph] = malloc(NPHENO*sizeof(double));
   }
  for(ph=0;ph<NPHENO;ph++) beta_ph1[ph]=0.0,beta_ph[ph]=0.0;

 if(MULTIPLE==1)
  {  T2=0.0;
       for(I=0;I<NS;I++)
        { 
          double Z=COEF[I][L][K];
				  for(ph=0;ph<NPHENO;ph++)
				  {
            G[ph][I]=GVALUE[ph][I]-Z*MAIN[ph][L][K];
             
	   		  	beta_ph1[ph]=beta_ph1[ph]+Z*(Y[ph][I]-AMU[ph]-G[ph][I]);
		  		}
		  		T2=T2+Z*Z;
        }

   for(ph=0;ph<NPHENO;ph++)
      for(ph1=0;ph1<NPHENO;ph1++)
          {
           sigma_beta_post[ph][ph1]=SIGMA[ph][ph1]*T2;
      if(ph1==ph)  sigma_beta_post[ph][ph1] = sigma_beta_post[ph][ph1] + 1/VMAIN[ph][K];
          }

         
     INVERSE(sigma_beta_post,NPHENO,tsigma1);
     XprimeY(tsigma1,SIGMA,NPHENO,NPHENO,NPHENO,tsigma2);

       for(ph=0;ph<NPHENO;ph++)
        for(ph1=0;ph1<NPHENO;ph1++)
          beta_ph[ph] = beta_ph[ph]+tsigma2[ph][ph1]*beta_ph1[ph1]; 

  
  for(ph=0;ph<NPHENO;ph++)     
   for(int k=0;k<NPHENO;k++) chol[ph][k]=0;
      Cholesky(tsigma1,NPHENO,chol);
      Multivariate_RNORM(chol,beta_ph,NPHENO,newbeta);
      for(int trait=0;trait<NPHENO;trait++) 
      {
        MAIN[trait][L][K]=newbeta[trait];
     		for(I=0;I<NS;I++) GVALUE[trait][I]=G[trait][I]+COEF[I][L][K]*MAIN[trait][L][K];
      }

  }
 for(ph=0;ph<NPHENO;ph++)
  {
   free(G[ph]);
   free(tsigma1[ph]);
   free(tsigma2[ph]);
   free(chol[ph]);
   free(sigma_beta_post[ph]);
  }
   free(G);
   free(tsigma1);
   free(tsigma2);
   free(chol);
   free(sigma_beta_post);
   free(beta_ph);
   free(beta_ph1);
   free(newbeta);
	return;


} */

//************************************************************************
//Update main effect indicators

void MT_MainEffectIndicator_GROUP0(int L,int K,int trait)    
{
	int I,K1,ph,ph1;
double **G,**G2,BF_10=0,GAMMA10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100);
  double T1=0,T2=0,T3=0,F1=0,F2=0;
  double **sigma_beta_post,**tsigma1,**tsigma2,**chol;
  double *beta_ph,*beta_ph1,*newbeta;
 
         beta_ph = malloc(NPHENO*sizeof(double));
         beta_ph1 = malloc(NPHENO*sizeof(double));
         newbeta = malloc(NPHENO*sizeof(double));
         G2 = malloc(NPHENO*sizeof(double *));
         G = malloc(NPHENO*sizeof(double *));
         for (ph=0;ph<NPHENO;ph++)
         {  
         G[ph] = malloc(NS*sizeof(double));
         G2[ph] = malloc(NS*sizeof(double));
          }

  tsigma1 = malloc(NPHENO*sizeof(double *));
  tsigma2 = malloc(NPHENO*sizeof(double *));
  chol = malloc(NPHENO*sizeof(double *));
  sigma_beta_post = malloc(NPHENO*sizeof(double *));
  for (ph=0;ph<NPHENO;ph++)
  {   tsigma1[ph] = malloc(NPHENO*sizeof(double));
      tsigma2[ph] = malloc(NPHENO*sizeof(double));
      chol[ph] = malloc(NPHENO*sizeof(double));
      sigma_beta_post[ph] = malloc(NPHENO*sizeof(double));
   }
  for(ph=0;ph<NPHENO;ph++) beta_ph1[ph]=0.0,beta_ph[ph]=0.0;


	if(MULTIPLE==0)
	{
      	 for(I=0;I<NS;I++)                     
        	{
        		double Z=COEF[trait][I][L][K];
        		G[0][I]=GVALUE[0][I]-Z*MAIN[0][L][K];
        		T1=T1+Z*(Y[0][I]-AMU[0]-G[0][I]);
        		T2=T2+Z*Z;
      	   }	 	
      	 T3=1/VMAIN[0][L][K]+T2/VE;	
	
        	for(I=0;I<NS;I++)                     
        	{
      		double Z=COEF[trait][I][L][K];
      		F1=F1+pow(Y[0][I]-AMU[0]-G[0][I],2)/VE;
      		F2=F2+pow(Y[0][I]-AMU[0]-G[0][I]-Z*T1/(VE*T3),2)/VE;
        	}
	BF_10=-0.5*log(VMAIN[0][L][K])-0.5*pow(T1/(VE*T3),2)/VMAIN[0][L][K]-0.5*log(T3);
	BF_10=-0.5*F2+0.5*F1+BF_10;      
        	
  }
  if(MULTIPLE==1)
  {  T2=0.0;
       for(I=0;I<NS;I++)
        { 
          double Z=COEF[trait][I][L][K];
				  for(ph=0;ph<NPHENO;ph++)
				  {
            G[ph][I]=GVALUE[ph][I]-Z*MAIN[ph][L][K];
             
	   		  	beta_ph1[ph]=beta_ph1[ph]+Z*(Y[ph][I]-AMU[ph]-G[ph][I]);
		  		}
		  		T2=T2+Z*Z;
        }

   for(ph=0;ph<NPHENO;ph++)
      for(ph1=0;ph1<NPHENO;ph1++)
          {
           sigma_beta_post[ph][ph1]=SIGMA[ph][ph1]*T2;
      if(ph1==ph)  sigma_beta_post[ph][ph1] = sigma_beta_post[ph][ph1] + 1/VMAIN[ph][L][K];
          }

         
     INVERSE(sigma_beta_post,NPHENO,tsigma1);
     XprimeY(tsigma1,SIGMA,NPHENO,NPHENO,NPHENO,tsigma2);
       for(ph=0;ph<NPHENO;ph++)
        for(ph1=0;ph1<NPHENO;ph1++)
          beta_ph[ph] = beta_ph[ph]+tsigma2[ph][ph1]*beta_ph1[ph1]; 

   F1 = LogLikelihood(Y,G,AMU,SIGMA);
      for(I=0;I<NS;I++)
       {
       double Z=COEF[trait][I][L][K]; 
        for(ph=0;ph<NPHENO;ph++)
          G2[ph][I] = G[ph][I] + Z*beta_ph[ph];
        }  
  F2 = LogLikelihood(Y,G2,AMU,SIGMA);      
 
  double tt;
  tt=Determinant(sigma_beta_post,NPHENO); 

  for(ph=0;ph<NPHENO;ph++) BF_10=BF_10-0.5*log(VMAIN[ph][L][K]);
  BF_10 = BF_10 - 0.5*log(tt);    
  for(ph=0;ph<NPHENO;ph++) BF_10 = BF_10 - 0.5*beta_ph[ph]*beta_ph[ph]/VMAIN[ph][L][K];
	BF_10=BF_10+F2-F1;      
  }

  
if(MULTIPLE==2)
  {  T1=0,T2=0;

       for(I=0;I<NS;I++)
        { 
          double Z=COEF[trait][I][L][K];
				  for(ph=0;ph<NPHENO;ph++)
				  {
            if(ph==trait) G[ph][I]=GVALUE[ph][I]-Z*MAIN[ph][L][K];
            else G[ph][I]=GVALUE[ph][I];
             
	   		  	T1=T1+Z*(Y[ph][I]-AMU[ph]-G[ph][I])*SIGMA[trait][ph];
  		    }
		  		  T2=T2+Z*Z;            								
          
        }

   F1 = LogLikelihood(Y,G,AMU,SIGMA);
   T3=T2*SIGMA[trait][trait] + 1.0/VMAIN[trait][L][K];
      for(I=0;I<NS;I++)
       {
        for(ph=0;ph<NPHENO;ph++)
          {
        double Z=COEF[trait][I][L][K]; 
             if(ph==trait) G2[ph][I] = G[ph][I] + Z*T1/T3;
             else G2[ph][I]=G[ph][I];
          }
        }  
  F2 = LogLikelihood(Y,G2,AMU,SIGMA);   
  BF_10 = -0.5*log(VMAIN[trait][L][K]) - 0.5*log(T3);  
  BF_10 = BF_10 - 0.5*pow(T1/T3,2)/VMAIN[trait][L][K];
	BF_10=BF_10+F2-F1;      
	
  }

	GAMMA10=log(W_MAIN)-log(1-W_MAIN);
	if(GIBBS==1) GAMMA_1=BF_10+GAMMA10-log(1+exp(BF_10+GAMMA10));

	if(GIBBS==0&&MAIN[trait][L][K]==0) GAMMA_1=BF_10;
	if(GIBBS==0&&MAIN[trait][L][K]!=0) GAMMA_0=-BF_10;  
	      

	double R=log(RANDOM());
  if(MULTIPLE==0)
  {
  	if(R<GAMMA_1) 
  	{
  		double U,U0;
  		ANORMAL(&U,&U0);
  		MAIN[0][L][K]=T1/(VE*T3)+U/sqrt(T3);
  		GAMMA_MAIN[0][L]=1;
  		GAMMA[0][L]=1;
  		for(I=0;I<NS;I++) GVALUE[0][I]=G[0][I]+COEF[0][I][L][K]*MAIN[0][L][K];
  	}
  	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
  	{
  		MAIN[0][L][K]=0;
  		double SUM=0;
  		for(K1=0;K1<NC;K1++) SUM=SUM+fabs(MAIN[0][L][K1]);
  		if(SUM==0) 
  		{
  			GAMMA_MAIN[L]=0; 
  			MT_ZeroEffect1(L,0);
  		}
  		for(I=0;I<NS;I++) GVALUE[0][I]=G[0][I];
  	}
  }
  if(MULTIPLE==1)
  {
  	if(R<GAMMA_1) 
  	{
      Cholesky(tsigma1,NPHENO,chol);
      Multivariate_RNORM(chol,beta_ph,NPHENO,newbeta);
      for(ph=0;ph<NPHENO;ph++) 
      {   
      MAIN[ph][L][K]=newbeta[ph];
 //     MT_MainEffect(L,ph,1,K);
  		for(I=0;I<NS;I++) GVALUE[ph][I]=G[ph][I]+COEF[ph][I][L][K]*MAIN[ph][L][K];
  		GAMMA_MAIN[ph][L]=1;
  		GAMMA[ph][L]=1;
       }
  	}
  	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
  	{
  		double SUM=0;
  	  for(ph=0;ph<NPHENO;ph++)
      {
      MAIN[ph][L][K]=0;
  		for(K1=0;K1<NC;K1++) SUM=SUM+fabs(MAIN[ph][L][K1]);
  		}
  		if(SUM==0) 
  		{
      for(ph=0;ph<NPHENO;ph++)  
       {
        GAMMA_MAIN[ph][L]=0; 
  			MT_ZeroEffect1(L,ph);		
  			}
  		}
  	for(I=0;I<NS;I++) 
      for(ph=0;ph<NPHENO;ph++)
            GVALUE[ph][I]=G[ph][I];
          
  	}
  }
 if(MULTIPLE==2)
  {
  	if(R<GAMMA_1) 
  	{
   		double U,U0;
  		ANORMAL(&U,&U0);
  		MAIN[trait][L][K]=T1/T3+U/sqrt(T3);
   	for(I=0;I<NS;I++) GVALUE[trait][I]=G[trait][I]+COEF[trait][I][L][K]*MAIN[trait][L][K];

  		GAMMA_MAIN[trait][L]=1;
  		GAMMA[trait][L]=1;     
  	}
  	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
  	{
  		double SUM=0;
      MAIN[trait][L][K]=0;
  		for(K1=0;K1<NC;K1++) SUM=SUM+fabs(MAIN[trait][L][K1]);
  		if(SUM==0) 
  		{
        GAMMA_MAIN[trait][L]=0; 
  			MT_ZeroEffect1(L,trait);		
  		}
  	for(I=0;I<NS;I++) 
        GVALUE[trait][I]=G[trait][I];
          
  	}
  }

  
 for(ph=0;ph<NPHENO;ph++)
  {
   free(G[ph]);
   free(G2[ph]);
   free(tsigma1[ph]);
   free(tsigma2[ph]);
   free(chol[ph]);
   free(sigma_beta_post[ph]);
  }
   free(G);
   free(G2);
   free(tsigma1);
   free(tsigma2);
   free(chol);
   free(sigma_beta_post);
   free(beta_ph);
   free(beta_ph1);
   free(newbeta);
	return;
}

/*
void MT_MainEffectIndicator_GROUP1(int L)   
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
		MT_ZeroEffect1(L);
	}

	return;
}
*/
//************************************************************************
//Update epistatic effect indicators


void MT_EpistasisIndicator_GROUP0(int L1,int L2,int K1,int K2,int trait)  
{
	int I,K01,K02,ph,ph1;
	double **G,**G2,BF_10=0.0,GAMMA10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100);
  double T1=0,T2=0,T3=0,F1=0,F2=0;
  double **sigma_beta_post,**tsigma1,**tsigma2,**chol;
  double *beta_ph,*beta_ph1,*newbeta;
 
         beta_ph = malloc(NPHENO*sizeof(double));
         beta_ph1 = malloc(NPHENO*sizeof(double));
         newbeta = malloc(NPHENO*sizeof(double));
         G2 = malloc(NPHENO*sizeof(double *));
         G = malloc(NPHENO*sizeof(double *));
         for (ph=0;ph<NPHENO;ph++)
         {  
         G[ph] = malloc(NS*sizeof(double));
         G2[ph] = malloc(NS*sizeof(double));
          }

  tsigma1 = malloc(NPHENO*sizeof(double *));
  tsigma2 = malloc(NPHENO*sizeof(double *));
  chol = malloc(NPHENO*sizeof(double *));
  sigma_beta_post = malloc(NPHENO*sizeof(double *));
  for (ph=0;ph<NPHENO;ph++)
  {   tsigma1[ph] = malloc(NPHENO*sizeof(double));
      tsigma2[ph] = malloc(NPHENO*sizeof(double));
      chol[ph] = malloc(NPHENO*sizeof(double));
      sigma_beta_post[ph] = malloc(NPHENO*sizeof(double));
   }
  for(ph=0;ph<NPHENO;ph++) beta_ph1[ph]=0.0,beta_ph[ph]=0.0;

	if(W_EPISTASIS!=0)
	{
if(MULTIPLE==0)
{	
	for(I=0;I<NS;I++)                     
	{
		double Z=COEF[0][I][L1][K1]*COEF[0][I][L2][K2];
		G[0][I]=GVALUE[0][I]-Z*EPISTATIC[0][L1][L2][K1][K2];
		T1=T1+Z*(Y[0][I]-AMU[0]-G[0][I]);
		T2=T2+Z*Z;
	}	 	
	T3=1/VEPISTASIS[0][L1][L2][K1][K2]+T2/VE;	
	
	for(I=0;I<NS;I++)                     
	{
		double Z=COEF[0][I][L1][K1]*COEF[0][I][L2][K2];
		F1=F1+pow(Y[0][I]-AMU[0]-G[0][I],2)/VE;
		F2=F2+pow(Y[0][I]-AMU[0]-G[0][I]-Z*T1/(VE*T3),2)/VE;
	}

	BF_10=-0.5*log(VEPISTASIS[0][L1][L2][K1][K2])-0.5*pow(T1/(VE*T3),2)/VEPISTASIS[0][L1][L2][K1][K2]-0.5*log(T3);
	BF_10=-0.5*F2+0.5*F1+BF_10;      
}
  if(MULTIPLE==1)
  {  T2=0.0;
       for(I=0;I<NS;I++)
        { 
      		double Z=COEF[trait][I][L1][K1]*COEF[trait][I][L2][K2];
				  for(ph=0;ph<NPHENO;ph++)
				  {
            G[ph][I]=GVALUE[ph][I]-Z*EPISTATIC[ph][L1][L2][K1][K2];
	   		  	beta_ph1[ph]=beta_ph1[ph]+Z*(Y[ph][I]-AMU[ph]-G[ph][I]);
		  		}
		  		T2=T2+Z*Z;
        }

   for(ph=0;ph<NPHENO;ph++)
      for(ph1=0;ph1<NPHENO;ph1++)
          {
           sigma_beta_post[ph][ph1]=SIGMA[ph][ph1]*T2;
 if(ph1==ph)  
 sigma_beta_post[ph][ph1] = sigma_beta_post[ph][ph1] + 1.0/VEPISTASIS[ph][L1][L2][K1][K2];
          }

         
     INVERSE(sigma_beta_post,NPHENO,tsigma1);
     XprimeY(tsigma1,SIGMA,NPHENO,NPHENO,NPHENO,tsigma2);
       for(ph=0;ph<NPHENO;ph++)
        for(ph1=0;ph1<NPHENO;ph1++)
          beta_ph[ph] = beta_ph[ph]+tsigma2[ph][ph1]*beta_ph1[ph1]; 

   F1 = LogLikelihood(Y,G,AMU,SIGMA);
      for(I=0;I<NS;I++)
       {
     		double Z=COEF[trait][I][L1][K1]*COEF[trait][I][L2][K2];
        for(ph=0;ph<NPHENO;ph++)
          G2[ph][I] = G[ph][I] + Z*beta_ph[ph];
        }  
  F2 = LogLikelihood(Y,G2,AMU,SIGMA);      
 
  double tt;
  tt=Determinant(sigma_beta_post,NPHENO); 

  for(ph=0;ph<NPHENO;ph++) BF_10=BF_10-0.5*log(VEPISTASIS[ph][L1][L2][K1][K2]);
  BF_10 = BF_10 - 0.5*log(tt);    
  for(ph=0;ph<NPHENO;ph++) 
  BF_10 = BF_10 - 0.5*beta_ph[ph]*beta_ph[ph]/VEPISTASIS[ph][L1][L2][K1][K2];
	
  BF_10=BF_10+F2-F1;      
  }

if(MULTIPLE==2)
  {  T1=0,T2=0;
       for(I=0;I<NS;I++)
        {
      		double Z=COEF[trait][I][L1][K1]*COEF[trait][I][L2][K2];
				  for(ph=0;ph<NPHENO;ph++)
				  {
            if(ph==trait) G[ph][I]=GVALUE[ph][I]-Z*EPISTATIC[ph][L1][L2][K1][K2];
            else G[ph][I]=GVALUE[ph][I];

	   		  	T1=T1+Z*(Y[ph][I]-AMU[ph]-G[ph][I])*SIGMA[trait][ph];
  		    }
		  		  T2=T2+Z*Z;

        }

   F1 = LogLikelihood(Y,G,AMU,SIGMA);
   T3=T2*SIGMA[trait][trait] + 1.0/VEPISTASIS[trait][L1][L2][K1][K2];

      for(I=0;I<NS;I++)
       {
        for(ph=0;ph<NPHENO;ph++)
          {
        		double Z=COEF[trait][I][L1][K1]*COEF[trait][I][L2][K2];
             if(ph==trait) G2[ph][I] = G[ph][I] + Z*T1/T3;
             else G2[ph][I]=G[ph][I];
          }
        }
  F2 = LogLikelihood(Y,G2,AMU,SIGMA);
  BF_10 = -0.5*log(VEPISTASIS[trait][L1][L2][K1][K2]) - 0.5*log(T3);
  BF_10 = BF_10 - 0.5*pow(T1/T3,2)/VEPISTASIS[trait][L1][L2][K1][K2];
	BF_10=BF_10+F2-F1;

  }

	GAMMA10=log(W_EPISTASIS)-log(1-W_EPISTASIS);
	if(GIBBS==1) GAMMA_1=BF_10+GAMMA10-log(1+exp(BF_10+GAMMA10));
	if(GIBBS==0 && EPISTATIC[trait][L1][L2][K1][K2]==0) GAMMA_1=BF_10;
	if(GIBBS==0 && EPISTATIC[trait][L1][L2][K1][K2]!=0) GAMMA_0=-BF_10;

	double R=log(RANDOM());

   if(MULTIPLE==0)
   {
    	if(R<GAMMA_1) 
    	{
    		double U,U0;
    		ANORMAL(&U,&U0);
    		EPISTATIC[0][L1][L2][K1][K2]=T1/(VE*T3)+U/sqrt(T3);
    		GAMMA_EPISTASIS[0][L1][L2]=1;
    		GAMMA[0][L1]=1,GAMMA[0][L2]=1;
    		for(I=0;I<NS;I++) 
        GVALUE[0][I]=G[0][I]+COEF[0][I][L1][K1]*COEF[0][I][L2][K2]*EPISTATIC[0][L1][L2][K1][K2];
    	}
    	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
    	{
    		EPISTATIC[0][L1][L2][K1][K2]=0;
    		double SUM=0;
    		for(K01=0;K01<NC;K01++) 
    			for(K02=0;K02<NC;K02++) SUM=SUM+fabs(EPISTATIC[0][L1][L2][K01][K02]);
    		if(SUM==0) 
    		{
    			GAMMA_EPISTASIS[0][L1][L2]=0;  
    			MT_ZeroEffect2(L1,L2,0);
    		}
    		for(I=0;I<NS;I++) GVALUE[0][I]=G[0][I];
    	}
    }

  if(MULTIPLE==1)
  {
  	if(R<GAMMA_1) 
  	{
      Cholesky(tsigma1,NPHENO,chol);
      Multivariate_RNORM(chol,beta_ph,NPHENO,newbeta);
      for(ph=0;ph<NPHENO;ph++) 
      {   
      EPISTATIC[ph][L1][L2][K1][K2]=newbeta[ph];
   	for(I=0;I<NS;I++) 
GVALUE[ph][I]=G[ph][I]+COEF[ph][I][L1][K1]*COEF[ph][I][L2][K2]*EPISTATIC[ph][L1][L2][K1][K2];

  		GAMMA_EPISTASIS[ph][L1][L2]=1;
   		GAMMA[ph][L1]=1,GAMMA[ph][L2]=1;
       }
  	}
  	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
  	{
  		double SUM=0;
  	  for(ph=0;ph<NPHENO;ph++)
      {
 	  	EPISTATIC[ph][L1][L2][K1][K2]=0;
    		for(K01=0;K01<NC;K01++) 
    			for(K02=0;K02<NC;K02++) SUM=SUM+fabs(EPISTATIC[ph][L1][L2][K01][K02]);
  		}
  		if(SUM==0) 
  		{
      for(ph=0;ph<NPHENO;ph++)  
       {
    			GAMMA_EPISTASIS[ph][L1][L2]=0;  
    			MT_ZeroEffect2(L1,L2,ph);
  			}
  		}
  	for(I=0;I<NS;I++) 
      for(ph=0;ph<NPHENO;ph++)
            GVALUE[ph][I]=G[ph][I];
          
  	}
  }

 if(MULTIPLE==2)
  {
  	if(R<GAMMA_1)
  	{
   		double U,U0;
  		ANORMAL(&U,&U0);
   		EPISTATIC[trait][L1][L2][K1][K2]=T1/T3+U/sqrt(T3);
   	for(I=0;I<NS;I++) 
     GVALUE[trait][I]=G[trait][I]+COEF[trait][I][L1][K1]*COEF[trait][I][L2][K2]*EPISTATIC[trait][L1][L2][K1][K2];

  		GAMMA_EPISTASIS[trait][L1][L2]=1;
   		GAMMA[trait][L1]=1,GAMMA[trait][L2]=1;
		
  	}
  	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
  	{
  		EPISTATIC[trait][L1][L2][K1][K2]=0;
   		double SUM=0;
    		for(K01=0;K01<NC;K01++) 
    			for(K02=0;K02<NC;K02++) SUM=SUM+fabs(EPISTATIC[trait][L1][L2][K01][K02]);
  		if(SUM==0)
  		{
    			GAMMA_EPISTASIS[trait][L1][L2]=0;  
    			MT_ZeroEffect2(L1,L2,trait);
  		}
  	for(I=0;I<NS;I++)
        GVALUE[trait][I]=G[trait][I];

  	}
  }
    
    
	}  //end of w_epistasis != 0
 /*
	if(W_EPISTASIS==0&&EPISTATIC[trait][L1][L2][K1][K2]!=0)
	{
		for(I=0;I<NS;I++) 
    GVALUE[trait][I]=GVALUE[trait][I]-COEF[trait][I][L1][K1]*COEF[trait][I][L2][K2]*EPISTATIC[trait][L1][L2][K1][K2];
		EPISTATIC[trait][L1][L2][K1][K2]=0;
		double SUM=0;
		for(K01=0;K01<NC;K01++) 
			for(K02=0;K02<NC;K02++) SUM=SUM+fabs(EPISTATIC[trait][L1][L2][K01][K02]);
		if(SUM==0) 
		{
			GAMMA_EPISTASIS[trait][L1][L2]=0;  
			MT_ZeroEffect2(L1,L2,trait);
		}
	}*/

 for(ph=0;ph<NPHENO;ph++)
  {
   free(G[ph]);
   free(G2[ph]);
   free(tsigma1[ph]);
   free(tsigma2[ph]);
   free(chol[ph]);
   free(sigma_beta_post[ph]);
  }
   free(G);
   free(G2);
   free(tsigma1);
   free(tsigma2);
   free(chol);
   free(sigma_beta_post);
   free(beta_ph);
   free(beta_ph1);
   free(newbeta);
	return;	
}

/*
void MT_EpistasisIndicator_GROUP1(int L1,int L2)    
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
		MT_ZeroEffect2(L1,L2);
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
		MT_ZeroEffect2(L1,L2);
	}

	return;	
}

*/
//************************************************************************
//Update g by e fixed effect indicators


void MT_GBYE_FIX_Indicator_GROUP0(int L1,int L2,int K,int trait)  //use the constraint model
{
	int I,K1,ph,ph1;
	double **G,**G2,BF_10=0.0,GAMMA10,GAMMA_1=-(1e+100),GAMMA_0=-(1e+100);

  double T1=0,T2=0,T3=0,F1=0,F2=0;
  double **sigma_beta_post,**tsigma1,**tsigma2,**chol;
  double *beta_ph,*beta_ph1,*newbeta;
 
         beta_ph = malloc(NPHENO*sizeof(double));
         beta_ph1 = malloc(NPHENO*sizeof(double));
         newbeta = malloc(NPHENO*sizeof(double));
         G2 = malloc(NPHENO*sizeof(double *));
         G = malloc(NPHENO*sizeof(double *));
         for (ph=0;ph<NPHENO;ph++)
         {  
         G[ph] = malloc(NS*sizeof(double));
         G2[ph] = malloc(NS*sizeof(double));
          }

  tsigma1 = malloc(NPHENO*sizeof(double *));
  tsigma2 = malloc(NPHENO*sizeof(double *));
  chol = malloc(NPHENO*sizeof(double *));
  sigma_beta_post = malloc(NPHENO*sizeof(double *));
  for (ph=0;ph<NPHENO;ph++)
  {   tsigma1[ph] = malloc(NPHENO*sizeof(double));
      tsigma2[ph] = malloc(NPHENO*sizeof(double));
      chol[ph] = malloc(NPHENO*sizeof(double));
      sigma_beta_post[ph] = malloc(NPHENO*sizeof(double));
   }
  for(ph=0;ph<NPHENO;ph++) beta_ph1[ph]=0.0,beta_ph[ph]=0.0;

  if(MULTIPLE==1)
  {  T2=0.0;
       for(I=0;I<NS;I++)
        { 
      		double Z=COEF_FIX[I][L1]*COEF[trait][I][L2][K];
				  for(ph=0;ph<NPHENO;ph++)
				  {
            G[ph][I]=GVALUE[ph][I]-Z*GBYE_FIX[ph][L1][L2][K];
             
	   		  	beta_ph1[ph]=beta_ph1[ph]+Z*(Y[ph][I]-AMU[ph]-G[ph][I]);
		  		}
		  		T2=T2+Z*Z;
        }

   for(ph=0;ph<NPHENO;ph++)
      for(ph1=0;ph1<NPHENO;ph1++)
          {
           sigma_beta_post[ph][ph1]=SIGMA[ph][ph1]*T2;
      if(ph1==ph)  
      sigma_beta_post[ph][ph1] = sigma_beta_post[ph][ph1] + 1.0/V_GBYE_FIX[ph][L1][L2][K];
          }

         
     INVERSE(sigma_beta_post,NPHENO,tsigma1);
     XprimeY(tsigma1,SIGMA,NPHENO,NPHENO,NPHENO,tsigma2);
       for(ph=0;ph<NPHENO;ph++)
        for(ph1=0;ph1<NPHENO;ph1++)
          beta_ph[ph] = beta_ph[ph]+tsigma2[ph][ph1]*beta_ph1[ph1]; 

   F1 = LogLikelihood(Y,G,AMU,SIGMA);
      for(I=0;I<NS;I++)
       {
    		double Z=COEF_FIX[I][L1]*COEF[trait][I][L2][K];
        for(ph=0;ph<NPHENO;ph++)
          G2[ph][I] = G[ph][I] + Z*beta_ph[ph];
        }  
  F2 = LogLikelihood(Y,G2,AMU,SIGMA);      
 
  double tt;
  tt=Determinant(sigma_beta_post,NPHENO); 

  for(ph=0;ph<NPHENO;ph++) BF_10=BF_10-0.5*log(V_GBYE_FIX[ph][L1][L2][K]);
  BF_10 = BF_10 - 0.5*log(tt);    
  for(ph=0;ph<NPHENO;ph++) 
  BF_10 = BF_10 - 0.5*beta_ph[ph]*beta_ph[ph]/V_GBYE_FIX[ph][L1][L2][K];
	
  BF_10=BF_10+F2-F1;      
  }


	if(MULTIPLE==2)
	{   
  T1=0,T2=0;
	for(I=0;I<NS;I++)                     
	{
		double Z=COEF_FIX[I][L1]*COEF[trait][I][L2][K];
	  for(ph=0;ph<NPHENO;ph++)
		  {
       if(ph==trait) G[ph][I]=GVALUE[ph][I]-Z*GBYE_FIX[ph][L1][L2][K];
       else G[ph][I]=GVALUE[ph][I];
	  	T1=T1+Z*(Y[ph][I]-AMU[ph]-G[ph][I])*SIGMA[trait][ph];
	    }
 		  T2=T2+Z*Z;            								
	}	 	

  F1 = LogLikelihood(Y,G,AMU,SIGMA);
  T3=T2*SIGMA[trait][trait] + 1.0/V_GBYE_FIX[trait][L1][L2][K];
	
  for(I=0;I<NS;I++)
    {
    for(ph=0;ph<NPHENO;ph++)
      {
  		double Z=COEF_FIX[I][L1]*COEF[trait][I][L2][K];
       if(ph==trait) G2[ph][I] = G[ph][I] + Z*T1/T3;
       else G2[ph][I]=G[ph][I];
       }
    }  
  F2 = LogLikelihood(Y,G2,AMU,SIGMA);   

  BF_10 = -0.5*log(V_GBYE_FIX[trait][L1][L2][K]) - 0.5*log(T3);  
  BF_10 = BF_10 - 0.5*pow(T1/T3,2)/V_GBYE_FIX[trait][L1][L2][K];
	BF_10=BF_10+F2-F1;      
 }
 
	GAMMA10=log(W_GBYE)-log(1-W_GBYE);
	if(GIBBS==1) GAMMA_1=BF_10+GAMMA10-log(1+exp(BF_10+GAMMA10));
	if(GIBBS==0&&GBYE_FIX[trait][L1][L2][K]==0) GAMMA_1=BF_10;
	if(GIBBS==0&&GBYE_FIX[trait][L1][L2][K]!=0) GAMMA_0=-BF_10;

	double R=log(RANDOM());
  if(MULTIPLE==1)
  {
  	if(R<GAMMA_1) 
  	{
      Cholesky(tsigma1,NPHENO,chol);
      Multivariate_RNORM(chol,beta_ph,NPHENO,newbeta);
      for(ph=0;ph<NPHENO;ph++) 
      {   
      GBYE_FIX[ph][L1][L2][K]=newbeta[ph];

  		for(I=0;I<NS;I++) 
      GVALUE[ph][I]=G[ph][I]+COEF_FIX[I][L1]*COEF[ph][I][L2][K]*GBYE_FIX[ph][L1][L2][K];

  		GAMMA_GBYE[ph][L1][L2]=1;
  		GAMMA[ph][L2]=1;
       }
  	}
  	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
  	{
  		double SUM=0;
  	  for(ph=0;ph<NPHENO;ph++)
      {
  		GBYE_FIX[ph][L1][L2][K]=0;
  		for(K1=0;K1<NC;K1++) SUM=SUM+fabs(GBYE_FIX[ph][L1][L2][K]);
  		}
  		if(SUM==0) 
  		{
      for(ph=0;ph<NPHENO;ph++)  
       {
        GAMMA_GBYE[ph][L1][L2]=0; 
  			MT_ZeroEffect1(L2,ph);		
  			}
  		}
  	for(I=0;I<NS;I++) 
      for(ph=0;ph<NPHENO;ph++)
            GVALUE[ph][I]=G[ph][I];
          
  	}
  }
  if(MULTIPLE==2)
	{
  	if(R<GAMMA_1) 
  	{
  		double U,U0;
  		ANORMAL(&U,&U0);
   		GBYE_FIX[trait][L1][L2][K]=T1/T3+U/sqrt(T3);		
  		GAMMA_GBYE[trait][L1][L2]=1;
  		GAMMA[trait][L2]=1;
  
  		for(I=0;I<NS;I++) 
  GVALUE[trait][I]=G[trait][I]+COEF_FIX[I][L1]*COEF[trait][I][L2][K]*GBYE_FIX[trait][L1][L2][K];
  	}
  	if( (R>=GAMMA_1&&GIBBS==1) || R<GAMMA_0 )
  	{
  		GBYE_FIX[trait][L1][L2][K]=0;
  		double SUM=0;
  		for(K1=0;K1<NC;K1++) SUM=SUM+fabs(GBYE_FIX[trait][L1][L2][K1]);
  		if(SUM==0) 
  		{
  			GAMMA_GBYE[trait][L1][L2]=0;  
  			MT_ZeroEffect1(L2,trait);
  		}
  		for(I=0;I<NS;I++) GVALUE[trait][I]=G[trait][I];
  	}
  }

 for(ph=0;ph<NPHENO;ph++)
  {
   free(G[ph]);
   free(G2[ph]);
   free(tsigma1[ph]);
   free(tsigma2[ph]);
   free(chol[ph]);
   free(sigma_beta_post[ph]);
  }
   free(G);
   free(G2);
   free(tsigma1);
   free(tsigma2);
   free(chol);
   free(sigma_beta_post);
   free(beta_ph);
   free(beta_ph1);
   free(newbeta);

	return;	
}
/*
void MT_GBYE_FIX_Indicator_GROUP1(int L1,int L2)  
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
		MT_ZeroEffect1(L2);
	}

	return;	
}
*/  
//*******************************************************************************    


// ************************************************************

void multipleTraitsMCMC()
{

int I,J,L,K,K1,K2,L0,L1,L2,PH;     
double **ybar,**sigma;
  
if(SEED==0)	srand(time(NULL));
else srand(SEED);
rand();    

//*************************************************************************
// calculating phenotypic mean and variance, and assign initial values

double YBAR[NPHENO],Y2BAR[NPHENO],VP[NPHENO];

sigma=malloc(NPHENO*sizeof(double *));
for(I=0;I<NPHENO;I++) sigma[I] = malloc(NPHENO*sizeof(double));

if(CATEGORY==1)
{
for(PH=0;PH<NPHENO;PH++)
{
	YBAR[PH]=0.0,Y2BAR[PH]=0.0;
	for(I=0;I<NS;I++)
	{
		YBAR[PH]=YBAR[PH]+Y[PH][I];
		Y2BAR[PH]=Y2BAR[PH]+Y[PH][I]*Y[PH][I];
	}
	YBAR[PH]=YBAR[PH]/NS;
	AMU[PH]=YBAR[PH];
	VP[PH]=(1.0/(NS-1))*(Y2BAR[PH]-NS*pow(YBAR[PH],2));
}	 
  VE=VP[0];  //Initializing Residual Variance

if(SPH==1)
	for(PH=0;PH<NPHENO;PH++)
  		for(I=0;I<NS;I++) Y[PH][I]=(Y[PH][I]-YBAR[PH])/sqrt(VP[PH]); //standardizing the phenotype

if(MULTIPLE>=1)
{
  ybar=malloc(NPHENO*sizeof(double *));
  for(I=0;I<NPHENO;I++) ybar[I] = malloc(NS1*sizeof(double));
  
 for(I=0;I<NPHENO;I++) AMU[I]=0;
  for(J=0;J<NPHENO;J++)
    for(I=0;I<NS;I++)
        AMU[J]=AMU[J]+ Y[J][I]/NS;

          for(J=0;J<NPHENO;J++) YBAR[J]=AMU[J];
  
  for(J=0;J<NPHENO;J++) 
    for(I=0;I<NS;I++) ybar[J][I] = (Y[J][I] - AMU[J])/sqrt(NS-1);
            
   XprimeX(ybar,NPHENO,NS,sigma);
   INVERSE(sigma,NPHENO,SIGMA);

  }


/*	AMU=YBAR, VE=VP;            //initial values for AMU and VE
	for(I=0;I<NRANCOVA;I++) VRAN[I]=VP;  //initial values for VRAN[I]
	*/
	
/*  for(I=0;I<NPHENO;I++)
  {
  for(J=0;J<NPHENO;J++) VP[J]=SIGMA[I][J];
  MT_Mean(YBAR,VP,I,fp);
  }*/
 
 
}


/*
if(CATEGORY!=1)
{
	AMU=0, VE=1.0;                        //initial values for AMU and VE
	for(I=0;I<NRANCOVA;I++) VRAN[I]=1.0;  //initial values for VRAN[I]
}
*/
//******************************************************************************
//For specify the prior variances of QTL effects and g by e interactions

double CC[NG];
if(CROSS==2)
{
	CC[0]=1.0/2, CC[1]=1.0/4;
}
else CC[0]=1.0/4;
for(K=0;K<NC;K++)
 for(PH=0;PH<NPHENO;PH++) 
  for(L=0;L<NQTL;L++)
    VMAIN[PH][L][K]=1.0/CC[K];

for(PH=0;PH<NPHENO;PH++)
{	
 		for(L1=0;L1<NQTL-1;L1++)
			for(L2=L1+1;L2<NQTL;L2++)
        for(K1=0;K1<NC;K1++)
      		for(K2=0;K2<NC;K2++) 
          VEPISTASIS[PH][L1][L2][K1][K2]=1.0/(CC[K1]*CC[K2]);
	}


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

  for(PH=0;PH<NPHENO;PH++)
 		for(L1=0;L1<NFIXCOVA;L1++)
			for(L2=0;L2<NQTL;L2++)
  			for(K=0;K<NC;K++) V_GBYE_FIX[PH][L1][L2][K]=1.0/(V_FIX[L1]*CC[K]);

//******************************************************************************
//Give initial threshold values for ordinal traits
/*
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
*/
//******************************************************************************
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

FILE *File7;
File7=fopen(sigmafile,"w");


// ***********************************************************                         
// ITERATION STARTS HERE

int ITER,ITER1;
for(ITER=0;ITER<NITER+(int)(1.0*NBURNIN/NTHIN);ITER++)
{ 
for(ITER1=0;ITER1<NTHIN;ITER1++)
{

//***********************************************************
//UPDATING THE VALUES OF THE LIABILITY

/// @warning: Review if this needs to be commented out. Also, call to TrunNormal doesnt have sufficient arguments
//if(CATEGORY!=1) 
	//for(I=0;I<NS;I++) Y[0][I]=TrunNormal(W[I],AMU[0]+GVALUE[0][I],VE);      

//********************************************************************************* 
//UPDATE PARAMETERS
 
for(PH=0;PH<NPHENO;PH++)  MT_Mean(YBAR[PH],VP[PH],PH);


//if(CATEGORY!=2 && MULTIPLE==0) ResidualVariance();

if(CATEGORY==1 && MULTIPLE>=1) ResidualVariance_MultipleTraits();

int NU=6; double H=0.1,S=2,TAU; //NU=degrees of freedom, H=heritability, TAU=scale

for(L=0;L<NQTL;L++)
    for(PH=0;PH<NPHENO;PH++)    
	      if(GAMMA_MAIN[PH][L]!=0.0)
        {
         MT_MainEffect(L,PH,0,0); // 0 0 not used here	
    			for(K=0;K<NC;K++)
    				if(MAIN[PH][L][K]!=0)
       			{
    					TAU=(NU-2)*H*VP[PH]/(NU*CC[K]);	
    					MT_MainVariance(L,K,NU,TAU,PH);
    	   		}		
         }

	
  /*  for(L=0;L<NQTL;L++)
	      if(GAMMA_MAIN[L]!=0.0) 
          for(K=0;K<NC;K++)    
            if(MAINTRUE==0) MT_MainEffect_Traditional(L,K); 
*/	

/*if(EPISTASIS==1)
{  for(PH=0;PH<NPHENO;PH++)
  	for(L1=0;L1<NQTL-1;L1++)
	   	for(L2=L1+1;L2<NQTL;L2++)
		   	if(GAMMA_EPISTASIS[PH][L1][L2]!=0)
			           MT_EpistaticEffect(L1,L2,PH);	
} */

	if(EPISTASIS==1)
	{
   for(PH=0;PH<NPHENO;PH++)
 		for(L1=0;L1<NQTL-1;L1++)
			for(L2=L1+1;L2<NQTL;L2++)
				if(GAMMA_EPISTASIS[PH][L1][L2]!=0.0) 
					{
           MT_EpistaticEffect(L1,L2,PH);	
              for(K1=0;K1<NC;K1++)
    						for(K2=0;K2<NC;K2++) 
    							if(EPISTATIC[PH][L1][L2][K1][K2]!=0.0)
    							{
    								TAU=(NU-2)*H*VP[PH]/(NU*CC[K1]*CC[K2]);
    								MT_EpistaticVariance(L1,L2,K1,K2,NU,TAU,PH);
    							}
					}		
	}

/*
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
} */

//for(K=0;K<NC;K++) VMAIN[K]=(float)NS;

/*
if(EPISTASIS==1) 
{
for(PH=0;PH<NPHENO;PH++)
{	for(K1=0;K1<NC;K1++)
		for(K2=0;K2<NC;K2++) VEPISTASIS[PH][K1][K2]=1.0/(CC[K1]*CC[K2]);
	MT_EpistaticVariance(VP[PH],PH);
	}
} */

/*
if(GBYE==1)
{  
	for(L1=0;L1<NFIXCOVA;L1++)
		if(GBYE_FIX_INDEX[L1]==1)
		{
		for(PH=0;PH<NPHENO;PH++)
		{
			for(K=0;K<NC;K++) V_GBYE_FIX[PH][L1][K]=1.0/(V_FIX[L1]*CC[K]);
			MT_GBYE_FixedCovariate_Variance(L1,VP[PH],PH);
			}
		}
} */  



if(ENV_FACTOR==1)
{ 
for(PH=0;PH<NPHENO;PH++)
	for(L=0;L<NRANCOVA;L++) 
	{
		MT_RandomCovariate(L,PH);		
		MT_RanVariance(L,VP[PH],PH);
	}
for(PH=0;PH<NPHENO;PH++)	
	for(L=0;L<NFIXCOVA;L++) MT_FixedCovariate(L,PH);		
}
/*
if(GBYE==1)
{ 
	for(L1=0;L1<NFIXCOVA;L1++)
	 for(PH=0;PH<NPHENO;PH++)
  	if(GBYE_FIX_INDEX[L1]==1)
			for(L2=0;L2<NQTL;L2++)
				if(GAMMA[PH][L2]!=0.0 && GAMMA_GBYE[PH][L1][L2]!=0.0) 
        MT_GBYE_FixedCovariate(L1,L2,PH);	
} */

	if(GBYE==1)
	{
	 for(PH=0;PH<NPHENO;PH++)
		for(L1=0;L1<NFIXCOVA;L1++)
			if(GBYE_FIX_INDEX[L1]==1)
				for(L2=0;L2<NQTL;L2++)
					if(GAMMA[PH][L2]!=0.0&&GAMMA_GBYE[PH][L1][L2]!=0.0)
					{
            MT_GBYE_FixedCovariate(L1,L2,PH);	
          	for(K=0;K<NC;K++) 
							if(GBYE_FIX[PH][L1][L2][K]!=0)
							{
								TAU=(NU-2)*H*VP[PH]/(NU*V_FIX[L1]*CC[K]);
				        MT_GBYE_FixedCovariate_Variance(L1,L2,K,NU,TAU,PH);
							}
					}		
	 }
  
//************************************************************************************
//UPDATE THE QTL INHERITANCE OF NON-FOUNDERS AND THE GENOTYPIC VALUES

if(UPDATEGENO==1)
{
  if(DiffLocation==1)
  {
    for(PH=0;PH<NPHENO;PH++)
      for(L=0;L<NQTL;L++) 
        	if(GAMMA[PH][L]!=0) 
        	{
        		PD1[L]=0.0,PD2[L]=0.0;
        		for(I=0;I<NS;I++)
        		{
         			int K,KK=0;
        			for(K=0;K<NG;K++)
        			   if(QPROB[QCHR[PH][L]][I][QLOC[PH][L]][K]>0.99) 
        			   {
        				   IBD=K;
        				   KK=1;
        			   } 
        			if(KK==0)
                {
                MT_QTLgenotype(PH,L,QCHR[PH][L],QLOC[PH][L],I); 
              	PD1[L]=PD1[L]+log(PDD1+1e-20);
                PD2[L]=PD2[L]+log(PDD2+1e-20);		    
                }
                if(IBD != GENO[PH][I][L])
                {
        				GVALUE[PH][I]=MT_GenotypeSampling(I,L,PH,IBD,0);
              	GENO[PH][I][L]=IBD;
              	Coefficient(IBD);
              	for(K=0;K<NC;K++) COEF[PH][I][L][K]=X[K];
                }
  
        		}
        	}
  }
  if(DiffLocation==0)
  {
    for(L=0;L<NQTL;L++) 
     {
      int T=0;
      for(PH=0;PH<NPHENO;PH++) T=T+GAMMA[PH][L];
        	if(T!=0) 
        	{
        		PD1[L]=0.0,PD2[L]=0.0;
        		for(I=0;I<NS;I++)
        		{
         			int K,KK=0;
        			for(K=0;K<NG;K++)
        			   if(QPROB[QCHR[0][L]][I][QLOC[0][L]][K]>0.99) 
        			   {
        				   IBD=K;
        				   KK=1;
        			   } 
        			if(KK==0)
                {
                MT_QTLgenotype(0,L,QCHR[0][L],QLOC[0][L],I); 
              	PD1[L]=PD1[L]+log(PDD1+1e-20);
                PD2[L]=PD2[L]+log(PDD2+1e-20);		    
                }
                if(IBD != GENO[0][I][L])
                {
              	for(PH=0;PH<NPHENO;PH++)	GVALUE[PH][I]=MT_GenotypeSampling(I,L,PH,IBD,0);
              	GENO[0][I][L]=IBD;
              	Coefficient(IBD);
              	for(K=0;K<NC;K++) COEF[0][I][L][K]=X[K];
                    for(PH=1;PH<NPHENO;PH++)
                       {
                          GENO[PH][I][L]=GENO[0][I][L];
                          for(K=0;K<NC;K++) COEF[PH][I][L][K]=COEF[0][I][L][K];     
                       }
              	
                }
  
        		}
        	}

     }
  }    	
}     
//**********************************************************************************
//UPDATING QTL POSITIONS
 

if(DiffLocation==1)
{
for(L=0;L<NQTL;L++) 
 {
  for(PH=0;PH<NPHENO;PH++) 
  {
      if(GAMMA[PH][L]!=0) 
    	{
		  int CAT[4],QLNEW,TEST=0;

  		double R=RANDOM();
	   	CAT[0]=(R>=0&&R<0.25),CAT[1]=(R>=0.25&&R<0.5),CAT[2]=(R>=0.5&&R<0.75),CAT[3]=(R>=0.75&&R<=1.0);
		  QLNEW=CAT[0]*(QLOC[PH][L]-2)+CAT[1]*(QLOC[PH][L]-1)+CAT[2]*(QLOC[PH][L]+1)+CAT[3]*(QLOC[PH][L]+2);
		  if(QLNEW<0) QLNEW=0;
		  if(QLNEW>=NGRID[QCHR[PH][L]]-1) QLNEW=NGRID[QCHR[PH][L]]-1;
		  if(QLNEW!=QLOC[PH][L])
		  {
			for(L0=0;L0<NQTL;L0++)
				if( QCHR[PH][L0]==QCHR[PH][L]&&L0!=L ) 
TEST=TEST+(fabs(GRID[QCHR[PH][L]][QLNEW]-GRID[QCHR[PH][L0]][QLOC[PH][L0]])<=DQQ[QCHR[PH][L]]);
// The above tests if proposed QTLNEW is not within the DQQ distance of 
// any existing QTLs in the model.	  
			if(TEST==0)   MT_QTLPOSITION(L,QLNEW,PH);
		  }
	   }
	}   
 } 
}	
  

if(DiffLocation==0)
{
for(L=0;L<NQTL;L++) 
 {
 	int T=0;
  for(PH=0;PH<NPHENO;PH++) T=T+GAMMA[PH][L];

      if(T!=0) 
    	{
		  int CAT[4],PH1,QLNEW,TEST=0;

  		double R=RANDOM();
	   	CAT[0]=(R>=0&&R<0.25),CAT[1]=(R>=0.25&&R<0.5),CAT[2]=(R>=0.5&&R<0.75),CAT[3]=(R>=0.75&&R<=1.0);
		  QLNEW=CAT[0]*(QLOC[0][L]-2)+CAT[1]*(QLOC[0][L]-1)+CAT[2]*(QLOC[0][L]+1)+CAT[3]*(QLOC[0][L]+2);
		  if(QLNEW<0) QLNEW=0;
		  if(QLNEW>=NGRID[QCHR[0][L]]-1) QLNEW=NGRID[QCHR[0][L]]-1;
		  if(QLNEW!=QLOC[0][L])
		  {
			for(L0=0;L0<NQTL;L0++)
				if( QCHR[0][L0]==QCHR[0][L]&&L0!=L ) 
    	TEST=TEST+(fabs(GRID[QCHR[0][L]][QLNEW]-GRID[QCHR[0][L0]][QLOC[0][L0]])<=DQQ[QCHR[0][L]]);
// The above tests if proposed QTLNEW is not within the DQQ distance of 
// any existing QTLs in the model.	  

// Have to test this in the future: QTLPOSITION_Samelocation might be the reason
//for anomalies in TMV: and QTLPOSITION_SameLocation is not coded for genoupdate
if(QTLLOC==0)			if(TEST==0)   QTLPOSITION_SameLocation(L,QLNEW);
if(QTLLOC==1)			if(TEST==0)   MT_QTLPOSITION(L,QLNEW,0);
		  }
         for(PH1=1;PH1<NPHENO;PH1++)
           {
               QCHR[PH1][L]=QCHR[0][L];
               QLOC[PH1][L]=QLOC[0][L];
               for(I=0;I<NS;I++)
               {
         				 GENO[PH1][I][L]=GENO[0][I][L];
                 for(K=0;K<NC;K++) COEF[PH1][I][L][K]=COEF[0][I][L][K];   
                }        
           }
  
//      if(T==NPHENO) PH=NPHENO-1;

	   }
	   
 } 
}	
      
//**********************************************************
// UPDATE MARGINAL EFFECT INDICATORS

for(L=0;L<NQTL;L++)
{
	int T,PH1,TT,K1;

	if(GROUP==0)
	{
		for(K=0;K<NC;K++)
		{
			T=1,TT=0;
			if(GIBBS==1) 
			{  
      if(DiffLocation==0) 
          {
          for(PH1=0;PH1<NPHENO;PH1++) TT=TT+GAMMA[PH1][L];        
                if(TT==0) 
                {
                T=MT_SamplingOnePosition(L,0);
                if(T!=0)   
                  for(PH1=1;PH1<NPHENO;PH1++)
                     {
                     QCHR[PH1][L]=QCHR[0][L];
                     QLOC[PH1][L]=QLOC[0][L];
                       for(I=0;I<NS;I++)
                				 {
                          GENO[PH1][I][L]=GENO[0][I][L];
                          for(K1=0;K1<NC;K1++) COEF[PH1][I][L][K1]=COEF[0][I][L][K1];     
                          }
                     }
                }   
          }

      			for(PH=0;PH<NPHENO;PH++)			
      			  {
             	if(DiffLocation==1)
             	   { 
                  T=1;
                  if(GAMMA[PH][L]==0) T=MT_SamplingOnePosition(L,PH);
                 }
                 
        				if(T!=0) MT_MainEffectIndicator_GROUP0(L,K,PH);
      
      			    if(MULTIPLE==1) 
                  {
                  for(PH1=0;PH1<NPHENO;PH1++) GAMMA[PH1][L]=GAMMA[PH][L];
                  PH=NPHENO-1;
                  }
        			}	
  				
  			
			}
	if(GIBBS==0)
			{
			int TTT[NPHENO],sumTTT=0,K1;
			for(PH=0;PH<NPHENO;PH++)			
			  {
    				double R=RANDOM();
    				TTT[PH]=0;
    				if( (MAIN[PH][L][K]==0&&R<=W_MAIN)||(MAIN[PH][L][K]!=0&&R>W_MAIN) ) 
    				{
             	if(DiffLocation==1)
             	  { T=1;
                  if(GAMMA[PH][L]==0) T=MT_SamplingOnePosition(L,PH);
              				if(T!=0) 
                				{
                					if(MAIN[PH][L][K]==0&&R<=W_MAIN)
                					{
                						double T0=0,U,U0;
                						for(J=0;J<NU;J++)
                						{
                							ANORMAL(&U,&U0);
                							T0=T0+U*U;
                						}
                						TAU=(NU-2)*H*VP[PH]/(NU*CC[K]);
                						VMAIN[PH][L][K]=NU*TAU/T0;
                					}
        				    MT_MainEffectIndicator_GROUP0(L,K,PH);
                				}
          			}
              else TTT[PH]=1,sumTTT=sumTTT+TTT[PH];  	
    				}
        }
         
           if(DiffLocation==0 && sumTTT!=0) 
            {
            for(PH1=0;PH1<NPHENO;PH1++) TT=TT+GAMMA[PH1][L];        
                  if(TT==0) 
                  {
                  T=MT_SamplingOnePosition(L,0);
                  if(T!=0)   
                    for(PH1=1;PH1<NPHENO;PH1++)
                       {
                       QCHR[PH1][L]=QCHR[0][L];
                       QLOC[PH1][L]=QLOC[0][L];
                       for(I=0;I<NS;I++)
                				 {
                          GENO[PH1][I][L]=GENO[0][I][L];
                          for(K1=0;K1<NC;K1++) COEF[PH1][I][L][K1]=COEF[0][I][L][K1];     
                          }
                       }
                  }   
            for(PH=0;PH<NPHENO;PH++) 
              {
                if(TTT[PH]==1 && T!=0) MT_MainEffectIndicator_GROUP0(L,K,PH);
              	    if(MULTIPLE==1) 
                      {
                      for(PH1=0;PH1<NPHENO;PH1++) GAMMA[PH1][L]=GAMMA[PH][L];
                      PH=NPHENO-1;
                      }
              }
                  
            }

			}

    
		}
	} 
  /*
	if(GROUP==1)
	{
		T=1;
		if(GIBBS==1) 
		{
			if(GAMMA[L]==0) T=MT_SamplingOnePosition(L);
			if(T!=0) MT_MainEffectIndicator_GROUP1(L);
		}
		if(GIBBS==0)
		{
			double R=RANDOM();
			if( (GAMMA_MAIN[L]==0&&R<=W_MAIN)||(GAMMA_MAIN[L]!=0&&R>W_MAIN) ) 
			{
				if(GAMMA[L]==0) T=MT_SamplingOnePosition(L);
				if(T!=0) MT_MainEffectIndicator_GROUP1(L);
			}
		}
	} */

 }     

  
    

//**********************************************************
// UPDATE TWO ORDER EPISTATIC EFFECT INDICATORS
		
if(EPISTASIS==1)
{
	int T,PH1;
	
	for(L1=0;L1<NQTL-1;L1++)
		for(L2=L1+1;L2<NQTL;L2++)
		{
			if(GROUP==0)
			{
				for(K1=0;K1<NC;K1++)
					for(K2=0;K2<NC;K2++)
					{
				/*	if(DEPENDENCE==1)
						{
							if(MAIN[L1][K1]!=0&&MAIN[L2][K2]!=0) W_EPISTASIS=C[0];
							if( (MAIN[L1][K1]!=0&&MAIN[L2][K2]==0)||(MAIN[L1][K1]==0&&MAIN[L2][K2]!=0) ) W_EPISTASIS=C[1];
							if(MAIN[L1][K1]==0&&MAIN[L2][K2]==0) W_EPISTASIS=C[2];
						}   */

						T=1;
		if(GIBBS==1) 
   	{
      	if(W_EPISTASIS!=0)
        {
            if(DiffLocation==0)
              {       int T1=0,T2=0;
      								for(PH=0;PH<NPHENO;PH++) 
                        {
                          T1=T1+GAMMA[PH][L1];
                          T2=T2+GAMMA[PH][L2];
                        }  
                      if(T1==0 && T2!=0) T=MT_SamplingOnePosition(L1,0);
      								if(T1!=0 && T2==0) T=MT_SamplingOnePosition(L2,0);
      								if(T1==0 && T2==0) 
      								{
      									T=MT_SamplingOnePosition(L1,0);
      									if(T!=0) T=MT_SamplingOnePosition(L2,0);
				if(CHRQTL[QCHR[0][L1]]+2>CHR_NQTL[QCHR[0][L1]]&&QCHR[0][L1]==QCHR[0][L2]
&&(GRID[QCHR[0][L1]][QLOC[0][L1]]-GRID[QCHR[0][L2]][QLOC[0][L2]])<=DQQ[QCHR[0][L1]]) T=0;
      								}
                      if(T!=0)
                        for(PH1=1;PH1<NPHENO;PH1++)
                           {
                           QCHR[PH1][L1]=QCHR[0][L1];
                           QLOC[PH1][L1]=QLOC[0][L1];
                           QCHR[PH1][L2]=QCHR[0][L2];
                           QLOC[PH1][L2]=QLOC[0][L2];
                           for(I=0;I<NS;I++)
                              {
                              GENO[PH1][I][L1]=GENO[0][I][L1];
                              GENO[PH1][I][L2]=GENO[0][I][L2];
                              COEF[PH1][I][L1][K1]=COEF[0][I][L1][K1];
                              COEF[PH1][I][L2][K2]=COEF[0][I][L2][K2];
                              }
                           }  
      								
							}
          }
        if(DiffLocation==0)
          for(PH=0;PH<NPHENO;PH++)
          {			if(T!=0) MT_EpistasisIndicator_GROUP0(L1,L2,K1,K2,PH);
      			    if(MULTIPLE==1)
                  {
                  for(PH1=0;PH1<NPHENO;PH1++){
                   GAMMA[PH1][L1]=GAMMA[PH][L1];
                   GAMMA[PH1][L2]=GAMMA[PH][L2];
                   }
                  PH=NPHENO-1;
                  }
          }
          if(DiffLocation==1)
  					for(PH=0;PH<NPHENO;PH++)
            {		
            	if(W_EPISTASIS!=0)
              {
								if(GAMMA[PH][L1]==0 && GAMMA[PH][L2]!=0) T=MT_SamplingOnePosition(L1,PH);
								if(GAMMA[PH][L1]!=0 && GAMMA[PH][L2]==0) T=MT_SamplingOnePosition(L2,PH);
								if(GAMMA[PH][L1]==0 && GAMMA[PH][L2]==0)
								{
									T=MT_SamplingOnePosition(L1,PH);
									if(T!=0) T=MT_SamplingOnePosition(L2,PH);
									if(CHRQTL[QCHR[PH][L1]]+2>CHR_NQTL[QCHR[PH][L1]]&&QCHR[PH][L1]==QCHR[PH][L2]
					&&(GRID[QCHR[PH][L1]][QLOC[PH][L1]]-GRID[QCHR[PH][L2]][QLOC[PH][L2]])<=DQQ[QCHR[PH][L1]]) T=0;

								}
						  }		
              	if(T!=0) MT_EpistasisIndicator_GROUP0(L1,L2,K1,K2,PH);

							
					  }		
		      			

       } 
	if(GIBBS==0)
	{
      			int TTT[NPHENO],sumTTT=0;
      			for(PH=0;PH<NPHENO;PH++)			
      			{	TTT[PH]=0;
 							double R=RANDOM();
if((EPISTATIC[PH][L1][L2][K1][K2]==0&&R<=W_EPISTASIS)||(EPISTATIC[PH][L1][L2][K1][K2]!=0&&R>W_EPISTASIS)) 
							{
								if(W_EPISTASIS!=0)
								{
								 if(DiffLocation==1)
								 	{ T=1;
                  if(GAMMA[PH][L1]==0&&GAMMA[PH][L2]!=0) T=MT_SamplingOnePosition(L1,PH);
									if(GAMMA[PH][L1]!=0&&GAMMA[PH][L2]==0) T=MT_SamplingOnePosition(L2,PH);
									if(GAMMA[PH][L1]==0&&GAMMA[PH][L2]==0) 
									{ T=MT_SamplingOnePosition(L1,PH);
										if(T!=0) T=MT_SamplingOnePosition(L2,PH);
if(CHRQTL[QCHR[PH][L1]]+2>CHR_NQTL[QCHR[PH][L1]]&&QCHR[PH][L1]==QCHR[PH][L2]
&&(GRID[QCHR[PH][L1]][QLOC[PH][L1]]-GRID[QCHR[PH][L2]][QLOC[PH][L2]])<=DQQ[QCHR[PH][L1]]) T=0;                }
						if(T!=0) 
							{
								if(EPISTATIC[PH][L1][L2][K1][K2]==0&&R<=W_EPISTASIS)
								{
									double T0=0,U,U0;
									for(J=0;J<NU;J++)
									{
										ANORMAL(&U,&U0);
										T0=T0+U*U;
									}
									TAU=(NU-2)*H*VP[PH]/(NU*CC[K1]*CC[K2]);
									VEPISTASIS[PH][L1][L2][K1][K2]=NU*TAU/T0;
								}
							  MT_EpistasisIndicator_GROUP0(L1,L2,K1,K2,PH);
							}

									} else TTT[PH]=1,sumTTT=sumTTT+TTT[PH];  	
                }
							}		
						}
            
           if(DiffLocation==0 && sumTTT!=0) 
            {
              int T1=0,T2=0;
      								for(PH1=0;PH1<NPHENO;PH1++) 
                        {
                          T1=T1+GAMMA[PH1][L1];
                          T2=T2+GAMMA[PH1][L2];
                        }  
                      if(T1==0 && T2!=0) T=MT_SamplingOnePosition(L1,0);
      								if(T1!=0 && T2==0) T=MT_SamplingOnePosition(L2,0);
      								if(T1==0 && T2==0) 
      								{
      									T=MT_SamplingOnePosition(L1,0);
      									if(T!=0) T=MT_SamplingOnePosition(L2,0);
				if(CHRQTL[QCHR[0][L1]]+2>CHR_NQTL[QCHR[0][L1]]&&QCHR[0][L1]==QCHR[0][L2]
&&(GRID[QCHR[0][L1]][QLOC[0][L1]]-GRID[QCHR[0][L2]][QLOC[0][L2]])<=DQQ[QCHR[0][L1]]) T=0;
      								}
                      if(T!=0)
                        for(PH1=1;PH1<NPHENO;PH1++)
                           {
                           QCHR[PH1][L1]=QCHR[0][L1];
                           QLOC[PH1][L1]=QLOC[0][L1];
                           QCHR[PH1][L2]=QCHR[0][L2];
                           QLOC[PH1][L2]=QLOC[0][L2];
                           for(I=0;I<NS;I++)
                              {
                              GENO[PH1][I][L1]=GENO[0][I][L1];
                              GENO[PH1][I][L2]=GENO[0][I][L2];
                              COEF[PH1][I][L1][K1]=COEF[0][I][L1][K1];
                              COEF[PH1][I][L2][K2]=COEF[0][I][L2][K2];
                              }
                           }  
            for(PH=0;PH<NPHENO;PH++) 
              {
                if(TTT[PH]==1 && T!=0) MT_EpistasisIndicator_GROUP0(L1,L2,K1,K2,PH);
              	    if(MULTIPLE==1) 
                      {
                      for(PH1=0;PH1<NPHENO;PH1++){ 
                      GAMMA[PH1][L1]=GAMMA[PH][L1];
                      GAMMA[PH1][L2]=GAMMA[PH][L2];
                      }
                      PH=NPHENO-1;
                      }
              }
                  
            }
            		
								
             
      }  
    }
  } 

/*			if(GROUP==1)
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
						if(GAMMA[L1]==0&&GAMMA[L2]!=0) T=MT_SamplingOnePosition(L1);
						if(GAMMA[L1]!=0&&GAMMA[L2]==0) T=MT_SamplingOnePosition(L2);
						if(GAMMA[L1]==0&&GAMMA[L2]==0)
						{
							T=MT_SamplingOnePosition(L1);
							if(T!=0) T=MT_SamplingOnePosition(L2);
							if(CHRQTL[QCHR[L1]]+2>CHR_NQTL[QCHR[L1]]&&QCHR[L1]==QCHR[L2]
								&&(GRID[QCHR[L1]][QLOC[L1]]-GRID[QCHR[L2]][QLOC[L2]])<=DQQ[QCHR[L1]]) T=0;
						}
					}
					if(T!=0) MT_EpistasisIndicator_GROUP1(L1,L2);
				}
				if(GIBBS==0)
				{
					double R=RANDOM();
					if( (GAMMA_EPISTASIS[L1][L2]==0&&R<=W_EPISTASIS)||(GAMMA_EPISTASIS[L1][L2]!=0&&R>W_EPISTASIS) )
					{
						if(W_EPISTASIS!=0)
						{
							if(GAMMA[L1]==0&&GAMMA[L2]!=0) T=MT_SamplingOnePosition(L1);
							if(GAMMA[L1]!=0&&GAMMA[L2]==0) T=MT_SamplingOnePosition(L2);
							if(GAMMA[L1]==0&&GAMMA[L2]==0) 
							{
								T=MT_SamplingOnePosition(L1);
								if(T!=0) T=MT_SamplingOnePosition(L2);
								if(CHRQTL[QCHR[L1]]+2>CHR_NQTL[QCHR[L1]]&&QCHR[L1]==QCHR[L2]
									&&(GRID[QCHR[L1]][QLOC[L1]]-GRID[QCHR[L2]][QLOC[L2]])<=DQQ[QCHR[L1]]) T=0;
							}
						}
						if(T!=0) MT_EpistasisIndicator_GROUP1(L1,L2);
					}
				}
			}  */

		} // L1 & L2 loop end here
} //if(epistasis) end here 

 
//************************************************************
//update g by e fixed effects INDICATORS

if(GBYE==1)
{ 
for(L1=0;L1<NFIXCOVA;L1++)
 if(GBYE_FIX_INDEX[L1]==1)
	for(PH=0;PH<NPHENO;PH++)
	 for(L2=0;L2<NQTL;L2++)		
		if(GAMMA[PH][L2]!=0)
		{	
			if(GROUP==0)
			{
				for(K=0;K<NC;K++)
				{
					if(GIBBS==1) MT_GBYE_FIX_Indicator_GROUP0(L1,L2,K,PH);
					if(GIBBS==0)
					{
						double R=RANDOM();
				if((GBYE_FIX[PH][L1][L2][K]==0&&R<=W_GBYE)||(GBYE_FIX[PH][L1][L2][K]!=0&&R>W_GBYE))
                      						{
                      							if(GBYE_FIX[PH][L1][L2][K]==0&&R<=W_GBYE)
                      							{
                      								double T0=0,U,U0;
                      								for(J=0;J<NU;J++)
                      								{
                      									ANORMAL(&U,&U0);
                      									T0=T0+U*U;
                      								}
                      								TAU=(NU-2)*H*VP[PH]/(NU*V_FIX[L1]*CC[K]);
                      								V_GBYE_FIX[PH][L1][L2][K]=NU*TAU/T0;
                      							}
                      	
                							MT_GBYE_FIX_Indicator_GROUP0(L1,L2,K,PH);
                      						}
							
					}
				}
     	    if(MULTIPLE==1) 
              {
              for(int PH1=0;PH1<NPHENO;PH1++) GAMMA[PH1][L2]=GAMMA[PH][L2];
              PH=NPHENO-1;
              }
			}

		/*	if(GROUP==1)
			{
				if(GIBBS==1) MT_GBYE_FIX_Indicator_GROUP1(L1,L2);

				if(GIBBS==0)
				{
					double R=RANDOM();
					if( (GAMMA_GBYE[L1][L2]==0&&R<=W_GBYE)||(GAMMA_GBYE[L1][L2]!=0&&R>W_GBYE) )
						MT_GBYE_FIX_Indicator_GROUP1(L1,L2);
				}
			}*/

		}
}

//**********************************************************
//update the threshold values
/*
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
*/
//***************************************************************
   
}    //ITER1 end here
   
//***************************************************************
//CALCULATE THE NUMBER OF QTL
/*
Rprintf("\n");
for(PH=0;PH<NPHENO;PH++)
{
Rprintf("\nGAMMA[%d]",PH+1);           
  for(L=0;L<NQTL;L++)
  Rprintf("\t%d",GAMMA[PH][L]);             
}
Rprintf("\n");
for(PH=0;PH<NPHENO;PH++)
{
Rprintf("\nGMAIN[%d]",PH+1);           
  for(L=0;L<NQTL;L++)
  Rprintf("\t%d",GAMMA_MAIN[PH][L]);             
}
Rprintf("\n");
for(PH=0;PH<NPHENO;PH++)
{
for(L1=0;L1<NFIXCOVA;L1++)
{
Rprintf("\nGBYE[%d,%d]",PH+1,L1+1);           
  for(L=0;L<NQTL;L++)
  Rprintf("\t%f",GAMMA_GBYE[PH][L1][L]);             
}  
} */




int QTL_INCLUDED[NPHENO];
for(PH=0;PH<NPHENO;PH++)
{ QTL_INCLUDED[PH]=0;
  for(L=0;L<NQTL;L++) QTL_INCLUDED[PH]=QTL_INCLUDED[PH]+(GAMMA[PH][L]==1);  
}
//*************************************************************** 
//SAVE THE RESULT 
 
if(ITER*ITER1>=NBURNIN)
{

if((ITER!=0)&&(ITER%10==0)&(VERBOSE>0))  
{
	Rprintf("%d",ITER);
	Rprintf("\n");
}
int PH1;
 if(MULTIPLE>=1)   INVERSE(SIGMA,NPHENO,sigma);


/*//CALCULATE DEVIANCE
double DEV=0;  
for(I=0;I<NS;I++) DEV=DEV+log(2*3.14159265358*VE)+pow(Y[I]-AMU-GVALUE[I],2)/VE;
*/
// save to "iterdiag"
for(PH1=0;PH1<NPHENO;PH1++)
{
fprintf(File1,"\n");
fprintf(File1,"%d\t",ITER+1);
fprintf(File1,"%d\t",PH1+1); /* yandell add */ 
fprintf(File1,"%d\t",QTL_INCLUDED[PH1]);
fprintf(File1,"%f\t",AMU[PH1]);
if(MULTIPLE>=1) VE=sigma[PH1][PH1];
fprintf(File1,"%f\t",VE);


// save to "covariates"
if(ENV_FACTOR==1)
{
	fprintf(File2,"\n");
  fprintf(File2,"%d\t",PH1+1); /* yandell add */ 		

	for(L=0;L<NFIXCOVA;L++) fprintf(File2,"%f\t",FIX[PH1][L]);
	for(L=0;L<NRANCOVA;L++) fprintf(File2,"%f\t",VRAN[PH1][L]);
} 
if(CATEGORY==3)
  for(L=2;L<=CN-2;L++) fprintf(File2,"%f\t",CUTPOINT[L]);
      
// save to "mainloci"
double VAR1[NG];
for(K=0;K<NC;K++) VAR1[K]=0;
for(L=0;L<NQTL;L++) 
	if(GAMMA[PH1][L]!=0) 
	{
		fprintf(File3,"\n");
		fprintf(File3,"%d\t",ITER+1);
    fprintf(File3,"%d\t",PH1+1); /* yandell add */ 		
        fprintf(File3,"%d\t",QTL_INCLUDED[PH1]);
        fprintf(File3,"%d\t",QCHR[PH1][L]+1);
//        fprintf(File3,"%f\t",GRID[QCHR[L]][QLOC[L]]);
		fprintf(File3,"%d\t",QLOC[PH1][L]);
		for(K=0;K<NC;K++) fprintf(File3,"%f\t",MAIN[PH1][L][K]);
		for(K=0;K<NC;K++)
		{
			double SS1=0,SS2=0,SS=0;
			if(MAIN[PH1][L][K]!=0)
			{
				for(I=0;I<NS;I++)
				{
					SS1=SS1+pow(COEF[PH1][I][L][K],2);
					SS2=SS2+COEF[PH1][I][L][K];
				}
				SS=1.0/(NS-1)*(SS1-NS*pow(SS2/NS,2))*pow(MAIN[PH1][L][K],2);
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
					if(EPISTATIC[PH1][L1][L2][K1][K2]!=0) N_EPIS=N_EPIS+1;

	for(L1=0;L1<NQTL;L1++)
		for(L2=L1+1;L2<NQTL;L2++) 
			if(GAMMA[PH1][L1]!=0&&GAMMA[PH1][L2]!=0&&GAMMA_EPISTASIS[PH1][L1][L2]!=0) 
			{
				fprintf(File4,"\n");
				fprintf(File4,"%d\t",ITER+1);
        fprintf(File4,"%d\t",PH1+1); /* yandell add */ 				
				fprintf(File4,"%d\t",N_EPIS);
				fprintf(File4,"%d\t",QCHR[PH1][L1]+1);
//				fprintf(File4,"%f\t",GRID[QCHR[L1]][QLOC[L1]]);
				fprintf(File4,"%d\t",QLOC[PH1][L1]);
				fprintf(File4,"%d\t",QCHR[PH1][L2]+1);
//				fprintf(File4,"%f\t",GRID[QCHR[L2]][QLOC[L2]]);
				fprintf(File4,"%d\t",QLOC[PH1][L2]);
				for(K1=0;K1<NC;K1++)
					for(K2=0;K2<NC;K2++) fprintf(File4,"%f\t",EPISTATIC[PH1][L1][L2][K1][K2]);

				for(K1=0;K1<NC;K1++)
				for(K2=0;K2<NC;K2++)
				{
					double SS1=0,SS2=0,SS=0;
					if(EPISTATIC[PH1][L1][L2][K1][K2]!=0)
					{
						for(I=0;I<NS;I++)
						{
							SS1=SS1+pow(COEF[PH1][I][L1][K1]*COEF[PH1][I][L2][K2],2);
							SS2=SS2+COEF[PH1][I][L1][K1]*COEF[PH1][I][L2][K2];
						}
						SS=1.0/(NS-1)*(SS1-NS*pow(SS2/NS,2))*pow(EPISTATIC[PH1][L1][L2][K1][K2],2);
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
				if(GBYE_FIX[PH1][L1][L2][K]!=0) N_GBYE=N_GBYE+1;
	
	for(L1=0;L1<NFIXCOVA;L1++)
		 if(GBYE_FIX_INDEX[L1]==1)
			for(L2=0;L2<NQTL;L2++)
				if(GAMMA[PH1][L2]!=0.0)
					if(GAMMA_GBYE[PH1][L1][L2]!=0.0)
					{
						fprintf(File5,"\n");
						fprintf(File5,"%d\t",ITER+1);
            fprintf(File5,"%d\t",PH1+1); /* yandell add */ 		
						fprintf(File5,"%d\t",N_GBYE);
						fprintf(File5,"%d\t",L1+1);
						fprintf(File5,"%d\t",QCHR[PH1][L2]+1);
					//	fprintf(File5,"%f\t",GRID[QCHR[L2]][QLOC[L2]]);
						fprintf(File5,"%d\t",QLOC[PH1][L2]);
						for(K=0;K<NC;K++) fprintf(File5,"%f\t",GBYE_FIX[PH1][L1][L2][K]);

						for(K=0;K<NC;K++)
						{
							double SS1=0,SS2=0,SS=0;
							if(GBYE_FIX[PH1][L1][L2][K]!=0)
							{
								for(I=0;I<NS;I++)
								{
									SS1=SS1+pow(COEF_FIX[I][L1]*COEF[PH1][I][L2][K],2);
									SS2=SS2+COEF_FIX[I][L1]*COEF[PH1][I][L2][K];
								}
								SS=1.0/(NS-1)*(SS1-NS*pow(SS2/NS,2))*pow(GBYE_FIX[PH1][L1][L2][K],2);
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
if(ENV_FACTOR==1) fprintf(File1,"%f\t",VP[PH1]-VAR-VE);
fprintf(File1,"%f\t",VAR);
/*
//save to "deviance"
fprintf(File6,"%d\t",ITER);
fprintf(File6,"%lf\t",DEV);
fprintf(File6,"\n");

}
 */
     
 
if(MULTIPLE>=1 & PH1 == 0) // only print sigma once per iteration.
 { 
    for(I=0;I<NPHENO;I++)
          for(J=0;J<NPHENO;J++)           fprintf(File7,"%f\t",sigma[I][J]);

        fprintf(File7,"\n");   
  }    
fflush(NULL);
 }
}

  
} //ITER end here


fclose(File1);  
fclose(File2); 
fclose(File3); 
fclose(File4);
fclose(File5);
fclose(File6);
fclose(File7);
  

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

      
