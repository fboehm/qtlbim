// QTLBIM - QTL Bayesian Interval Mapping
// Defining all the global variables used in the code
// This should only be a temporary file, a stepping stone in the process of moving away from global variables
//********************************************************************

#include "GlobalVars.h"

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
int UPDATEPOS;        // 1: update QTL positions; 0: doesn't update QTL positions

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

double  **VMAIN;         // prior variance of main effects
double ****VEPISTASIS;    // prior variance of epistatic effects
double ***V_GBYE_FIX;    // prior variance of g-by-e effects

double  *VMAIN1;         // prior variance of main effects
double **VEPISTASIS1;    // prior variance of epistatic effects
double **V_GBYE_FIX1;    // prior variance of g-by-e effects

//*********************************************************************************
// genetic model parameters

double AMU;             // overall mean
double *VE;              // residual variance
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

//**********************************************************************************

double PDD1, PDD2;
double *PD1, *PD2;

int IBD;

double  *X;

//**********************************************************************************

char iterfile[100];
char pairfile[100];
char mainfile[100];
char gbyefile[100];
char covfile[100];
char devfile[100];

//**********************************************************************************

double *CENSOR_LO;
double *CENSOR_HI;   
