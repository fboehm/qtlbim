
/* 
 * This header file has global data structures that are used in Yi's C code genoprob.c. 
 * This code calculates the grid points and genotypic probabilities.
 * It also has function declarations required for setting the analysis parameters
 * and functions in genoprob.c
 * 
 */  
//********************************************************************************

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Random.h>
#include <stdio.h>
#include <time.h>

//********************************************************************************
// global parameters

int CROSS;         // cross type 0:RILs, 1:BC, 2:F2

int NG;            // number of genotypes

int NS;            // number of individuals 

int CHL;           // number of grids

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

//***************************************************************************************

void R_GenoGrid(int *nind, int *nmar, double *map, double *chrlen, int *max_chrlen,
	            int *genodata, int *type, int *mapfunc, double *step,double *chrloci,double *chrprob);


void MULTIPOINT_BC(double *MMK,int MKQ,double *PAA,double *PQ);

void MULTIPOINT_F2(double *MMK,int MKQ,double *PAA,double *PQ);

void getGenoGrid(double *chrloci,double *chrprob);

//***************************************************************************************

void R_OutputManager(char **iterFile,char **covFile,char **mainFile,char **pairFile,char **gbyeFile);
