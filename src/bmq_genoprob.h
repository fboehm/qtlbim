
/* 
 * This header file has global data structures that are used in Yi's C code genoprob.c. 
 * This code calculates the grid points and genotypic probabilities.
 * It also has function declarations required for setting the analysis parameters
 * and functions in genoprob.c
 * 
 */  
//********************************************************************************
//********************************************************************************

#ifndef BMQ_GENOPROB_H
#define BMQ_GENOPROB_H


// global parameters

extern int CROSS;         // cross type 0:RILs, 1:BC, 2:F2
extern int NG;            // number of genotypes
extern int NS;            // number of individuals 
extern int CHL;           // number of grids
extern int M;             // number of marker   

//********************************************************************************

extern int HALDANE;      // 1: Haldane distance; 0: Kosambi distance
extern double CSTEP;     // length of the grid
extern double *RECM;     // marker map
extern double LCHR;      // chr length
extern int **MK;         // marker genotype
extern int KCHR;         // number of markers

extern double ***GenoProb;
extern double *CHR_GRID;

//***************************************************************************************

void R_GenoGrid(int *nind, int *nmar, double *map, double *chrlen, int *max_chrlen,
	            int *genodata, int *type, int *mapfunc, double *step,double *chrloci,double *chrprob);


void MULTIPOINT_BC(double *MMK,int MKQ,double *PAA,double *PQ);

void MULTIPOINT_F2(double *MMK,int MKQ,double *PAA,double *PQ);

void getGenoGrid(double *chrloci,double *chrprob);

//***************************************************************************************

void R_OutputManager(char **iterFile,char **covFile,char **mainFile,char **pairFile,char **gbyeFile);

#endif // BMQ_GENOPROB_H

