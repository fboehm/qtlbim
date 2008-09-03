// QTLBIM - QTL Bayesian Interval Mapping
// Functions that are called from the R environment
//********************************************************************

#ifndef R_INTERFACES_H
#define R_INTERFACES_H


void RSingleTraitMCMCSetup(int *nind,int *nchr,int *ngen, int *nloci,double *loci,double *prob,
					  double *yvalue,int *traittype,int *ncategory,
					  int *iter,int *thin,int *burnin,int *genoupdate,
                      int *epis,int *emainqtl,int *eqtl,int *mnqtl,
					  double *interval,int *chrnqtl,
					  int *envi,int *qtl_envi,int *nrancov,int *nfixcov,int *intcov,double *rancoef,double *fixcoef,int *nran,
					  int *depen,double *prop,int *contrast,double *censor_lo,double *censor_hi,
					  int *seed,int *verbose);


void RMultipleTraitsMCMCSetup(int *nind,int *nchr,int *ngen,int *npheno,int *nloci,double *loci,double *prob,
					  double *yvalue,int *multipletrait,int *traittype,int *ncategory,
					  int *iter,int *thin,int *burnin,int *algorithm,int *genoupdate,
					  int *epis,int *emainqtl,int *eqtl,int *mnqtl,
					  double *interval,int *chrnqtl,
					  int *envi,int *qtl_envi,int *nrancov,int *nfixcov,int *intcov,double *rancoef,double *fixcoef,int *nran,
					  int *depen,double *prop,
					  int *seed,int *verbose,
					  int *diffloc, int *qtlloc);  


void RBayesianAnovaSetup( int *nind,int *nchr,int *ngen, int *nloci,double *loci,double *prob,
					 double *yvalue,int *traittype,int *ncategory,
					 int *iter,int *thin,int *burnin,
                     int *genoupdate,int *posupdate,
                     int *epis,int *mnqtl,
					 double *interval,int *chrnqtl,
					 int *envi,int *qtl_envi,int *nrancov,int *nfixcov,int *intcov,double *rancoef,double *fixcoef,int *nran,
				     int *contrast,
					 int *qchr,int *qloc,int *gamma_main,int *gamma_epis,double *gamma_gbye,
					 int *seed,int *verbose );


void ROutputManager(char **iterFile,char **covFile,char **mainFile,char **pairFile,char **gbyeFile,char **devFile,char **sigmaFile);

void R_CheckUserInterrupt(void);

#endif // R_INTERFACES_H
