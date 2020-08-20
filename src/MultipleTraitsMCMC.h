// QTLBIM - QTL Bayesian Interval Mapping
// Functions relevant to MCMC Model Selection for multiple trait analyses
//********************************************************************

#ifndef MULTIPLE_TRAITS_MCMC_H
#define MULTIPLE_TRAITS_MCMC_H

/// Main MCMC analysis driver for multiple traits
void multipleTraitsMCMC();

/// @note: Temporary redeclaration here so that the compiler does not throw a warning. Definition in SingleTraitsMCMCSamplingRoutines.c
void Coefficient(int GENOTYPE);

#endif // MULTIPLE_TRAITS_MCMC_H
