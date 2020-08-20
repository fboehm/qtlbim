#include "GlobalVars_SingleTrait.h"

//***************************************************************************************
// group qbOp variable definitions
/// @{


/// coefficients of QTL main effectsf
double ***COEF;
/// QTL position indicators, position is GRID[QCHR[L]][QLOC[L]]
int *QLOC;
/// chromosomes that QTL locate
int *QCHR;

    /// group qbOpCovariate variable definitions
    /// @{
    /// effects of fixed covariates
    double  *FIX;
    /// effects of random covariates
    double  **RAN;
    /// variances of random covariates
    double  *VRAN;
    /// interactions of QTL main effects and fixed covariates
    double ***GBYE_FIX;
    /// @}

    /// group qbOpIterdiag variable definitions
    /// @{
    /// overall mean
    double AMU;
    /// QTL indicators
    int *GAMMA;
    /// residual variance
    double *VE;
    /// @}

    /// group qbOpMainloci variable definitions
    /// @{
    /// main effects
    double  **MAIN;
    /// main effects indicators
    int *GAMMA_MAIN;
    /// prior variance of main effects
    double  **VMAIN;
    /// prior variance of main effects
    double  *VMAIN1;
    /// @}

    /// group qbOpPairloci variable definitions
    /// @{
    /// epistatic effects
    double ****EPISTATIC;
    /// epistatic effects indicators
    int **GAMMA_EPISTASIS;
    /// prior variance of epistatic effects
    double ****VEPISTASIS;
    /// prior variance of epistatic effects
    double **VEPISTASIS1;
    /// @}

    /// group qbOpGbyE variable definitions
    /// @{
    /// g-by-e indicators
    double  **GAMMA_GBYE;
    /// prior variance of g-by-e effects
    double ***V_GBYE_FIX;
    /// prior variance of g-by-e effects
    double **V_GBYE_FIX1;
    /// @}

/// @}


/// QTL genotypes
int **GENO;
/// phenotypic data
double *Y;
/// genotypic valyes
double *GVALUE;
