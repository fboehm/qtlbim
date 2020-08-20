#ifndef GLOBAL_VARS_SINGLE_TRAIT_H
#define GLOBAL_VARS_SINGLE_TRAIT_H

/// @addtogroup qbGeno
/// @{

/// QTL genotypes
/// @note: RV contests the inclusion of this variable into a Genotype data structure. It is the only var in Geno that gets extra dimensions for multiple traits. Does that indicate that it belongs elsewhere?
extern int **GENO;

/// @}

/// @addtogroup qbPheno
/// @{

/// phenotypic data
extern double *Y;

/// @}

/// genotypic valyes
extern double *GVALUE;


//***************************************************************************************
/// @defgroup qbOp Variables that should reside in the qp.op structure
/// @{


/// coefficients of QTL main effects
extern double ***COEF;
/// QTL position indicators, position is GRID[QCHR[L]][QLOC[L]]
extern int *QLOC;
/// chromosomes that QTL locate
extern int *QCHR;

    /// @defgroup qbOpCovariate Variables in qb.op.covariate
    /// @{
    /// effects of fixed covariates
    extern double  *FIX;
    /// effects of random covariates
    extern double  **RAN;
    /// variances of random covariates
    extern double  *VRAN;
    /// interactions of QTL main effects and fixed covariates
    extern double ***GBYE_FIX;
    /// @}

    /// @defgroup qbOpIterdiag Variables in qb.op.iterdiag
    /// @{
    /// overall mean
    extern double AMU;
    /// QTL indicators
    extern int *GAMMA;
    /// residual variance
    extern double *VE;
    /// @}

    /// @defgroup qbOpMainloci Variables in qb.op.mainloci
    /// @{
    /// main effects
    extern double  **MAIN;
    /// main effects indicators
    extern int *GAMMA_MAIN;
    /// prior variance of main effects
    extern double  **VMAIN;
    /// prior variance of main effects
    extern double  *VMAIN1;
    /// @}

    /// @defgroup qbOpPairloci Variables in qb.op.pairloci
    /// @{
    /// epistatic effects
    extern double ****EPISTATIC;
    /// epistatic effects indicators
    extern int **GAMMA_EPISTASIS;
    /// prior variance of epistatic effects
    extern double ****VEPISTASIS;
    /// prior variance of epistatic effects
    extern double **VEPISTASIS1;
    /// @}

    /// @defgroup qbOpGbyE Variables in qb.op.gbye
    /// @{
    /// g-by-e indicators
    extern double  **GAMMA_GBYE;
    /// prior variance of g-by-e effects
    extern double ***V_GBYE_FIX;
    /// prior variance of g-by-e effects
    extern double **V_GBYE_FIX1;
    /// @}

/// @}

//*************************************************************************************************************
#endif // GLOBAL_VARS_SINGLE_TRAIT_H
