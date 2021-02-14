/** 
@brief This file holds all the variables that are differently dimensioned for the multiple traits code.
 *  
 * @note None of the global variables in this file have been extern-ed. This is so that they do not conflict with
 * similarly named single-trait versions of these variables that have different dimensions.
 */ 
#ifndef GLOBAL_VARS_MULTIPLE_TRAIT_H
#define GLOBAL_VARS_MULTIPLE_TRAIT_H


/// 0=uses QTLPOSITION_samelocation 1= uses QTLPOSITION
extern int QTLLOC;
/// Whether different locations are considered for each trait for each indicator or same locations are considered
extern int DiffLocation;
/// Residual variance matrix for multiple triats
extern double **SIGMA;
/// Main Effect Variance scale for beta ~ N(0,cI) c=scale;
extern double SCALE;
//*******************************************************************************************


/// @addtogroup qbGeno
/// @{

/// QTL genotypes for the MT code
extern int ***GENO;

/// @}

/// @addtogroup qbPheno
/// @{

/// phenotypic data for the MT code
extern double **Y;
/// The number of phenotypes to be jointly analyzed
extern int NPHENO;
/// 
extern int MULTIPLE;

/// @}

/// genotypic valyes for the MT code
extern double **GVALUE;


//***************************************************************************************
/// @defgroup qbOp_MT Multiple Traits variants of the qb.op variables
/// @{


/// coefficients of QTL main effects
extern double ****COEF;
/// QTL position indicators, position is GRID[QCHR[L]][QLOC[L]]
extern int **QLOC;
/// chromosomes that QTL locate
extern int **QCHR;

    /// @defgroup qbOpCovariate_MT Multiple Traits variant of qb.op.covariate
    /// @{
    /// effects of fixed covariates
    extern double  **FIX;
    /// effects of random covariates
    extern double  ***RAN;
    /// variances of random covariates
    extern double  **VRAN;
    /// interactions of QTL main effects and fixed covariates
    extern double ****GBYE_FIX;
    /// @}

    /// @defgroup qbOpIterdiag_MT Multiple Traits variant of qb.op.iterdiag
    /// @{
    /// overall mean
    extern double *AMU;
    /// QTL indicators
    extern int **GAMMA;
    /// residual variance
    extern double VE;
    /// @}

    /// @defgroup qbOpMainloci_MT Multiple Traits variant of qb.op.mainloci
    /// @{
    /// main effects
    extern double  ***MAIN;
    /// main effects indicators
    extern int **GAMMA_MAIN;
    /// prior variance of main effects
    extern double  ***VMAIN;
    /// @}

    /// @defgroup qbOpPairloci_MT Multiple Traits variant of qb.op.pairloci
    /// @{
    /// epistatic effects
    extern double *****EPISTATIC;
    /// epistatic effects indicators
    extern int ***GAMMA_EPISTASIS;
    /// prior variance of epistatic effects
    extern double *****VEPISTASIS;
    /// @}

    /// @defgroup qbOpGbyE_MT Multiple Traits variant of qb.op.gbye
    /// @{
    /// g-by-e indicators
    extern double  ***GAMMA_GBYE;
    /// prior variance of g-by-e effects
    extern double ****V_GBYE_FIX;
    /// @}

/// @}

#endif // GLOBAL_VARS_MULTIPLE_TRAIT_H
