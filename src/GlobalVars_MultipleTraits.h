/** 
@brief This file holds all the variables that are differently dimensioned for the multiple traits code.
 *  
 * @note None of the global variables in this file have been extern-ed. This is so that they do not conflict with
 * similarly named single-trait versions of these variables that have different dimensions.
 */ 
#ifndef GLOBAL_VARS_MULTIPLE_TRAIT_H
#define GLOBAL_VARS_MULTIPLE_TRAIT_H


/// 0=uses QTLPOSITION_samelocation 1= uses QTLPOSITION
int QTLLOC;
/// Whether different locations are considered for each trait for each indicator or same locations are considered
int DiffLocation;       
/// Residual variance matrix for multiple triats
double **SIGMA;
/// Main Effect Variance scale for beta ~ N(0,cI) c=scale;
double SCALE;
//*******************************************************************************************


/// @addtogroup qbGeno
/// @{

/// QTL genotypes for the MT code
int ***GENO;

/// @}

/// @addtogroup qbPheno
/// @{

/// phenotypic data for the MT code
double **Y;
/// The number of phenotypes to be jointly analyzed
int NPHENO;            
/// 
int MULTIPLE;

/// @}

/// genotypic valyes for the MT code
double **GVALUE;


//***************************************************************************************
/// @defgroup qbOp_MT Multiple Traits variants of the qb.op variables
/// @{


/// coefficients of QTL main effects
double ****COEF;
/// QTL position indicators, position is GRID[QCHR[L]][QLOC[L]]
int **QLOC;
/// chromosomes that QTL locate
int **QCHR;

    /// @defgroup qbOpCovariate_MT Multiple Traits variant of qb.op.covariate
    /// @{
    /// effects of fixed covariates
    double  **FIX;
    /// effects of random covariates
    double  ***RAN;
    /// variances of random covariates
    double  **VRAN;
    /// interactions of QTL main effects and fixed covariates
    double ****GBYE_FIX;
    /// @}

    /// @defgroup qbOpIterdiag_MT Multiple Traits variant of qb.op.iterdiag
    /// @{
    /// overall mean
    double *AMU;
    /// QTL indicators
    int **GAMMA;
    /// residual variance
    double VE;
    /// @}

    /// @defgroup qbOpMainloci_MT Multiple Traits variant of qb.op.mainloci
    /// @{
    /// main effects
    double  ***MAIN;
    /// main effects indicators
    int **GAMMA_MAIN;
    /// prior variance of main effects
    double  ***VMAIN;
    /// @}

    /// @defgroup qbOpPairloci_MT Multiple Traits variant of qb.op.pairloci
    /// @{
    /// epistatic effects
    double *****EPISTATIC;
    /// epistatic effects indicators
    int ***GAMMA_EPISTASIS;
    /// prior variance of epistatic effects
    double *****VEPISTASIS;
    /// @}

    /// @defgroup qbOpGbyE_MT Multiple Traits variant of qb.op.gbye
    /// @{
    /// g-by-e indicators
    double  ***GAMMA_GBYE;
    /// prior variance of g-by-e effects
    double ****V_GBYE_FIX;
    /// @}

/// @}

#endif // GLOBAL_VARS_MULTIPLE_TRAIT_H
