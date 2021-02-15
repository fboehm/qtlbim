
<!-- README.md is generated from README.Rmd. Please edit that file -->

# qtlbim

<!-- badges: start -->
<!-- badges: end -->

The goal of `qtlbim` is to perform multivariate bayesian QTL mapping
after the methods developed in several journal articles by Yi, et
al. Currently, it works only with two-parent crosses, but we’d like to
develop methods for multiparent crosses like the Diversity Outbred mice.

This particular repository was forked from the CRAN Github repository
and subsequently updated to work with new versions of gcc.

## Installation

`qtlbim` isn’t on CRAN right now.

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fboehm/qtlbim")
```

Detailed use instructions are available in the vignettes. We reproduce
below part of the text from one vignette.

## Example

This is a basic example:

``` r
library(qtlbim)
#> Loading required package: qtl
#> Loading required package: lattice
#> Loading required package: coda
#> Loading required package: tools
#> Loading required package: MASS
## basic example code
```

This document focuses on the `hyper` dataset from the `qtl` (Broman et
al. 2003) R package, which was initially studied in Sugiyama et
al. (2001). The `hyper` dataset is stored in `qtl` as a `cross` object.
The `R/qtlbim` package processes this `cross` object to create a `qb`
object called `qbHyper`, containing the MCMC samples. The `hyper` demo
shows how this is done.

It is possible to directly load the already saved `qb` object with the
`data` command. Following this by a call to `qb.cross` extracts a
version of the `cross` object used to create the `qb` object.

``` r
data(qbHyper)
hyper <- qb.cross(qbHyper)
```

Alternatively, a object can be created by the following sequence of
commands. First load the hyper data set from `R/qtl`, and subset on the
autosomes, as `R/qtlbim` does not yet handle the X chromosome properly.

``` r
data(hyper)
hyper <- subset(hyper, chr=1:19)
```

To run the MCMC sampler on the `hyper` data we use the command

``` r
hyper <- qb.genoprob(hyper, step=2) 
# Now run the MCMC model selection algorithm.
# This can take several minutes.
qbHyper <- qb.mcmc(hyper, pheno.col = 1, seed = 1616)
#> qb.mcmc is running 60,600 iterations. The current iterations are saved: 
#> 200
#> 400
#> 600
#> 800
#> 1000
#> 1200
#> 1400
#> 1600
#> 1800
#> 2000
#> 2200
#> 2400
#> 2600
#> 2800
#> 3000
#> MCMC sample has been saved to: ./bp_Feb-14-190121.
#> Bayesian MCMC took 0.74 minutes.
#> Warning: Removing external directory ./bp_Feb-14-190121
```

The option `seed` sets the pseudorandom number seed so that this run can
be repeated exactly. The `qb` object called `qbHyper` is used throughout
this vignette.

## Creating Bayesian interval mapping MCMC samples

This section describes in more detail how to create Markov chain Monte
Carlo (MCMC) samples from the Bayesian posterior to be used for QTL
mapping. The next step to mapping with the `R/qtl` package would be to
use the function `calc.genoprob` to create genotype probabilities based
on a Hidden Markov model. However, for Bayesian model selection, we
replace `calc.genoprob` with the `R/qtlbim` function `qb.genoprob`. The
function `qb.genoprob` performs some bookkeeping before calling
`calc.genoprob` with the variable stepwidth option for pseudomarker
positions. The probabilities for genotypes at pseudomarkers and at
markers with missing data are calculated by `calc.genoprob` from the
observed marker data using the multipoint method (Jiang and Zeng, 1997).

The MCMC samples are created by `qb.mcmc` after running `qb.genoprob`.
In the simplest case, MCMC samples are created with the following two
calls:

``` r
hyper <- qb.genoprob(hyper, step=2) 
qbHyper <- qb.mcmc(hyper, pheno.col = 1)
#> qb.mcmc is running 60,600 iterations. The current iterations are saved: 
#> 200
#> 400
#> 600
#> 800
#> 1000
#> 1200
#> 1400
#> 1600
#> 1800
#> 2000
#> 2200
#> 2400
#> 2600
#> 2800
#> 3000
#> MCMC sample has been saved to: ./bp_Feb-14-190206.
#> Bayesian MCMC took 0.75 minutes.
#> Warning: Removing external directory ./bp_Feb-14-190206
```

By default the `qb.mcmc` function prints out progress messages of the
number of iterations completed. These progress messages can be
suppressed by setting `verbose=FALSE`. Arguments for the routines
`qb.data` and `qb.model`, described below, can be passed through
`qb.mcmc`. Otherwise, default values are used. The detail below for
`qb.data.`qb.model`, and`qb.mcmc\` routines could be skipped in favor of
default settings.

The function `qb.data` specifies the traits to be analyzed, their
underlying distribution, the random and/or fixed covariates and whether
to standardize or to use a boxcox transformation. Note that, the cross
object can have several phenotypes and some of which could be used as
covariates.

``` r
qbData <- qb.data(hyper, pheno.col = 1, trait = "normal", fixcov = 0, rancov = 0)
```

The `R/qtlbim` routines handle normal, binary and ordinal data. In
addition, the user can specify fixed (`fixcov`) and random (`rancov`)
covariate(s). \[The `pheno.col`, `fixcov`, and `rancov` values can be
numeric indices to the phenotype names, or character strings with exact
phenotype names.\] Fixed covariates can be included as interacting
covariates with the `intcov` option to `qb.model` (see below).

The function `qb.model` defines the model parameters, using defaults
that work well in most settings. Users are probably most interested in
specifying if `epistasis` is considered, the prior expected number of
main effect QTLs (`main.nqtl`), and the prior expected total number of
QTLs (`mean.nqtl`), which includes additional QTLs with only epistatic
effects. A user may set `main.nqtl` and `mean.nqtl` based on previous
QTL analysis, for example using `R/qtl`. Setting the maximum number of
QTLs overall (`max.nqtl`) or per chromosome (`chr.nqtl`), and setting
the minimum `interval` between linked QTL, can be used to restrict
sampling as needed.

Typically a real data set has several traits which can be considered as
covariates. The `intcov` option specifies which covariate(s) can
interact with QTLs, or equivalently, which environmental factors may
have GxE interactions. The `intcov` should be a vector of 0s and 1s of
the same length as the `fixcov` option specified for `qb.data` (see
above).

``` r
qbModel <- qb.model(hyper, epistasis = TRUE, main.nqtl = 3, 
interval = rep(5,nchr(hyper)), chr.nqtl = rep(2,nchr(hyper)), 
depen = FALSE, prop = c(0.5, 0.1, 0.05))
```

The function `qb.mcmc` creates MCMC samples on the data and model
specified. The results are initially saved in a unique directory under
`mydir`, which is removed at completion of the command. Options for
`qb.data` and `qb.model` can be passed directly to `qb.mcmc`, or as the
objects created above.

``` r
qb <- qb.mcmc(hyper, data = qbData, model = qbModel, mydir = ".", 
              n.iter = 3000, n.thin = 20, genoupdate = TRUE)
#> qb.mcmc is running 60,600 iterations. The current iterations are saved: 
#> 200
#> 400
#> 600
#> 800
#> 1000
#> 1200
#> 1400
#> 1600
#> 1800
#> 2000
#> 2200
#> 2400
#> 2600
#> 2800
#> 3000
#> MCMC sample has been saved to: ./bp_Feb-14-190251.
#> Bayesian MCMC took 0.75 minutes.
#> Warning: Removing external directory ./bp_Feb-14-190251
```

The `genoupdate` option simulates pseudomarker and missing marker
genotypes if `TRUE`, or uses a Haley-Knott (1992) type approach if
`FALSE`; the latter is faster, but not generally recommended if there
are many missing genotypes or selective genotyping. `n.iter` samples are
saved, thinning to one in `n.thin` from the MCMC samples to reduce
serial correlation. That is, `n.iter * n.thin` samples are drawn, after
an initial `n.burnin` samples (1% of total by default) are discarded to
allow the chain to converge closer to the posterior distribution.
