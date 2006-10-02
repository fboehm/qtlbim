
//********************************************************************

#include <stdio.h>
#include <time.h>
#include <math.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>

#include "bmq_mcmc.h"
#include "bmq_genoprob.h"

//***************************************************************************
void R_OutputManager(char **iterFile,char **covFile,char **mainFile,char **pairFile,char **gbyeFile){
    
   strcpy(covfile,covFile[0]);
   strcpy(mainfile,mainFile[0]);
   strcpy(pairfile,pairFile[0]);
   strcpy(gbyefile,gbyeFile[0]);
   strcpy(iterfile,iterFile[0]);
}
