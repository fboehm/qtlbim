#include "GlobalVars.h"
#include "GlobalVars_SingleTrait.h"
#include "RInterfaces.h"

#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>

#include <stdio.h>

//***************************************************************************
void ROutputManager(char **iterFile,char **covFile,char **mainFile,char **pairFile,char **gbyeFile,char **devFile,char **sigmaFile)
{
   strcpy(covfile,covFile[0]);
   strcpy(mainfile,mainFile[0]);
   strcpy(pairfile,pairFile[0]);
   strcpy(gbyefile,gbyeFile[0]);
   strcpy(iterfile,iterFile[0]);
   strcpy(devfile,devFile[0]);
   strcpy(sigmafile,sigmaFile[0]);
   
}
