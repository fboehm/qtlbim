
//********************************************************************

#include <stdio.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>

#include "GlobalVars.h"
#include "RInterfaces.h"

//***************************************************************************
void ROutputManager(char **iterFile,char **covFile,char **mainFile,char **pairFile,char **gbyeFile,char **devFile)
{
   strcpy(covfile,covFile[0]);
   strcpy(mainfile,mainFile[0]);
   strcpy(pairfile,pairFile[0]);
   strcpy(gbyefile,gbyeFile[0]);
   strcpy(iterfile,iterFile[0]);
   strcpy(devfile,devFile[0]);
}
