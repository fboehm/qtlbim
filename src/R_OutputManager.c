

#include "bmq_mcmc.h"
#include "stdio.h"
#include "bmq_genoprob.h"

void R_OutputManager(char **iterFile,char **covFile,char **mainFile,char **pairFile,char **gbyeFile){
    
   strcpy(covfile,covFile[0]);
   strcpy(mainfile,mainFile[0]);
   strcpy(pairfile,pairFile[0]);
   strcpy(gbyefile,gbyeFile[0]);
   strcpy(iterfile,iterFile[0]);
}
