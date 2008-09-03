#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>

void outbredF2(char **genfile,int *F0inds,int *F1inds,int *F2inds) {
	FILE *fp1,*fp2;
	int i,j,k,l,m,I,J;
	int nmarkers,unique;
	char popcode[4][16],missing[16],***gen,***dat1,***dat2,***dat3;
	char str1[16],str2[16],str3[16],str4[16];
	int **genotype;
	int F0,F1,F2;
	int aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii;
	char geno[2][14][16];

	F0=F0inds[0];
	F1=F1inds[0];
	F2=F2inds[0];
	if((fp1=fopen(genfile[0],"r"))==NULL) exit(1);
	fscanf(fp1,"%d",&nmarkers);
	for (i=0;i<nmarkers;i++) fscanf(fp1,"%*s");
	for (i=0;i<4;i++) fscanf(fp1,"%s",popcode[i]);
	fscanf(fp1,"%*s %*s");
	fscanf(fp1,"%s",missing);
	gen=(char ***)S_alloc((F0+F1+F2),sizeof(char **));
	for (i=0;i<(F0+F1+F2);i++) {
		gen[i]=(char **)S_alloc((nmarkers*2+5),sizeof(char *));
		for (j=0;j<(nmarkers*2+5);j++) gen[i][j]=(char *)S_alloc(16,sizeof(char));
	}
	for (i=0;i<(F0+F1+F2);i++) {
		fscanf(fp1,"%s %s %s %s %s",gen[i][0],gen[i][1],gen[i][2],gen[i][3],gen[i][4]);
		for (j=0;j<nmarkers;j++) fscanf(fp1,"%s %s",gen[i][((j*2)+5)],gen[i][((j*2)+6)]);
	}
	if(fclose(fp1)==EOF) exit(1);

	dat1=(char ***)S_alloc(F0,sizeof(char **));
	for (i=0;i<F0;i++) {
		dat1[i]=(char **)S_alloc((2*nmarkers+2),sizeof(char *));
		for (j=0;j<(2*nmarkers+2);j++) dat1[i][j]=(char *)S_alloc(16,sizeof(char));
	}
	dat2=(char ***)S_alloc(F1,sizeof(char **));
	for (i=0;i<F1;i++) {
		dat2[i]=(char **)S_alloc((2*nmarkers+1),sizeof(char *));
		for (j=0;j<(2*nmarkers+1);j++) dat2[i][j]=(char *)S_alloc(16,sizeof(char));
	}
	dat3=(char ***)S_alloc(F2,sizeof(char **));
	for (i=0;i<F2;i++) {
		dat3[i]=(char **)S_alloc((2*nmarkers+1),sizeof(char *));
		for (j=0;j<(2*nmarkers+1);j++) dat3[i][j]=(char *)S_alloc(16,sizeof(char));
	}
	k=0;
	l=0;
	m=0;
	for (i=0;i<(F0+F1+F2);i++) {
		if ((!strcmp(gen[i][4],popcode[0])) || (!strcmp(gen[i][4],popcode[1]))) {
			strcpy(dat1[k][0],gen[i][0]);
			strcpy(dat1[k][1],gen[i][4]);
			for (j=0;j<nmarkers;j++) {
				strcpy(dat1[k][((j*2)+2)],gen[i][((j*2)+5)]);
				strcpy(dat1[k][((j*2)+3)],gen[i][((j*2)+6)]);
			}
			k++;
		}
		if (!strcmp(gen[i][4],popcode[2])) {
			strcpy(dat2[l][0],gen[i][0]);
			for (j=0;j<nmarkers;j++) {
				strcpy(dat2[l][((j*2)+1)],gen[i][((j*2)+5)]);
				strcpy(dat2[l][((j*2)+2)],gen[i][((j*2)+6)]);
			}
			l++;
		}
		if (!strcmp(gen[i][4],popcode[3])) {
			strcpy(dat3[m][0],gen[i][0]);
			for (j=0;j<nmarkers;j++) {
				strcpy(dat3[m][((j*2)+1)],gen[i][((j*2)+5)]);
				strcpy(dat3[m][((j*2)+2)],gen[i][((j*2)+6)]);
			}
			m++;
		}
	}
	genotype=(int **)S_alloc(F2,sizeof(int *));
	for (i=0;i<F2;i++) genotype[i]=(int *)S_alloc(nmarkers,sizeof(int));
	for (i=0;i<F2;i++)
		for (j=0;j<nmarkers;j++)
			genotype[i][j]=0;
	/*for each individual*/
	for (I=0;I<F2;I++) {
		/*row corresponding to the F2 individual*/
		aaa=0;
		while (strcmp(dat3[I][0],gen[aaa][0])) aaa++;
		/*row for the F1 sire*/
		bbb=0;
		while (strcmp(gen[aaa][1],gen[bbb][0])) bbb++;
		/*row for the F1 dam*/
		ccc=0;
		while (strcmp(gen[aaa][2],gen[ccc][0])) ccc++;
		/*row for the F1 sire*/
		ddd=0;
		while (strcmp(gen[aaa][1],dat2[ddd][0])) ddd++;
		/*row for the F1 dam*/
		eee=0;
		while (strcmp(gen[aaa][2],dat2[eee][0])) eee++;
		/*row for the paternal grandfather*/
		fff=0;
		while (strcmp(gen[bbb][1],dat1[fff][0])) fff++;
		/*row for the paternal grandmother*/
		ggg=0;
		while (strcmp(gen[bbb][2],dat1[ggg][0])) ggg++;
		/*row for the maternal grandfather*/
		hhh=0;
		while (strcmp(gen[ccc][1],dat1[hhh][0])) hhh++;
		/*row for the maternal grandmother*/
		iii=0;
		while (strcmp(gen[ccc][2],dat1[iii][0])) iii++;
		/*for each marker*/
		for (J=0;J<nmarkers;J++) {
			/*order is F2 individual, F1 sire, F1 dam, pgf, pgm, mgf, mgm*/
			for (i=0;i<=1;i++)
				for (j=0;j<=13;j++)
					strcpy(geno[i][j],missing);
			/*look up the 14 alleles*/
			strcpy(geno[0][0],dat3[I][((2*J)+1)]);
			strcpy(geno[0][1],dat3[I][((2*J)+2)]);
			strcpy(geno[0][2],dat2[ddd][((2*J)+1)]);
			strcpy(geno[0][3],dat2[ddd][((2*J)+2)]);
			strcpy(geno[0][4],dat2[eee][((2*J)+1)]);
			strcpy(geno[0][5],dat2[eee][((2*J)+2)]);
			strcpy(geno[0][6],dat1[fff][((2*J)+2)]);
			strcpy(geno[0][7],dat1[fff][((2*J)+3)]);
			strcpy(geno[0][8],dat1[ggg][((2*J)+2)]);
			strcpy(geno[0][9],dat1[ggg][((2*J)+3)]);
			strcpy(geno[0][10],dat1[hhh][((2*J)+2)]);
			strcpy(geno[0][11],dat1[hhh][((2*J)+3)]);
			strcpy(geno[0][12],dat1[iii][((2*J)+2)]);
			strcpy(geno[0][13],dat1[iii][((2*J)+3)]);
			/*look up the founders' lines*/
			strcpy(geno[1][6],dat1[fff][1]);
			strcpy(geno[1][7],dat1[fff][1]);
			strcpy(geno[1][8],dat1[ggg][1]);
			strcpy(geno[1][9],dat1[ggg][1]);
			strcpy(geno[1][10],dat1[hhh][1]);
			strcpy(geno[1][11],dat1[hhh][1]);
			strcpy(geno[1][12],dat1[iii][1]);
			strcpy(geno[1][13],dat1[iii][1]);
			/*force the order of founders to be line A, line B, line B, line A*/
			if (!strcmp(dat1[fff][1],popcode[1])) {
				strcpy(str1,geno[0][6]);
				strcpy(str2,geno[0][7]);
				strcpy(str3,geno[1][6]);
				strcpy(str4,geno[1][7]);
				strcpy(geno[0][6],geno[0][8]);
				strcpy(geno[0][7],geno[0][9]);
				strcpy(geno[1][6],geno[1][8]);
				strcpy(geno[1][7],geno[1][9]);
				strcpy(geno[0][8],str1);
				strcpy(geno[0][9],str2);
				strcpy(geno[1][8],str3);
				strcpy(geno[1][9],str4);
			}
			if (!strcmp(dat1[hhh][1],popcode[0])) {
				strcpy(str1,geno[0][10]);
				strcpy(str2,geno[0][11]);
				strcpy(str3,geno[1][10]);
				strcpy(str4,geno[1][11]);
				strcpy(geno[0][10],geno[0][12]);
				strcpy(geno[0][11],geno[0][13]);
				strcpy(geno[1][10],geno[1][12]);
				strcpy(geno[1][11],geno[1][13]);
				strcpy(geno[0][12],str1);
				strcpy(geno[0][13],str2);
				strcpy(geno[1][12],str3);
				strcpy(geno[1][13],str4);
			}

/************************************************************************************************************************************/

			/*determine F0 to F1 sire descent for two alleles and no missingness*/
			unique=1;
			if (strcmp(geno[0][6],geno[0][7])) unique++;
			if (strcmp(geno[0][6],geno[0][8]) && strcmp(geno[0][7],geno[0][8])) unique++;
			if (strcmp(geno[0][6],geno[0][9]) && strcmp(geno[0][7],geno[0][9]) && strcmp(geno[0][8],geno[0][9])) unique++;
			if ((unique==2) && strcmp(geno[0][6],missing) && strcmp(geno[0][7],missing) && strcmp(geno[0][8],missing) && strcmp(geno[0][9],missing) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing)) {
				if (!strcmp(geno[0][2],geno[0][3])) {	/*if the F1 sire is homozygous*/
					strcpy(geno[1][2],geno[1][6]);
					strcpy(geno[1][3],geno[1][8]);
				} else {
					if ((!strcmp(geno[0][6],geno[0][7])) && strcmp(geno[0][8],geno[0][9]) && (!strcmp(geno[0][6],geno[0][8])) && strcmp(geno[0][7],geno[0][9])) {	/*AAxAB*/
						if (!strcmp(geno[0][2],geno[0][6])) {	/*F1 sire is AB*/
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][9]);
						} else {	/*F1 sire is BA*/
							strcpy(geno[1][2],geno[1][9]);
							strcpy(geno[1][3],geno[1][6]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						}
					}
					if ((!strcmp(geno[0][6],geno[0][7])) && strcmp(geno[0][8],geno[0][9]) && strcmp(geno[0][6],geno[0][8]) && (!strcmp(geno[0][7],geno[0][9]))) {	/*AAxBA*/
						if (!strcmp(geno[0][2],geno[0][6])) {	/*F1 sire is AB*/
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][8]);
						} else {	/*F1 sire is BA*/
							strcpy(geno[1][2],geno[1][8]);
							strcpy(geno[1][3],geno[1][6]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						}
					}
					if ((!strcmp(geno[0][6],geno[0][7])) && (!strcmp(geno[0][8],geno[0][9])) && strcmp(geno[0][6],geno[0][8]) && strcmp(geno[0][7],geno[0][9])) {	/*AAxBB*/
						if (!strcmp(geno[0][2],geno[0][6])) {	/*F1 sire is AB*/
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][8]);
						} else {	/*F1 sire is BA*/
							strcpy(geno[1][2],geno[1][8]);
							strcpy(geno[1][3],geno[1][6]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						}
					}
					if (strcmp(geno[0][6],geno[0][7]) && (!strcmp(geno[0][8],geno[0][9])) && (!strcmp(geno[0][6],geno[0][8])) && strcmp(geno[0][7],geno[0][9])) {	/*ABxAA*/
						if (!strcmp(geno[0][2],geno[0][6])) {	/*F1 sire is AB*/
							strcpy(geno[1][2],geno[1][8]);
							strcpy(geno[1][3],geno[1][7]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						} else {	/*F1 sire is BA*/
							strcpy(geno[1][2],geno[1][7]);
							strcpy(geno[1][3],geno[1][8]);
						}	
					}
					if (strcmp(geno[0][6],geno[0][7]) && strcmp(geno[0][8],geno[0][9]) && (!strcmp(geno[0][6],geno[0][8])) && (!strcmp(geno[0][7],geno[0][9]))) {	/*ABxAB or BAxBA*/
						strcpy(geno[1][2],missing);
						strcpy(geno[1][3],missing);
					}
					if (strcmp(geno[0][6],geno[0][7]) && strcmp(geno[0][8],geno[0][9]) && strcmp(geno[0][6],geno[0][8]) && strcmp(geno[0][7],geno[0][9])) {	/*ABxBA or BAxAB*/
						strcpy(geno[1][2],missing);
						strcpy(geno[1][3],missing);
					}
					if (strcmp(geno[0][6],geno[0][7]) && (!strcmp(geno[0][8],geno[0][9])) && strcmp(geno[0][6],geno[0][8]) && (!strcmp(geno[0][7],geno[0][9]))) {	/*ABxBB*/
						if (!strcmp(geno[0][2],geno[0][6])) {	/*F1 sire is AB*/
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][8]);
						} else {	/*F1 sire is BA*/
							strcpy(geno[1][2],geno[1][8]);
							strcpy(geno[1][3],geno[1][6]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						}	
					}
					if (strcmp(geno[0][6],geno[0][7]) && (!strcmp(geno[0][8],geno[0][9])) && strcmp(geno[0][6],geno[0][8]) && (!strcmp(geno[0][7],geno[0][9]))) {	/*BAxAA*/
						if (!strcmp(geno[0][2],geno[0][6])) {	/*F1 sire is BA*/
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][8]);
						} else {	/*F1 sire is AB*/
							strcpy(geno[1][2],geno[1][8]);
							strcpy(geno[1][3],geno[1][6]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						}	
					}
					if (strcmp(geno[0][6],geno[0][7]) && (!strcmp(geno[0][8],geno[0][9])) && (!strcmp(geno[0][6],geno[0][8])) && strcmp(geno[0][7],geno[0][9])) {	/*BAxBB*/
						if (!strcmp(geno[0][2],geno[0][6])) {	/*F1 sire is BA*/
							strcpy(geno[1][2],geno[1][8]);
							strcpy(geno[1][3],geno[1][7]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						} else {	/*F1 sire is AB*/
							strcpy(geno[1][2],geno[1][7]);
							strcpy(geno[1][3],geno[1][8]);
						}	
					}
					if ((!strcmp(geno[0][6],geno[0][7])) && (!strcmp(geno[0][8],geno[0][9])) && strcmp(geno[0][6],geno[0][8]) && strcmp(geno[0][7],geno[0][9])) {	/*BBxAA*/
						if (!strcmp(geno[0][2],geno[0][6])) {	/*F1 sire is BA*/
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][8]);
						} else {	/*F1 sire is AB*/
							strcpy(geno[1][2],geno[1][8]);
							strcpy(geno[1][3],geno[1][6]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						}	
					}
					if ((!strcmp(geno[0][6],geno[0][7])) && strcmp(geno[0][8],geno[0][9]) && strcmp(geno[0][6],geno[0][8]) && (!strcmp(geno[0][7],geno[0][9]))) {	/*BBxAB*/
						if (!strcmp(geno[0][2],geno[0][6])) {	/*F1 sire is BA*/
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][8]);
						} else {	/*F1 sire is AB*/
							strcpy(geno[1][2],geno[1][8]);
							strcpy(geno[1][3],geno[1][6]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						}	
					}
					if ((!strcmp(geno[0][6],geno[0][7])) && strcmp(geno[0][8],geno[0][9]) && (!strcmp(geno[0][6],geno[0][8])) && strcmp(geno[0][7],geno[0][9])) {	/*BBxBA*/
						if (!strcmp(geno[0][2],geno[0][6])) {	/*F1 sire is BA*/
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][9]);
						} else {	/*F1 sire is AB*/
							strcpy(geno[1][2],geno[1][9]);
							strcpy(geno[1][3],geno[1][6]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						}	
					}
				}
			}

/************************************************************************************************************************************/

			/*determine F0 to F1 dam descent for two alleles and no missingness*/
			unique=1;
			if (strcmp(geno[0][10],geno[0][11])) unique++;
			if (strcmp(geno[0][10],geno[0][12]) && strcmp(geno[0][11],geno[0][12])) unique++;
			if (strcmp(geno[0][10],geno[0][13]) && strcmp(geno[0][11],geno[0][13]) && strcmp(geno[0][12],geno[0][13])) unique++;
			if ((unique==2) && strcmp(geno[0][10],missing) && strcmp(geno[0][11],missing) && strcmp(geno[0][12],missing) && strcmp(geno[0][13],missing) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing)) {
				if (!strcmp(geno[0][4],geno[0][5])) {	/*if the F1 dam is homozygous*/
					strcpy(geno[1][4],geno[1][10]);
					strcpy(geno[1][5],geno[1][12]);
				} else {
					if ((!strcmp(geno[0][10],geno[0][11])) && strcmp(geno[0][12],geno[0][13]) && (!strcmp(geno[0][10],geno[0][12])) && strcmp(geno[0][11],geno[0][13])) {	/*AAxAB*/
						if (!strcmp(geno[0][4],geno[0][10])) {	/*F1 dam is AB*/
							strcpy(geno[1][4],geno[1][10]);
							strcpy(geno[1][5],geno[1][13]);
						} else {	/*F1 dam is BA*/
							strcpy(geno[1][4],geno[1][13]);
							strcpy(geno[1][5],geno[1][10]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						}
					}
					if ((!strcmp(geno[0][10],geno[0][11])) && strcmp(geno[0][12],geno[0][13]) && strcmp(geno[0][10],geno[0][12]) && (!strcmp(geno[0][11],geno[0][13]))) {	/*AAxBA*/
						if (!strcmp(geno[0][4],geno[0][10])) {	/*F1 dam is AB*/
							strcpy(geno[1][4],geno[1][10]);
							strcpy(geno[1][5],geno[1][12]);
						} else {	/*F1 dam is BA*/
							strcpy(geno[1][4],geno[1][12]);
							strcpy(geno[1][5],geno[1][10]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						}
					}
					if ((!strcmp(geno[0][10],geno[0][11])) && (!strcmp(geno[0][12],geno[0][13])) && strcmp(geno[0][10],geno[0][12]) && strcmp(geno[0][11],geno[0][13])) {	/*AAxBB*/
						if (!strcmp(geno[0][4],geno[0][10])) {	/*F1 dam is AB*/
							strcpy(geno[1][4],geno[1][10]);
							strcpy(geno[1][5],geno[1][12]);
						} else {	/*F1 dam is BA*/
							strcpy(geno[1][4],geno[1][12]);
							strcpy(geno[1][5],geno[1][10]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str1);
						}
					}
					if (strcmp(geno[0][10],geno[0][11]) && (!strcmp(geno[0][12],geno[0][13])) && (!strcmp(geno[0][10],geno[0][12])) && strcmp(geno[0][11],geno[0][13])) {	/*ABxAA*/
						if (!strcmp(geno[0][4],geno[0][10])) {	/*F1 dam is AB*/
							strcpy(geno[1][4],geno[1][12]);
							strcpy(geno[1][5],geno[1][11]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						} else {	/*F1 dam is BA*/
							strcpy(geno[1][4],geno[1][11]);
							strcpy(geno[1][5],geno[1][12]);
						}	
					}
					if (strcmp(geno[0][10],geno[0][11]) && strcmp(geno[0][12],geno[0][13]) && (!strcmp(geno[0][10],geno[0][12])) && (!strcmp(geno[0][11],geno[0][13]))) {	/*ABxAB or BAxBA*/
						strcpy(geno[1][4],missing);
						strcpy(geno[1][5],missing);
					}
					if (strcmp(geno[0][10],geno[0][11]) && strcmp(geno[0][12],geno[0][13]) && strcmp(geno[0][10],geno[0][12]) && strcmp(geno[0][11],geno[0][13])) {	/*ABxBA or BAxAB*/
						strcpy(geno[1][4],missing);
						strcpy(geno[1][5],missing);
					}
					if (strcmp(geno[0][10],geno[0][11]) && (!strcmp(geno[0][12],geno[0][13])) && strcmp(geno[0][10],geno[0][12]) && (!strcmp(geno[0][11],geno[0][13]))) {	/*ABxBB*/
						if (!strcmp(geno[0][4],geno[0][10])) {	/*F1 dam is AB*/
							strcpy(geno[1][4],geno[1][10]);
							strcpy(geno[1][5],geno[1][12]);
						} else {	/*F1 dam is BA*/
							strcpy(geno[1][4],geno[1][12]);
							strcpy(geno[1][5],geno[1][10]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						}	
					}
					if (strcmp(geno[0][10],geno[0][11]) && (!strcmp(geno[0][12],geno[0][13])) && strcmp(geno[0][10],geno[0][12]) && (!strcmp(geno[0][11],geno[0][13]))) {	/*BAxAA*/
						if (!strcmp(geno[0][4],geno[0][10])) {	/*F1 dam is BA*/
							strcpy(geno[1][4],geno[1][10]);
							strcpy(geno[1][5],geno[1][12]);
						} else {	/*F1 dam is AB*/
							strcpy(geno[1][4],geno[1][12]);
							strcpy(geno[1][5],geno[1][10]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						}	
					}
					if (strcmp(geno[0][10],geno[0][11]) && (!strcmp(geno[0][12],geno[0][13])) && (!strcmp(geno[0][10],geno[0][12])) && strcmp(geno[0][11],geno[0][13])) {	/*BAxBB*/
						if (!strcmp(geno[0][4],geno[0][10])) {	/*F1 dam is BA*/
							strcpy(geno[1][4],geno[1][12]);
							strcpy(geno[1][5],geno[1][11]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						} else {	/*F1 dam is AB*/
							strcpy(geno[1][4],geno[1][11]);
							strcpy(geno[1][5],geno[1][12]);
						}	
					}
					if ((!strcmp(geno[0][10],geno[0][11])) && (!strcmp(geno[0][12],geno[0][13])) && strcmp(geno[0][10],geno[0][12]) && strcmp(geno[0][11],geno[0][13])) {	/*BBxAA*/
						if (!strcmp(geno[0][4],geno[0][10])) {	/*F1 dam is BA*/
							strcpy(geno[1][4],geno[1][10]);
							strcpy(geno[1][5],geno[1][12]);
						} else {	/*F1 dam is AB*/
							strcpy(geno[1][4],geno[1][12]);
							strcpy(geno[1][5],geno[1][10]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						}	
					}
					if ((!strcmp(geno[0][10],geno[0][11])) && strcmp(geno[0][12],geno[0][13]) && strcmp(geno[0][10],geno[0][12]) && (!strcmp(geno[0][11],geno[0][13]))) {	/*BBxAB*/
						if (!strcmp(geno[0][4],geno[0][10])) {	/*F1 dam is BA*/
							strcpy(geno[1][4],geno[1][10]);
							strcpy(geno[1][5],geno[1][12]);
						} else {	/*F1 dam is AB*/
							strcpy(geno[1][4],geno[1][12]);
							strcpy(geno[1][5],geno[1][10]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						}	
					}
					if ((!strcmp(geno[0][10],geno[0][11])) && strcmp(geno[0][12],geno[0][13]) && (!strcmp(geno[0][10],geno[0][12])) && strcmp(geno[0][11],geno[0][13])) {	/*BBxBA*/
						if (!strcmp(geno[0][4],geno[0][10])) {	/*F1 dam is BA*/
							strcpy(geno[1][4],geno[1][10]);
							strcpy(geno[1][5],geno[1][13]);
						} else {	/*F1 dam is AB*/
							strcpy(geno[1][4],geno[1][13]);
							strcpy(geno[1][5],geno[1][10]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						}	
					}
				}
			}

/************************************************************************************************************************************/

			/*determine F0 to F1 sire descent for two alleles and missingness*/
			if ((!strcmp(geno[0][6],missing)) && (!strcmp(geno[0][7],missing)) && strcmp(geno[0][8],missing) && strcmp(geno[0][9],missing) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing)) {
				if (!strcmp(geno[0][8],geno[0][9])) {
					if (!strcmp(geno[0][2],geno[0][3])) {
						strcpy(geno[1][2],geno[1][6]);
						strcpy(geno[1][3],geno[1][8]);
					} else {
						if (!strcmp(geno[0][2],geno[0][8])) {
							strcpy(geno[1][2],geno[1][8]);
							strcpy(geno[1][3],geno[1][6]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						} else {
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][8]);
						}
					}
				} else {
					if (!strcmp(geno[0][2],geno[0][3])) {
						strcpy(geno[1][2],geno[1][6]);
						strcpy(geno[1][3],geno[1][8]);
					} else {
						if (((!strcmp(geno[0][2],geno[0][8])) && strcmp(geno[0][3],geno[0][9])) || (strcmp(geno[0][3],geno[0][8]) && (!strcmp(geno[0][2],geno[0][9])))) {
							strcpy(geno[1][2],geno[1][8]);
							strcpy(geno[1][3],geno[1][6]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						}
						if (((!strcmp(geno[0][3],geno[0][8])) && strcmp(geno[0][2],geno[0][9])) || (strcmp(geno[0][2],geno[0][8]) && (!strcmp(geno[0][3],geno[0][9])))) {
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][8]);
						}
					}
				}
			}
			if (strcmp(geno[0][6],missing) && strcmp(geno[0][7],missing) && (!strcmp(geno[0][8],missing)) && (!strcmp(geno[0][9],missing)) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing)) {
				if (!strcmp(geno[0][6],geno[0][7])) {
					if (!strcmp(geno[0][2],geno[0][3])) {
						strcpy(geno[1][2],geno[1][6]);
						strcpy(geno[1][3],geno[1][8]);
					} else {
						if (!strcmp(geno[0][2],geno[0][6])) {
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][8]);
						} else {
							strcpy(geno[1][2],geno[1][8]);
							strcpy(geno[1][3],geno[1][6]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						}
					}
				} else {
					if (!strcmp(geno[0][2],geno[0][3])) {
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][8]);
					} else {
						if (((!strcmp(geno[0][2],geno[0][6])) && strcmp(geno[0][3],geno[0][7])) || (strcmp(geno[0][3],geno[0][6]) && (!strcmp(geno[0][2],geno[0][7])))) {
							strcpy(geno[1][2],geno[1][6]);
							strcpy(geno[1][3],geno[1][8]);
						}
						if (((!strcmp(geno[0][3],geno[0][6])) && strcmp(geno[0][2],geno[0][7])) || (strcmp(geno[0][2],geno[0][6]) && (!strcmp(geno[0][3],geno[0][7])))) {
							strcpy(geno[1][2],geno[1][8]);
							strcpy(geno[1][3],geno[1][6]);
							strcpy(str1,geno[0][2]);
							strcpy(str2,geno[1][2]);
							strcpy(geno[0][2],geno[0][3]);
							strcpy(geno[1][2],geno[1][3]);
							strcpy(geno[0][3],str1);
							strcpy(geno[1][3],str2);
						}
					}
				}
			}

/************************************************************************************************************************************/

			/*determine F0 to F1 dam descent for two alleles and missingness*/
			if ((!strcmp(geno[0][10],missing)) && (!strcmp(geno[0][11],missing)) && strcmp(geno[0][12],missing) && strcmp(geno[0][13],missing) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing)) {
				if (!strcmp(geno[0][12],geno[0][13])) {
					if (!strcmp(geno[0][4],geno[0][5])) {
						strcpy(geno[1][4],geno[1][10]);
						strcpy(geno[1][5],geno[1][12]);
					} else {
						if (!strcmp(geno[0][4],geno[0][12])) {
							strcpy(geno[1][4],geno[1][12]);
							strcpy(geno[1][5],geno[1][10]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						} else {
							strcpy(geno[1][4],geno[1][10]);
							strcpy(geno[1][5],geno[1][12]);
						}
					}
				} else {
					if (!strcmp(geno[0][4],geno[0][5])) {
						strcpy(geno[1][4],geno[1][10]);
						strcpy(geno[1][5],geno[1][12]);
					} else {
						if (((!strcmp(geno[0][4],geno[0][12])) && strcmp(geno[0][5],geno[0][13])) || (strcmp(geno[0][5],geno[0][12]) && (!strcmp(geno[0][4],geno[0][13])))) {
							strcpy(geno[1][4],geno[1][12]);
							strcpy(geno[1][5],geno[1][10]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						}
						if (((!strcmp(geno[0][5],geno[0][12])) && strcmp(geno[0][4],geno[0][13])) || (strcmp(geno[0][4],geno[0][12]) && (!strcmp(geno[0][5],geno[0][13])))) {
							strcpy(geno[1][4],geno[1][10]);
							strcpy(geno[1][5],geno[1][12]);
						}
					}
				}
			} else if (strcmp(geno[0][10],missing) && strcmp(geno[0][11],missing) && (!strcmp(geno[0][12],missing)) && (!strcmp(geno[0][13],missing)) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing)) {
				if (!strcmp(geno[0][10],geno[0][11])) {
					if (!strcmp(geno[0][4],geno[0][5])) {
						strcpy(geno[1][4],geno[1][10]);
						strcpy(geno[1][5],geno[1][12]);
					} else {
						if (!strcmp(geno[0][4],geno[0][10])) {
							strcpy(geno[1][4],geno[1][10]);
							strcpy(geno[1][5],geno[1][12]);
						} else {
							strcpy(geno[1][4],geno[1][12]);
							strcpy(geno[1][5],geno[1][10]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[1][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						}
					}
				} else {
					if (!strcmp(geno[0][4],geno[0][5])) {
						strcpy(geno[1][4],geno[1][10]);
						strcpy(geno[1][5],geno[1][12]);
					} else {
						if (((!strcmp(geno[0][4],geno[0][10])) && strcmp(geno[0][5],geno[0][11])) || (strcmp(geno[0][5],geno[0][10]) && (!strcmp(geno[0][4],geno[0][11])))) {
							strcpy(geno[1][4],geno[1][10]);
							strcpy(geno[1][5],geno[1][12]);
						}
						if (((!strcmp(geno[0][5],geno[0][10])) && strcmp(geno[0][4],geno[0][11])) || (strcmp(geno[0][4],geno[0][10]) && (!strcmp(geno[0][5],geno[0][11])))) {
							strcpy(geno[1][4],geno[1][12]);
							strcpy(geno[1][5],geno[1][10]);
							strcpy(str1,geno[0][4]);
							strcpy(str2,geno[1][4]);
							strcpy(geno[0][4],geno[0][5]);
							strcpy(geno[1][4],geno[0][5]);
							strcpy(geno[0][5],str1);
							strcpy(geno[1][5],str2);
						}
					}
				}
			}

/************************************************************************************************************************************/

			/*determine F0 to F1 sire descent for four alleles*/
			unique=1;
			if (strcmp(geno[0][6],geno[0][7])) unique++;
			if (strcmp(geno[0][6],geno[0][8]) && strcmp(geno[0][7],geno[0][8])) unique++;
			if (strcmp(geno[0][6],geno[0][9]) && strcmp(geno[0][7],geno[0][9]) && strcmp(geno[0][8],geno[0][9])) unique++;
			if ((unique==4) && strcmp(geno[0][6],missing) && strcmp(geno[0][7],missing) && strcmp(geno[0][8],missing) && strcmp(geno[0][9],missing) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing)) {
				if ((!strcmp(geno[0][2],geno[0][6])) || (!strcmp(geno[0][2],geno[0][7]))) {
					strcpy(geno[1][2],geno[1][6]);
					strcpy(geno[1][3],geno[1][8]);
				} else if ((!strcmp(geno[0][2],geno[0][8])) || (!strcmp(geno[0][2],geno[0][9]))) {
					strcpy(geno[1][2],geno[1][8]);
					strcpy(geno[1][3],geno[1][6]);
					strcpy(str1,geno[0][2]);
					strcpy(str2,geno[1][2]);
					strcpy(geno[0][2],geno[0][3]);
					strcpy(geno[1][2],geno[1][3]);
					strcpy(geno[0][3],str1);
					strcpy(geno[1][3],str2);
				}
			}

/************************************************************************************************************************************/

			/*determine F0 to F1 dam descent for four alleles*/
			unique=1;
			if (strcmp(geno[0][10],geno[0][11])) unique++;
			if (strcmp(geno[0][10],geno[0][12]) && strcmp(geno[0][11],geno[0][12])) unique++;
			if (strcmp(geno[0][10],geno[0][13]) && strcmp(geno[0][11],geno[0][13]) && strcmp(geno[0][12],geno[0][13])) unique++;
			if ((unique==4) && strcmp(geno[0][10],missing) && strcmp(geno[0][11],missing) && strcmp(geno[0][12],missing) && strcmp(geno[0][13],missing) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing)) {
				if ((!strcmp(geno[0][4],geno[0][10])) || (!strcmp(geno[0][4],geno[0][11]))) {
					strcpy(geno[1][4],geno[1][10]);
					strcpy(geno[1][5],geno[1][12]);
				} else if ((!strcmp(geno[0][4],geno[0][12])) || (!strcmp(geno[0][4],geno[0][13]))) {
					strcpy(geno[1][4],geno[1][12]);
					strcpy(geno[1][5],geno[1][10]);
					strcpy(str1,geno[0][4]);
					strcpy(str2,geno[1][4]);
					strcpy(geno[0][4],geno[0][5]);
					strcpy(geno[1][4],geno[1][5]);
					strcpy(geno[0][5],str1);
					strcpy(geno[1][5],str2);
				}
			}

/************************************************************************************************************************************/

			/*determine F0 to F1 sire descent for three alleles*/
			unique=1;
			if (strcmp(geno[0][6],geno[0][7])) unique++;
			if (strcmp(geno[0][6],geno[0][8]) && strcmp(geno[0][7],geno[0][8])) unique++;
			if (strcmp(geno[0][6],geno[0][9]) && strcmp(geno[0][7],geno[0][9]) && strcmp(geno[0][8],geno[0][9])) unique++;
			if ((unique==3) && strcmp(geno[0][6],missing) && strcmp(geno[0][7],missing) && strcmp(geno[0][8],missing) && strcmp(geno[0][9],missing) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing)) {
				if (!strcmp(geno[0][6],geno[0][7])) {
					if (!strcmp(geno[0][2],geno[0][6])) {
						strcpy(geno[1][2],geno[1][6]);
						strcpy(geno[1][3],geno[1][8]);
					} else {
						strcpy(geno[1][2],geno[1][8]);
						strcpy(geno[1][3],geno[1][6]);
						strcpy(str1,geno[0][2]);
						strcpy(str2,geno[1][2]);
						strcpy(geno[0][2],geno[0][3]);
						strcpy(geno[1][2],geno[1][3]);
						strcpy(geno[0][3],str1);
						strcpy(geno[1][3],str2);
					}
				} else if (!strcmp(geno[0][8],geno[0][9])) {
					if (!strcmp(geno[0][2],geno[0][8])) {
						strcpy(geno[1][2],geno[1][8]);
						strcpy(geno[1][3],geno[1][6]);
						strcpy(str1,geno[0][2]);
						strcpy(str2,geno[1][2]);
						strcpy(geno[0][2],geno[0][3]);
						strcpy(geno[1][2],geno[1][3]);
						strcpy(geno[0][3],str1);
						strcpy(geno[1][3],str2);
					} else {
						strcpy(geno[1][2],geno[1][6]);
						strcpy(geno[1][3],geno[1][8]);
					}
				} else if (strcmp(geno[0][6],geno[0][7]) && strcmp(geno[0][8],geno[0][9])) {
					if (!strcmp(geno[0][2],geno[0][3])) {
						strcpy(geno[1][2],geno[1][6]);
						strcpy(geno[1][3],geno[1][8]);
					} else {
						if ((!strcmp(geno[0][6],geno[0][8])) || (!strcmp(geno[0][7],geno[0][8]))) {	/*dup==3*/
							if (!strcmp(geno[0][6],geno[0][8])) {
								if ((!strcmp(geno[0][2],geno[0][6])) && (!strcmp(geno[0][3],geno[0][9]))) {
									strcpy(geno[1][2],geno[1][6]);
									strcpy(geno[1][3],geno[1][8]);
								} else if ((!strcmp(geno[0][2],geno[0][9])) && (!strcmp(geno[0][3],geno[0][6]))) {
									strcpy(geno[1][2],geno[1][9]);
									strcpy(geno[1][3],geno[1][6]);
									strcpy(str1,geno[0][2]);
									strcpy(str2,geno[1][2]);
									strcpy(geno[0][2],geno[0][3]);
									strcpy(geno[1][2],geno[1][3]);
									strcpy(geno[0][3],str1);
									strcpy(geno[1][3],str2);
								} else if ((!strcmp(geno[0][2],geno[0][7])) && (!strcmp(geno[0][3],geno[0][8]))) {
									strcpy(geno[1][2],geno[1][6]);
									strcpy(geno[1][3],geno[1][8]);
								} else if ((!strcmp(geno[0][2],geno[0][8])) && (!strcmp(geno[0][3],geno[0][7]))) {
									strcpy(geno[1][2],geno[1][8]);
									strcpy(geno[1][3],geno[1][6]);
									strcpy(str1,geno[0][2]);
									strcpy(str2,geno[1][2]);
									strcpy(geno[0][2],geno[0][3]);
									strcpy(geno[1][2],geno[1][3]);
									strcpy(geno[0][3],str1);
									strcpy(geno[1][3],str2);
								} else if ((!strcmp(geno[0][2],geno[0][7])) && (!strcmp(geno[0][3],geno[0][9]))) {
									strcpy(geno[1][2],geno[1][6]);
									strcpy(geno[1][3],geno[1][8]);
								} else if ((!strcmp(geno[0][2],geno[0][9])) && (!strcmp(geno[0][3],geno[0][7]))) {
									strcpy(geno[1][2],geno[1][8]);
									strcpy(geno[1][3],geno[1][6]);
									strcpy(str1,geno[0][2]);
									strcpy(str2,geno[1][2]);
									strcpy(geno[0][2],geno[0][3]);
									strcpy(geno[1][2],geno[1][3]);
									strcpy(geno[0][3],str1);
									strcpy(geno[1][3],str2);
								}
							} else {
								if ((!strcmp(geno[0][2],geno[0][7])) && (!strcmp(geno[0][3],geno[0][9]))) {
									strcpy(geno[1][2],geno[1][7]);
									strcpy(geno[1][3],geno[1][8]);
								} else if ((!strcmp(geno[0][2],geno[0][9])) && (!strcmp(geno[0][3],geno[0][7]))) {
									strcpy(geno[1][2],geno[1][9]);
									strcpy(geno[1][3],geno[1][7]);
									strcpy(str1,geno[0][2]);
									strcpy(str2,geno[1][2]);
									strcpy(geno[0][2],geno[0][3]);
									strcpy(geno[1][2],geno[1][3]);
									strcpy(geno[0][3],str1);
									strcpy(geno[1][3],str2);
								} else if ((!strcmp(geno[0][2],geno[0][6])) && (!strcmp(geno[0][3],geno[0][8]))) {
									strcpy(geno[1][2],geno[1][6]);
									strcpy(geno[1][3],geno[1][8]);
								} else if ((!strcmp(geno[0][2],geno[0][8])) && (!strcmp(geno[0][3],geno[0][6]))) {
									strcpy(geno[1][2],geno[1][8]);
									strcpy(geno[1][3],geno[1][6]);
									strcpy(str1,geno[0][2]);
									strcpy(str2,geno[1][2]);
									strcpy(geno[0][2],geno[0][3]);
									strcpy(geno[1][2],geno[1][3]);
									strcpy(geno[0][3],str1);
									strcpy(geno[1][3],str2);
								} else if ((!strcmp(geno[0][2],geno[0][6])) && (!strcmp(geno[0][3],geno[0][9]))) {
									strcpy(geno[1][2],geno[1][6]);
									strcpy(geno[1][3],geno[1][9]);
								} else if ((!strcmp(geno[0][2],geno[0][9])) && (!strcmp(geno[0][3],geno[0][6]))) {
									strcpy(geno[1][2],geno[1][9]);
									strcpy(geno[1][3],geno[1][6]);
									strcpy(str1,geno[0][2]);
									strcpy(str2,geno[1][2]);
									strcpy(geno[0][2],geno[0][3]);
									strcpy(geno[1][2],geno[1][3]);
									strcpy(geno[0][3],str1);
									strcpy(geno[1][3],str2);
								}
							}
						} else { /*dup==4*/
							if ((!strcmp(geno[0][6],geno[0][9]))) {
								if ((!strcmp(geno[0][2],geno[0][6])) && (!strcmp(geno[0][3],geno[0][8]))) {
									strcpy(geno[1][2],geno[1][6]);
									strcpy(geno[1][3],geno[1][8]);
								} else if ((!strcmp(geno[0][2],geno[0][8])) && (!strcmp(geno[0][3],geno[0][6]))) {
									strcpy(geno[1][2],geno[1][8]);
									strcpy(geno[1][3],geno[1][6]);
									strcpy(str1,geno[0][2]);
									strcpy(str2,geno[1][2]);
									strcpy(geno[0][2],geno[0][3]);
									strcpy(geno[1][2],geno[1][3]);
									strcpy(geno[0][3],str1);
									strcpy(geno[1][3],str2);
								} else if ((!strcmp(geno[0][2],geno[0][7])) && (!strcmp(geno[0][3],geno[0][9]))) {
									strcpy(geno[1][2],geno[1][6]);
									strcpy(geno[1][3],geno[1][8]);
								} else if ((!strcmp(geno[0][2],geno[0][9])) && (!strcmp(geno[0][3],geno[0][7]))) {
									strcpy(geno[1][2],geno[1][9]);
									strcpy(geno[1][3],geno[1][7]);
									strcpy(str1,geno[0][2]);
									strcpy(str2,geno[1][2]);
									strcpy(geno[0][2],geno[0][3]);
									strcpy(geno[1][2],geno[1][3]);
									strcpy(geno[0][3],str1);
									strcpy(geno[1][3],str2);
								} else if ((!strcmp(geno[0][2],geno[0][7])) && (!strcmp(geno[0][3],geno[0][8]))) {
									strcpy(geno[1][2],geno[1][7]);
									strcpy(geno[1][3],geno[1][8]);
								} else if ((!strcmp(geno[0][2],geno[0][8])) && (!strcmp(geno[0][3],geno[0][7]))) {
									strcpy(geno[1][2],geno[1][8]);
									strcpy(geno[1][3],geno[1][7]);
									strcpy(str1,geno[0][2]);
									strcpy(str2,geno[1][2]);
									strcpy(geno[0][2],geno[0][3]);
									strcpy(geno[1][2],geno[1][3]);
									strcpy(geno[0][3],str1);
									strcpy(geno[1][3],str2);
								}
							} else {
								if ((!strcmp(geno[0][2],geno[0][7])) && (!strcmp(geno[0][3],geno[0][8]))) {
									strcpy(geno[1][2],geno[1][7]);
									strcpy(geno[1][3],geno[1][8]);
								} else if ((!strcmp(geno[0][2],geno[0][8])) && (!strcmp(geno[0][3],geno[0][7]))) {
									strcpy(geno[1][2],geno[1][8]);
									strcpy(geno[1][3],geno[1][7]);
									strcpy(str1,geno[0][2]);
									strcpy(str2,geno[1][2]);
									strcpy(geno[0][2],geno[0][3]);
									strcpy(geno[1][2],geno[1][3]);
									strcpy(geno[0][3],str1);
									strcpy(geno[1][3],str2);
								} else if ((!strcmp(geno[0][2],geno[0][6])) && (!strcmp(geno[0][3],geno[0][9]))) {
									strcpy(geno[1][2],geno[1][6]);
									strcpy(geno[1][3],geno[1][9]);
								} else if ((!strcmp(geno[0][2],geno[0][9])) && (!strcmp(geno[0][3],geno[0][6]))) {
									strcpy(geno[1][2],geno[1][9]);
									strcpy(geno[1][3],geno[1][6]);
									strcpy(str1,geno[0][2]);
									strcpy(str2,geno[1][2]);
									strcpy(geno[0][2],geno[0][3]);
									strcpy(geno[1][2],geno[1][3]);
									strcpy(geno[0][3],str1);
									strcpy(geno[1][3],str2);
								} else if ((!strcmp(geno[0][2],geno[0][6])) && (!strcmp(geno[0][3],geno[0][8]))) {
									strcpy(geno[1][2],geno[1][6]);
									strcpy(geno[1][3],geno[1][8]);
								} else if ((!strcmp(geno[0][2],geno[0][8])) && (!strcmp(geno[0][3],geno[0][6]))) {
									strcpy(geno[1][2],geno[1][8]);
									strcpy(geno[1][3],geno[1][6]);
									strcpy(str1,geno[0][2]);
									strcpy(str2,geno[1][2]);
									strcpy(geno[0][2],geno[0][3]);
									strcpy(geno[1][2],geno[1][3]);
									strcpy(geno[0][3],str1);
									strcpy(geno[1][3],str2);
								}
							}
						}
					}
				}
			}

/************************************************************************************************************************************/

			/*determine F0 to F1 dam descent for three allleles*/
			unique=1;
			if (strcmp(geno[0][10],geno[0][11])) unique++;
			if (strcmp(geno[0][10],geno[0][12]) && strcmp(geno[0][11],geno[0][12])) unique++;
			if (strcmp(geno[0][10],geno[0][13]) && strcmp(geno[0][11],geno[0][13]) && strcmp(geno[0][12],geno[0][13])) unique++;
			if ((unique==3) && strcmp(geno[0][10],missing) && strcmp(geno[0][11],missing) && strcmp(geno[0][12],missing) && strcmp(geno[0][13],missing) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing)) {
				if (!strcmp(geno[0][10],geno[0][11])) {
					if (!strcmp(geno[0][4],geno[0][10])) {
						strcpy(geno[1][4],geno[1][10]);
						strcpy(geno[1][5],geno[1][12]);
					} else {
						strcpy(geno[1][4],geno[1][12]);
						strcpy(geno[1][5],geno[1][10]);
						strcpy(str1,geno[0][4]);
						strcpy(str2,geno[1][4]);
						strcpy(geno[0][4],geno[0][5]);
						strcpy(geno[1][4],geno[1][5]);
						strcpy(geno[0][5],str1);
						strcpy(geno[1][5],str2);
					}
				} else if (!strcmp(geno[0][12],geno[0][13])) {
					if (!strcmp(geno[0][4],geno[0][12])) {
						strcpy(geno[1][4],geno[1][12]);
						strcpy(geno[1][5],geno[1][10]);
						strcpy(str1,geno[0][4]);
						strcpy(str2,geno[1][4]);
						strcpy(geno[0][4],geno[0][5]);
						strcpy(geno[1][4],geno[1][5]);
						strcpy(geno[0][5],str1);
						strcpy(geno[1][5],str2);
					} else {
						strcpy(geno[1][4],geno[1][10]);
						strcpy(geno[1][5],geno[1][12]);
					}
				} else if (strcmp(geno[0][10],geno[0][11]) && strcmp(geno[0][12],geno[0][13])) {
					if (!strcmp(geno[0][4],geno[0][5])) {
						strcpy(geno[1][4],geno[1][10]);
						strcpy(geno[1][5],geno[1][12]);
					} else {
						if ((!strcmp(geno[0][10],geno[0][12])) || (!strcmp(geno[0][11],geno[0][12]))) {	/*dup==3*/
							if (!strcmp(geno[0][10],geno[0][12])) {
								if ((!strcmp(geno[0][4],geno[0][10])) && (!strcmp(geno[0][5],geno[0][13]))) {
									strcpy(geno[1][4],geno[1][10]);
									strcpy(geno[1][5],geno[1][12]);
								} else if ((!strcmp(geno[0][4],geno[0][13])) && (!strcmp(geno[0][5],geno[0][10]))) {
									strcpy(geno[1][4],geno[1][13]);
									strcpy(geno[1][5],geno[1][10]);
									strcpy(str1,geno[0][4]);
									strcpy(str2,geno[1][4]);
									strcpy(geno[0][4],geno[0][5]);
									strcpy(geno[1][4],geno[1][5]);
									strcpy(geno[0][5],str1);
									strcpy(geno[1][5],str2);
								} else if ((!strcmp(geno[0][4],geno[0][11])) && (!strcmp(geno[0][5],geno[0][12]))) {
									strcpy(geno[1][4],geno[1][10]);
									strcpy(geno[1][5],geno[1][12]);
								} else if ((!strcmp(geno[0][4],geno[0][12])) && (!strcmp(geno[0][5],geno[0][11]))) {
									strcpy(geno[1][4],geno[1][12]);
									strcpy(geno[1][5],geno[1][10]);
									strcpy(str1,geno[0][4]);
									strcpy(str2,geno[1][4]);
									strcpy(geno[0][4],geno[0][5]);
									strcpy(geno[1][4],geno[1][5]);
									strcpy(geno[0][5],str1);
									strcpy(geno[1][5],str2);
								} else if ((!strcmp(geno[0][4],geno[0][11])) && (!strcmp(geno[1][5],geno[1][13]))) {
									strcpy(geno[1][4],geno[1][10]);
									strcpy(geno[1][5],geno[1][12]);
								} else if ((!strcmp(geno[0][4],geno[0][13])) && (!strcmp(geno[0][5],geno[0][11]))) {
									strcpy(geno[1][4],geno[1][12]);
									strcpy(geno[1][5],geno[1][10]);
									strcpy(str1,geno[0][4]);
									strcpy(str2,geno[1][4]);
									strcpy(geno[0][4],geno[0][5]);
									strcpy(geno[1][4],geno[1][5]);
									strcpy(geno[0][5],str1);
									strcpy(geno[1][5],str2);
								}
							} else {
								if ((!strcmp(geno[0][4],geno[0][11])) && (!strcmp(geno[0][5],geno[0][13]))) {
									strcpy(geno[1][4],geno[1][11]);
									strcpy(geno[1][5],geno[1][13]);
								} else if ((!strcmp(geno[0][4],geno[0][13])) && (!strcmp(geno[0][5],geno[0][11]))) {
									strcpy(geno[1][4],geno[1][13]);
									strcpy(geno[1][5],geno[1][11]);
									strcpy(str1,geno[0][4]);
									strcpy(str2,geno[1][4]);
									strcpy(geno[0][4],geno[0][5]);
									strcpy(geno[1][4],geno[1][5]);
									strcpy(geno[0][5],str1);
									strcpy(geno[1][5],str2);
								} else if ((!strcmp(geno[0][4],geno[0][10])) && (!strcmp(geno[0][5],geno[0][12]))) {
									strcpy(geno[1][4],geno[1][10]);
									strcpy(geno[1][5],geno[1][12]);
								} else if ((!strcmp(geno[0][4],geno[0][12])) && (!strcmp(geno[0][5],geno[0][10]))) {
									strcpy(geno[1][4],geno[1][12]);
									strcpy(geno[1][5],geno[1][10]);
									strcpy(str1,geno[0][4]);
									strcpy(str2,geno[1][4]);
									strcpy(geno[0][4],geno[0][5]);
									strcpy(geno[1][4],geno[1][5]);
									strcpy(geno[0][5],str1);
									strcpy(geno[1][5],str2);
								} else if ((!strcmp(geno[0][4],geno[0][10])) && (!strcmp(geno[0][5],geno[0][13]))) {
									strcpy(geno[1][4],geno[1][10]);
									strcpy(geno[1][5],geno[1][13]);
								} else if ((!strcmp(geno[0][4],geno[0][13])) && (!strcmp(geno[0][5],geno[0][10]))) {
									strcpy(geno[1][4],geno[1][13]);
									strcpy(geno[1][5],geno[1][10]);
									strcpy(str1,geno[0][4]);
									strcpy(str2,geno[1][4]);
									strcpy(geno[0][4],geno[0][5]);
									strcpy(geno[1][4],geno[1][5]);
									strcpy(geno[0][5],str1);
									strcpy(geno[1][5],str2);
								}
							}
						} else { /*dup==4*/
							if (!strcmp(geno[0][10],geno[0][13])) {
								if ((!strcmp(geno[0][4],geno[0][10])) && (!strcmp(geno[0][5],geno[0][12]))) {
									strcpy(geno[1][4],geno[1][10]);
									strcpy(geno[1][5],geno[1][12]);
								} else if ((!strcmp(geno[0][4],geno[0][12])) && (!strcmp(geno[0][5],geno[0][10]))) {
									strcpy(geno[1][4],geno[1][12]);
									strcpy(geno[1][5],geno[1][10]);
									strcpy(str1,geno[0][4]);
									strcpy(str2,geno[1][4]);
									strcpy(geno[0][4],geno[0][5]);
									strcpy(geno[1][4],geno[1][5]);
									strcpy(geno[0][5],str1);
									strcpy(geno[1][5],str2);
								} else if ((!strcmp(geno[0][4],geno[0][11])) && (!strcmp(geno[0][5],geno[0][13]))) {
									strcpy(geno[1][4],geno[1][11]);
									strcpy(geno[1][5],geno[1][13]);
								} else if ((!strcmp(geno[0][4],geno[0][13])) && (!strcmp(geno[0][5],geno[0][11]))) {
									strcpy(geno[1][4],geno[1][13]);
									strcpy(geno[1][5],geno[1][11]);
									strcpy(str1,geno[0][4]);
									strcpy(str2,geno[1][4]);
									strcpy(geno[0][4],geno[0][5]);
									strcpy(geno[1][4],geno[1][5]);
									strcpy(geno[0][5],str1);
									strcpy(geno[1][5],str2);
								} else if ((!strcmp(geno[0][4],geno[0][11])) && (!strcmp(geno[0][5],geno[0][12]))) {
									strcpy(geno[1][4],geno[1][11]);
									strcpy(geno[1][5],geno[1][12]);
								} else if ((!strcmp(geno[0][4],geno[0][12])) && (!strcmp(geno[0][5],geno[0][11]))) {
									strcpy(geno[1][4],geno[1][12]);
									strcpy(geno[1][5],geno[1][11]);
									strcpy(str1,geno[0][4]);
									strcpy(str2,geno[1][4]);
									strcpy(geno[0][4],geno[0][5]);
									strcpy(geno[1][4],geno[1][5]);
									strcpy(geno[0][5],str1);
									strcpy(geno[1][5],str2);
								}
							} else {
								if ((!strcmp(geno[0][4],geno[0][11])) && (!strcmp(geno[0][5],geno[0][12]))) {
									strcpy(geno[1][4],geno[1][11]);
									strcpy(geno[1][5],geno[1][12]);
								} else if ((!strcmp(geno[0][4],geno[0][12])) && (!strcmp(geno[0][5],geno[0][11]))) {
									strcpy(geno[1][4],geno[1][12]);
									strcpy(geno[1][5],geno[1][11]);
									strcpy(str1,geno[0][4]);
									strcpy(str2,geno[1][4]);
									strcpy(geno[0][4],geno[0][5]);
									strcpy(geno[1][4],geno[1][5]);
									strcpy(geno[0][5],str1);
									strcpy(geno[1][5],str2);
								} else if ((!strcmp(geno[0][4],geno[0][10])) && (!strcmp(geno[0][5],geno[0][13]))) {
									strcpy(geno[1][4],geno[1][10]);
									strcpy(geno[1][5],geno[1][13]);
								} else if ((!strcmp(geno[0][4],geno[0][13])) && (!strcmp(geno[0][5],geno[0][10]))) {
									strcpy(geno[1][4],geno[1][13]);
									strcpy(geno[1][5],geno[1][10]);
									strcpy(str1,geno[0][4]);
									strcpy(str2,geno[1][4]);
									strcpy(geno[0][4],geno[0][5]);
									strcpy(geno[1][4],geno[1][5]);
									strcpy(geno[0][5],str1);
									strcpy(geno[1][5],str2);
								} else if ((!strcmp(geno[0][4],geno[0][10])) && (!strcmp(geno[0][5],geno[0][12]))) {
									strcpy(geno[1][4],geno[1][10]);
									strcpy(geno[1][5],geno[1][12]);
								} else if ((!strcmp(geno[0][4],geno[0][12])) && (!strcmp(geno[0][5],geno[0][10]))) {
									strcpy(geno[1][4],geno[1][12]);
									strcpy(geno[1][5],geno[1][10]);
									strcpy(str1,geno[0][4]);
									strcpy(str2,geno[1][4]);
									strcpy(geno[0][4],geno[0][5]);
									strcpy(geno[1][4],geno[1][5]);
									strcpy(geno[0][5],str1);
									strcpy(geno[1][5],str2);
								}
							}
						}
					}
				}
			}

/************************************************************************************************************************************/

			/*determine F1 to F2 descent for two alleles*/
			unique=1;
			if (strcmp(geno[0][2],geno[0][3])) unique++;
			if (strcmp(geno[0][2],geno[0][4]) && strcmp(geno[0][3],geno[0][4])) unique++;
			if (strcmp(geno[0][2],geno[0][5]) && strcmp(geno[0][3],geno[0][5]) && strcmp(geno[0][4],geno[0][5])) unique++;
			if ((unique==2) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing) && strcmp(geno[0][0],missing) && strcmp(geno[0][1],missing) && strcmp(geno[1][2],missing) && strcmp(geno[1][3],missing) && strcmp(geno[1][4],missing) && strcmp(geno[1][5],missing)) {
				if (strcmp(geno[0][2],geno[0][3]) && strcmp(geno[0][4],geno[0][5])) {
					if (!strcmp(geno[0][2],geno[0][4])) {
						if (!strcmp(geno[0][0],geno[0][1])) genotype[I][J]=2;
						/*if (strcmp(geno[0][0],geno[0][1])) genotype[I][J]=0;*/
					} else if (strcmp(geno[0][2],geno[0][4])) {
						if ((!strcmp(geno[0][0],geno[0][1])) && (!strcmp(geno[0][0],geno[0][2]))) genotype[I][J]=1;
						else if ((!strcmp(geno[0][0],geno[0][1])) && (!strcmp(geno[0][0],geno[0][3]))) genotype[I][J]=3;
						else if (strcmp(geno[0][0],geno[0][1])) genotype[I][J]=2;
					}
				} else if ((!strcmp(geno[0][2],geno[0][3])) && strcmp(geno[0][4],geno[0][5])) {
					if (!strcmp(geno[0][2],geno[0][4])) {
						if (!strcmp(geno[0][0],geno[0][1])) genotype[I][J]=5;
						else if (strcmp(geno[0][0],geno[0][1])) genotype[I][J]=4;
					} else if (strcmp(geno[0][2],geno[0][4])) {
						if (!strcmp(geno[0][0],geno[0][1])) genotype[I][J]=4;
						else if (strcmp(geno[0][0],geno[0][1])) genotype[I][J]=5;
					}
				} else if (strcmp(geno[0][2],geno[0][3]) && (!strcmp(geno[0][4],geno[0][5]))) {
					if ((!strcmp(geno[0][2],geno[0][4]))) {
						if (!strcmp(geno[0][0],geno[0][1])) genotype[I][J]=4;
						else if (strcmp(geno[0][0],geno[0][1])) genotype[I][J]=5;
					} else if (strcmp(geno[0][2],geno[0][4])) {
						if (!strcmp(geno[0][0],geno[0][1])) genotype[I][J]=5;
						else if (strcmp(geno[0][0],geno[0][1])) genotype[I][J]=4;
					}
				}
				/*if ((!strcmp(geno[0][2],geno[0][3])) && (!strcmp(geno[0][4],geno[0][5]))) {
					if (!strcmp(geno[0][2],geno[0][4])) genotype[I][J]=0;
					if (strcmp(geno[0][2],geno[0][4])) genotype[I][J]=0;
				}*/
			} else if ((unique==2) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing) && strcmp(geno[0][0],missing) && strcmp(geno[0][1],missing) && (!strcmp(geno[1][2],missing)) && (!strcmp(geno[1][3],missing)) && strcmp(geno[1][4],missing) && strcmp(geno[1][5],missing)) {
				if (strcmp(geno[0][2],geno[0][3]) && strcmp(geno[0][4],geno[0][5])) {
					if (!strcmp(geno[0][0],geno[0][1])) {
						if (!strcmp(geno[0][0],geno[0][4])) genotype[I][J]=5;
						else if (!strcmp(geno[0][0],geno[0][5])) genotype[I][J]=4;
					}
					/*if (strcmp(geno[0][0],geno[0][1])) genotype[I][J]=0;*/
				} else if ((!strcmp(geno[0][2],geno[0][3])) && (!strcmp(geno[0][2],geno[0][4]))) {
					if (!strcmp(geno[0][0],geno[0][1])) genotype[I][J]=5;
					else if (strcmp(geno[0][0],geno[0][1])) genotype[I][J]=4;
				} else if ((!strcmp(geno[0][2],geno[0][3])) && (!strcmp(geno[0][2],geno[0][5]))) {
					if (!strcmp(geno[0][0],geno[0][1])) genotype[I][J]=4;
					else if (strcmp(geno[0][0],geno[0][1])) genotype[I][J]=5;
				}
				/*if (!strcmp(geno[0][4],geno[0][5])) genotype[I][J]=0;*/
			} else if ((unique==2) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing) && strcmp(geno[0][0],missing) && strcmp(geno[0][1],missing) && strcmp(geno[1][2],missing) && strcmp(geno[1][3],missing) && (!strcmp(geno[1][4],missing)) && (!strcmp(geno[1][5],missing))) {
				if (strcmp(geno[0][2],geno[0][3]) && strcmp(geno[0][4],geno[0][5])) {
					if (!strcmp(geno[0][0],geno[0][1])) {
						if (!strcmp(geno[0][0],geno[0][2])) genotype[I][J]=4;
						else if (!strcmp(geno[0][0],geno[0][3])) genotype[I][J]=5;
					}
					/*if (strcmp(geno[0][0],geno[0][1])) genotype[I][J]=0;*/
				} else if ((!strcmp(geno[0][4],geno[0][5])) && (!strcmp(geno[0][2],geno[0][4]))) {
					if (!strcmp(geno[0][0],geno[0][1])) genotype[I][J]=4;
					else if (strcmp(geno[0][0],geno[0][1])) genotype[I][J]=5;
				} else if ((!strcmp(geno[0][4],geno[0][5])) && (!strcmp(geno[0][3],geno[0][4]))) {
					if (!strcmp(geno[0][0],geno[0][1])) genotype[I][J]=5;
					else if (strcmp(geno[0][0],geno[0][1])) genotype[I][J]=4;
				}
				/*if (!strcmp(geno[0][2],geno[0][3])) genotype[I][J]=0;*/
			}

/************************************************************************************************************************************/

			/*determine F1 to F2 descent for two alleles*/
			if (strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing) && (!strcmp(geno[0][4],missing)) && (!strcmp(geno[0][5],missing)) && strcmp(geno[0][0],missing) && strcmp(geno[0][1],missing) && strcmp(geno[1][2],missing) && strcmp(geno[1][3],missing)) {
				/*if (!strcmp(geno[0][2],geno[0][3])) genotype[I][J]=0;*/
				if (strcmp(geno[0][2],geno[0][3])) {
					if ((!strcmp(geno[0][0],geno[0][2])) && strcmp(geno[0][1],geno[0][3])) genotype[I][J]=4;
					/*if ((!strcmp(geno[0][0],geno[0][2])) && (!strcmp(geno[0][1],geno[0][3]))) genotype[I][J]=0;*/
					else if ((!strcmp(geno[0][1],geno[0][2])) && strcmp(geno[0][0],geno[0][3])) genotype[I][J]=4;
					/*if ((!strcmp(geno[0][1],geno[0][2])) && (!strcmp(geno[0][0],geno[0][3]))) genotype[I][J]=0;*/
					else if ((!strcmp(geno[0][0],geno[0][3])) && strcmp(geno[0][1],geno[0][2])) genotype[I][J]=5;
					/*if ((!strcmp(geno[0][0],geno[0][3])) && (!strcmp(geno[0][1],geno[0][2]))) genotype[I][J]=0;*/
					else if ((!strcmp(geno[0][1],geno[0][3])) && strcmp(geno[0][0],geno[0][2])) genotype[I][J]=5;
					/*if ((!strcmp(geno[0][1],geno[0][3])) && (!strcmp(geno[0][0],geno[0][2]))) genotype[I][J]=0;*/
				}
			} else if ((!strcmp(geno[0][2],missing)) && (!strcmp(geno[0][3],missing)) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing) && strcmp(geno[0][0],missing) && strcmp(geno[0][1],missing) && strcmp(geno[1][4],missing) && strcmp(geno[1][5],missing)) {
				/*if (!strcmp(geno[0][4],geno[0][5])) genotype[I][J]=0;*/
				if (strcmp(geno[0][4],geno[0][5])) {
					if ((!strcmp(geno[0][0],geno[0][4])) && strcmp(geno[0][1],geno[0][5])) genotype[I][J]=5;
					/*if ((!strcmp(geno[0][0],geno[0][4])) && (!strcmp(geno[0][1],geno[0][5]))) genotype[I][J]=0;*/
					else if ((!strcmp(geno[0][1],geno[0][4])) && strcmp(geno[0][0],geno[0][5])) genotype[I][J]=5;
					/*if ((!strcmp(geno[0][1],geno[0][4])) && (!strcmp(geno[0][0],geno[0][5]))) genotype[I][J]=0;*/
					else if ((!strcmp(geno[0][0],geno[0][5])) && strcmp(geno[0][1],geno[0][4])) genotype[I][J]=4;
					/*if ((!strcmp(geno[0][0],geno[0][5])) && (!strcmp(geno[0][1],geno[0][4]))) genotype[I][J]=0;*/
					else if ((!strcmp(geno[0][1],geno[0][5])) && strcmp(geno[0][0],geno[0][4])) genotype[I][J]=4;
					/*if ((!strcmp(geno[0][1],geno[0][5])) && (!strcmp(geno[0][0],geno[0][4]))) genotype[I][J]=0;*/
				}
			}

/************************************************************************************************************************************/

			/*determine F1 to F2 descent for four alleles*/
			unique=1;
			if (strcmp(geno[0][2],geno[0][3])) unique++;
			if (strcmp(geno[0][2],geno[0][4]) && strcmp(geno[0][3],geno[0][4])) unique++;
			if (strcmp(geno[0][2],geno[0][5]) && strcmp(geno[0][3],geno[0][5]) && strcmp(geno[0][4],geno[0][5])) unique++;
			if ((unique==4) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing) && strcmp(geno[0][0],missing) && strcmp(geno[0][1],missing) && strcmp(geno[1][2],missing) && strcmp(geno[1][3],missing) && strcmp(geno[1][4],missing) && strcmp(geno[1][5],missing)) {
				if (!strcmp(geno[0][0],geno[0][2])) {
					if (!strcmp(geno[0][1],geno[0][4])) genotype[I][J]=2;
					else if (!strcmp(geno[0][1],geno[0][5])) genotype[I][J]=1;
				} else if (!strcmp(geno[0][0],geno[0][3])) {
					if (!strcmp(geno[0][1],geno[0][4])) genotype[I][J]=3;
					else if (!strcmp(geno[0][1],geno[0][5])) genotype[I][J]=2;
				} else if (!strcmp(geno[0][0],geno[0][4])) {
					if (!strcmp(geno[0][1],geno[0][2])) genotype[I][J]=2;
					else if (!strcmp(geno[0][1],geno[0][3])) genotype[I][J]=3;
				} else if (!strcmp(geno[0][0],geno[0][5])) {
					if (!strcmp(geno[0][1],geno[0][2])) genotype[I][J]=1;
					else if (!strcmp(geno[0][1],geno[0][3])) genotype[I][J]=2;
				}
			} else if ((unique==4) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing) && strcmp(geno[0][0],missing) && strcmp(geno[0][1],missing) && (((!strcmp(geno[1][2],missing)) && (!strcmp(geno[1][3],missing))) || ((!strcmp(geno[1][4],missing)) && (!strcmp(geno[1][5],missing))))) {
				if ((!strcmp(geno[1][2],missing)) && (!strcmp(geno[1][3],missing))) {
					if ((!strcmp(geno[0][0],geno[0][4])) || (!strcmp(geno[0][1],geno[0][4]))) genotype[I][J]=5;
					else if ((!strcmp(geno[0][0],geno[0][5])) || (!strcmp(geno[0][1],geno[0][5]))) genotype[I][J]=4;
				} else if ((!strcmp(geno[1][4],missing)) && (!strcmp(geno[1][5],missing))) {
					if ((!strcmp(geno[0][0],geno[0][2])) || (!strcmp(geno[0][1],geno[0][2]))) genotype[I][J]=4;
					else if ((!strcmp(geno[0][0],geno[0][3])) || (!strcmp(geno[0][1],geno[0][3]))) genotype[I][J]=5;
				}
			}

/************************************************************************************************************************************/

			/*determine F1 to F2 descent for three allleles*/
			unique=1;
			if (strcmp(geno[0][2],geno[0][3])) unique++;
			if (strcmp(geno[0][2],geno[0][4]) && strcmp(geno[0][3],geno[0][4])) unique++;
			if (strcmp(geno[0][2],geno[0][5]) && strcmp(geno[0][3],geno[0][5]) && strcmp(geno[0][4],geno[0][5])) unique++;
			if ((unique==3) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing) && strcmp(geno[0][0],missing) && strcmp(geno[0][1],missing) && strcmp(geno[1][2],missing) && strcmp(geno[1][3],missing) && strcmp(geno[1][4],missing) && strcmp(geno[1][5],missing)) {
				if (!strcmp(geno[0][2],geno[0][3])) {
					if (!strcmp(geno[0][0],geno[0][2])) {
						if (!strcmp(geno[0][1],geno[0][4])) genotype[I][J]=5;
						else if (!strcmp(geno[0][1],geno[0][5])) genotype[I][J]=4;
					} else {
						if (!strcmp(geno[0][0],geno[0][4])) genotype[I][J]=5;
						else if (!strcmp(geno[0][0],geno[0][5])) genotype[I][J]=4;
					}
				} else if (!strcmp(geno[0][4],geno[0][5])) {
					if (!strcmp(geno[0][0],geno[0][4])) {
						if (!strcmp(geno[0][1],geno[0][2])) genotype[I][J]=4;
						else if (!strcmp(geno[0][1],geno[0][3])) genotype[I][J]=5;
					} else {
						if (!strcmp(geno[0][0],geno[0][2])) genotype[I][J]=4;
						else if (!strcmp(geno[0][0],geno[0][3])) genotype[I][J]=5;
					}
				} else if (strcmp(geno[0][2],geno[0][3]) && strcmp(geno[0][4],geno[0][5])) {
					if (!strcmp(geno[0][2],geno[0][4])) {
						if ((!strcmp(geno[0][0],geno[0][2])) && (!strcmp(geno[0][1],geno[0][4]))) genotype[I][J]=2;
						else if ((!strcmp(geno[0][0],geno[0][2])) && (!strcmp(geno[0][1],geno[0][5]))) genotype[I][J]=1;
						else if ((!strcmp(geno[0][0],geno[0][5])) && (!strcmp(geno[0][1],geno[0][2]))) genotype[I][J]=1;
						else if ((!strcmp(geno[0][0],geno[0][3])) && (!strcmp(geno[0][1],geno[0][4]))) genotype[I][J]=3;
						else if ((!strcmp(geno[0][0],geno[0][4])) && (!strcmp(geno[0][1],geno[0][3]))) genotype[I][J]=3;
						else if ((!strcmp(geno[0][0],geno[0][3])) && (!strcmp(geno[0][1],geno[0][5]))) genotype[I][J]=2;
						else if ((!strcmp(geno[0][0],geno[0][5])) && (!strcmp(geno[0][1],geno[0][3]))) genotype[I][J]=2;
					} else if (!strcmp(geno[0][3],geno[0][4])) {
						if ((!strcmp(geno[0][0],geno[0][3])) && (!strcmp(geno[0][1],geno[0][4]))) genotype[I][J]=3;
						else if ((!strcmp(geno[0][0],geno[0][3])) && (!strcmp(geno[0][1],geno[0][5]))) genotype[I][J]=2;
						else if ((!strcmp(geno[0][0],geno[0][5])) && (!strcmp(geno[0][1],geno[0][3]))) genotype[I][J]=2;
						else if ((!strcmp(geno[0][0],geno[0][2])) && (!strcmp(geno[0][1],geno[0][4]))) genotype[I][J]=2;
						else if ((!strcmp(geno[0][0],geno[0][4])) && (!strcmp(geno[0][1],geno[0][2]))) genotype[I][J]=2;
						else if ((!strcmp(geno[0][0],geno[0][2])) && (!strcmp(geno[0][1],geno[0][5]))) genotype[I][J]=1;
						else if ((!strcmp(geno[0][0],geno[0][5])) && (!strcmp(geno[0][1],geno[0][2]))) genotype[I][J]=1;
					} else if (!strcmp(geno[0][2],geno[0][5])) {
						if ((!strcmp(geno[0][0],geno[0][2])) && (!strcmp(geno[0][1],geno[0][5]))) genotype[I][J]=1;
						else if ((!strcmp(geno[0][0],geno[0][2])) && (!strcmp(geno[0][1],geno[0][4]))) genotype[I][J]=2;
						else if ((!strcmp(geno[0][0],geno[0][4])) && (!strcmp(geno[0][1],geno[0][2]))) genotype[I][J]=2;
						else if ((!strcmp(geno[0][0],geno[0][3])) && (!strcmp(geno[0][1],geno[0][5]))) genotype[I][J]=2;
						else if ((!strcmp(geno[0][0],geno[0][5])) && (!strcmp(geno[0][1],geno[0][3]))) genotype[I][J]=2;
						else if ((!strcmp(geno[0][0],geno[0][3])) && (!strcmp(geno[0][1],geno[0][4]))) genotype[I][J]=3;
						else if ((!strcmp(geno[0][0],geno[0][4])) && (!strcmp(geno[0][1],geno[0][3]))) genotype[I][J]=3;
					} else if (!strcmp(geno[0][3],geno[0][5])) {
						if ((!strcmp(geno[0][0],geno[0][3])) && (!strcmp(geno[0][1],geno[0][5]))) genotype[I][J]=2;
						else if ((!strcmp(geno[0][0],geno[0][3])) && (!strcmp(geno[0][1],geno[0][4]))) genotype[I][J]=3;
						else if ((!strcmp(geno[0][0],geno[0][4])) && (!strcmp(geno[0][1],geno[0][3]))) genotype[I][J]=3;
						else if ((!strcmp(geno[0][0],geno[0][2])) && (!strcmp(geno[0][1],geno[0][5]))) genotype[I][J]=1;
						else if ((!strcmp(geno[0][0],geno[0][5])) && (!strcmp(geno[0][1],geno[0][2]))) genotype[I][J]=1;
						else if ((!strcmp(geno[0][0],geno[0][2])) && (!strcmp(geno[0][1],geno[0][4]))) genotype[I][J]=2;
						else if ((!strcmp(geno[0][0],geno[0][4])) && (!strcmp(geno[0][1],geno[0][2]))) genotype[I][J]=2;
					}
				}
			} else if ((unique==3) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing) && strcmp(geno[0][0],missing) && strcmp(geno[0][1],missing) && (!strcmp(geno[1][2],missing)) && (!strcmp(geno[1][3],missing))) {
				/*if (!strcmp(geno[0][4],geno[0][5])) genotype[I][J]=0;*/
				if (strcmp(geno[0][4],geno[0][5])) {
					if ((!strcmp(geno[0][0],geno[0][4])) || (!strcmp(geno[0][1],geno[0][4]))) genotype[I][J]=5;
					else if ((!strcmp(geno[0][0],geno[0][5])) || (!strcmp(geno[0][1],geno[0][5]))) genotype[I][J]=4;
				}
			} else if ((unique==3) && strcmp(geno[0][2],missing) && strcmp(geno[0][3],missing) && strcmp(geno[0][4],missing) && strcmp(geno[0][5],missing) && strcmp(geno[0][0],missing) && strcmp(geno[0][1],missing) && (!strcmp(geno[1][4],missing)) && (!strcmp(geno[1][5],missing))) {
				/*if (!strcmp(geno[0][2],geno[0][3])) genotype[I][J]=0;*/
				if (strcmp(geno[0][2],geno[0][3])) {
					if ((!strcmp(geno[0][0],geno[0][2])) || (!strcmp(geno[0][1],geno[0][2]))) genotype[I][J]=4;
					else if ((!strcmp(geno[0][0],geno[0][3])) || (!strcmp(geno[0][1],geno[0][3]))) genotype[I][J]=5;
				}
			}

/************************************************************************************************************************************/
		}	/*end J*/
	}	/*end I*/

	if((fp2=fopen("newgenfile.txt","w"))==NULL) exit(1);
	for (I=0;I<F2;I++) {
		for (J=0;J<nmarkers;J++) fprintf(fp2,"%d ",genotype[I][J]);
		fprintf(fp2,"\n");
	}
	if(fclose(fp2)==EOF) exit(1);
}
