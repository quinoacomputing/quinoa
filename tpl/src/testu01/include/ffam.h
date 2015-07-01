 
/* ffam.h  for ANSI C */
#ifndef FFAM_H
#define FFAM_H
 
#include "unif01.h"



typedef struct {
   unif01_Gen **Gen;
   int *LSize;
   int *Resol;
   int Ng;
   char *name;
} ffam_Fam;



ffam_Fam * ffam_CreateFam (int Ng, char *name);



void ffam_DeleteFam (ffam_Fam *fam);



void ffam_PrintFam (ffam_Fam *fam);



void ffam_ReallocFam (ffam_Fam *fam, int Ng);



ffam_Fam * ffam_CreateSingle (unif01_Gen *gen, int resol, int i1, int i2);



void ffam_DeleteSingle (ffam_Fam *fam);

 
#include <stdio.h>

FILE * ffam_OpenFile (char *filename, char *defaultfile);


#endif
 

