 
/* ffsr.h  for ANSI C */
#ifndef FFSR_H
#define FFSR_H
 
#include "ffam.h"


ffam_Fam * ffsr_CreateLFSR1 (char *fname, int i1, int i2, int istep);



ffam_Fam * ffsr_CreateLFSR2 (char *fname, int i1, int i2, int istep);



ffam_Fam * ffsr_CreateLFSR3 (char *fname, int i1, int i2, int istep);



ffam_Fam * ffsr_CreateGFSR3 (char *fname, int i1, int i2, int istep);



ffam_Fam * ffsr_CreateGFSR5 (char *fname, int i1, int i2, int istep);



ffam_Fam * ffsr_CreateTGFSR1 (char *fname, int i1, int i2, int istep);



ffam_Fam * ffsr_CreateTausLCG2 (char *fname, int i1, int i2, int istep);


void ffsr_DeleteLFSR1 (ffam_Fam *fam);
void ffsr_DeleteLFSR2 (ffam_Fam *fam);
void ffsr_DeleteLFSR3 (ffam_Fam *fam);
void ffsr_DeleteTausLCG2 (ffam_Fam *fam);


 
#endif
 

