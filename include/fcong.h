 
/* fcong.h  for ANSI C */
#ifndef FCONG_H
#define FCONG_H
 
#include "ffam.h"


ffam_Fam * fcong_CreateLCG (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateLCGPow2 (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateMRG2 (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateMRG3 (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateCombL2 (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateCombWH2 (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateInvImpl (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateInvImpl2a (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateInvImpl2b (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateInvExpl (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateInvExpl2a (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateInvExpl2b (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateInvMRG2 (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateCubic1 (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateCombCubic2 (char *fname, int i1, int i2, int istep);



ffam_Fam * fcong_CreateCombCubLCG (char *fname, int i1, int i2, int istep);


void fcong_DeleteLCG     (ffam_Fam *fam);
void fcong_DeleteLCGPow2 (ffam_Fam *fam);
void fcong_DeleteMRG2    (ffam_Fam *fam);
void fcong_DeleteMRG3    (ffam_Fam *fam);
void fcong_DeleteCombL2  (ffam_Fam *fam);
void fcong_DeleteCombWH2 (ffam_Fam *fam);
void fcong_DeleteInvImpl   (ffam_Fam *fam);
void fcong_DeleteInvImpl2a (ffam_Fam *fam);
void fcong_DeleteInvImpl2b (ffam_Fam *fam);
void fcong_DeleteInvExpl   (ffam_Fam *fam);
void fcong_DeleteInvExpl2a (ffam_Fam *fam);
void fcong_DeleteInvExpl2b (ffam_Fam *fam);
void fcong_DeleteInvMRG2   (ffam_Fam *fam);
void fcong_DeleteCubic1    (ffam_Fam *fam);
void fcong_DeleteCombCubic2  (ffam_Fam *fam);
void fcong_DeleteCombCubLCG  (ffam_Fam *fam);


 
#endif
 

