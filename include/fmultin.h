 
/* fmultin.h for ANSI C */
#ifndef FMULTIN_H
#define FMULTIN_H
 
#include "gdef.h"
#include "ftab.h"
#include "ffam.h"
#include "fres.h"
#include "fcho.h"
#include "smultin.h"


extern long fmultin_Maxn;


typedef struct {
   smultin_Param *Par;
   fres_Cont *PowDiv[smultin_MAX_DELTA];
   fres_Poisson *Coll;
   fres_Poisson *Empty;
   fres_Poisson *Balls[1 + smultin_MAXB];
   ftab_Table *COApprox;
} fmultin_Res;



fmultin_Res * fmultin_CreateRes (smultin_Param *par);



void fmultin_DeleteRes (fmultin_Res *res);



fcho_Cho * fmultin_CreateEC_DT (long N, int t, double EC);
fcho_Cho * fmultin_CreateEC_2HT (long N, int t, double EC);
fcho_Cho * fmultin_CreateEC_2L (long N, double EC);
fcho_Cho * fmultin_CreateEC_T (long N, double EC);



void fmultin_DeleteEC (fcho_Cho *cho);



fcho_Cho * fmultin_CreateDens_DT (int t, double R);
fcho_Cho * fmultin_CreateDens_2HT (int t, double R);
fcho_Cho * fmultin_CreateDens_2L (double R);
fcho_Cho * fmultin_CreateDens_T (double R);



void fmultin_DeleteDens (fcho_Cho *cho);



fcho_Cho * fmultin_CreatePer_DT (int t, double R);
fcho_Cho * fmultin_CreatePer_2HT (int t, double R);
fcho_Cho * fmultin_CreatePer_2L (double R);
fcho_Cho * fmultin_CreatePer_T (double R);



void fmultin_DeletePer (fcho_Cho *cho);


void fmultin_Serial1 (ffam_Fam *fam, smultin_Param *par,
                      fmultin_Res *res, fcho_Cho2 *cho,
                      long N, int r, int t, lebool Sparse,
                      int Nr, int j1, int j2, int jstep);



void fmultin_SerialOver1 (ffam_Fam *fam, smultin_Param *par,
                          fmultin_Res *res, fcho_Cho2 *cho,
                          long N, int r, int t, lebool Sparse,
                          int Nr, int j1, int j2, int jstep);



void fmultin_SerialBits1 (ffam_Fam *fam, smultin_Param *par,
                          fmultin_Res *res, fcho_Cho2 *cho,
                          long N, int r, int s, lebool Sparse,
                          int Nr, int j1, int j2, int jstep);



void fmultin_SerialBitsOver1 (ffam_Fam *fam, smultin_Param *par,
                              fmultin_Res *res, fcho_Cho2 *cho,
                              long N, int r, int s, lebool Sparse,
                              int Nr, int j1, int j2, int jstep);



void fmultin_Permut1 (ffam_Fam *fam, smultin_Param *par,
                      fmultin_Res *res, fcho_Cho2 *cho,
                      long N, int r, lebool Sparse,
                      int Nr, int j1, int j2, int jstep);

 
#endif
 

