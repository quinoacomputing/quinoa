 
/* fnpair.h  for ANSI C */
#ifndef FNPAIR_H
#define FNPAIR_H
 
#include "gdef.h"
#include "ffam.h"
#include "fres.h"
#include "fcho.h"
#include "snpair.h"


extern long fnpair_Maxn;


typedef struct {
   ftab_Table *PVal[snpair_StatType_N];
} fnpair_Res1;



fnpair_Res1 * fnpair_CreateRes1 (void);



void fnpair_DeleteRes1 (fnpair_Res1 *res);



fcho_Cho *fnpair_CreateM1 (int maxm);



void fnpair_DeleteM1 (fcho_Cho * cho);



void fnpair_ClosePairs1 (ffam_Fam *fam, fnpair_Res1 *res, fcho_Cho2 *cho,
                         long N, int r, int t, int p, int m,
                         int Nr, int j1, int j2, int jstep);



void fnpair_Bickel1 (ffam_Fam *fam, fnpair_Res1 *res, fcho_Cho *cho,
                     long N, int r, int t, int p, lebool Torus,
                     int Nr, int j1, int j2, int jstep);



void fnpair_BitMatch1 (ffam_Fam *fam, fnpair_Res1 *res, fcho_Cho *cho,
                       long N, int r, int t,
                       int Nr, int j1, int j2, int jstep);

 
#endif
 

