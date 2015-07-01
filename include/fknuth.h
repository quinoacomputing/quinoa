 
/* fknuth.h for ANSI C */
#ifndef FKNUTH_H
#define FKNUTH_H
 
#include "gdef.h"
#include "ffam.h"
#include "fres.h"
#include "fcho.h"


extern long fknuth_Maxn;


typedef struct {
   fres_Cont *Chi;
   fres_Cont *AD;
} fknuth_Res1;



fknuth_Res1 * fknuth_CreateRes1 (void);



void fknuth_DeleteRes1 (fknuth_Res1 *res);


void fknuth_Serial1 (void);



void fknuth_SerialSparse1 (void);



void fknuth_Collision1 (void);



void fknuth_Permutation1 (void);



void fknuth_CollisionPermut1 (void);



void fknuth_Gap1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                  long N, int r, double Alpha, double Beta,
                  int Nr, int j1, int j2, int jstep);



void fknuth_SimpPoker1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                        long N, int r, int d, int k,
                        int Nr, int j1, int j2, int jstep);



void fknuth_CouponCollector1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                              long N, int r, int d,
                              int Nr, int j1, int j2, int jstep);



void fknuth_Run1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                  long N, int r, lebool Up, lebool Indep,
                  int Nr, int j1, int j2, int jstep);




void fknuth_MaxOft1 (ffam_Fam *fam, fknuth_Res1 *res, fcho_Cho *cho,
                     long N, int r, int d, int t,
                     int Nr, int j1, int j2, int jstep);
 
#endif
 
