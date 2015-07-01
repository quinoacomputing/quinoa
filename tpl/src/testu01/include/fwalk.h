 
/* fwalk.h  for ANSI C */
#ifndef FWALK_H
#define FWALK_H
 
#include "ffam.h"
#include "fres.h"
#include "fcho.h"


extern long fwalk_Maxn;
extern long fwalk_MaxL;
extern double fwalk_MinMu;


typedef struct {
   fres_Cont *H;
   fres_Cont *M;
   fres_Cont *J;
   fres_Cont *R;
   fres_Cont *C;
} fwalk_Res1;



fwalk_Res1 * fwalk_CreateRes1 (void);



void fwalk_DeleteRes1 (fwalk_Res1 *res);


void fwalk_RWalk1 (ffam_Fam *fam, fwalk_Res1 *res, fcho_Cho2 *cho,
                   long N, long n, int r, int s, long L,
                   int Nr, int j1, int j2, int jstep);



void fwalk_VarGeoP1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
                     long N, long n, int r, double Mu,
                     int Nr, int j1, int j2, int jstep);



void fwalk_VarGeoN1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
                     long N, long n, int r, double Mu,
                     int Nr, int j1, int j2, int jstep);

  
#endif
 

