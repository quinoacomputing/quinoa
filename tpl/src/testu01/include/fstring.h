 
/* fstring.h for ANSI C */
#ifndef FSTRING_H
#define FSTRING_H
 
#include "ffam.h"
#include "fres.h"
#include "fcho.h"


extern long fstring_Maxn, fstring_MaxL;


typedef struct {
   fres_Cont *BLen;
   fres_Disc *GLen;
} fstring_Res1;



fstring_Res1 * fstring_CreateRes1 (void);



void fstring_DeleteRes1 (fstring_Res1 *res);


typedef struct {
   fres_Cont *NBits;
   fres_Cont *NRuns;
} fstring_Res2;



fstring_Res2 * fstring_CreateRes2 (void);



void fstring_DeleteRes2 (fstring_Res2 *res);


void fstring_Period1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                      long N, int r, int s,
                      int Nr, int j1, int j2, int jstep);



void fstring_Run1 (ffam_Fam *fam, fstring_Res2 *res, fcho_Cho *cho,
                   long N, int r, int s,
                   int Nr, int j1, int j2, int jstep);



void fstring_AutoCor1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                       long N, int r, int s, int d,
                       int Nr, int j1, int j2, int jstep);



void fstring_LongHead1 (ffam_Fam *fam, fstring_Res1 *res, fcho_Cho2 *cho,
                        long N, long n, int r, int s, long L,
                        int Nr, int j1, int j2, int jstep);



void fstring_HamWeight1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
                         long N, long n, int r, int s, long L,
                         int Nr, int j1, int j2, int jstep);



void fstring_HamWeight2 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
                         long N, long n, int r, int s, long L,
                         int Nr, int j1, int j2, int jstep);



void fstring_HamCorr1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
                       long N, long n, int r, int s, long L,
                       int Nr, int j1, int j2, int jstep);



void fstring_HamIndep1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
                        long N, long n, int r, int s, long L,
                        int Nr, int j1, int j2, int jstep);


 
#endif
 

