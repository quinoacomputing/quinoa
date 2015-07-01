 
/* fvaria.h  for ANSI C */
#ifndef FVARIA_H
#define FVARIA_H
 
#include "ffam.h"
#include "fres.h"
#include "fcho.h"


extern long fvaria_MaxN;
extern long fvaria_Maxn;
extern long fvaria_Maxk;
extern long fvaria_MaxK;


void fvaria_SampleMean1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                         long n, int r,
                         int Nr, int j1, int j2, int jstep);



void fvaria_SampleCorr1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                         long N, int r, int k,
                         int Nr, int j1, int j2, int jstep);



void fvaria_SampleProd1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                         long N, int r, int t,
                         int Nr, int j1, int j2, int jstep);



void fvaria_SumLogs1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                      long N, int r,
                      int Nr, int j1, int j2, int jstep);



void fvaria_SumCollector1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                           long N, int r, double g,
                           int Nr, int j1, int j2, int jstep);



void fvaria_Appearance1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                         long N, int r, int s, int L,
                         int Nr, int j1, int j2, int jstep); 



void fvaria_WeightDistrib1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
                            long N, long n, int r, long k,
                            double alpha, double beta,
                            int Nr, int j1, int j2, int jstep);

  
#endif
 

