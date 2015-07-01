 
/* sstring.h for ANSI C */
#ifndef SSTRING_H
#define SSTRING_H
 
#include "tables.h" 
#include "unif01.h"
#include "sres.h"



#define sstring_MAXD 8



extern lebool sstring_Counters;



extern lebool sstring_CorrFlag;


typedef struct {

   sres_Chi2 *Chi;
   sres_Disc *Disc;


} sstring_Res2;



sstring_Res2 * sstring_CreateRes2 (void);



void sstring_DeleteRes2 (sstring_Res2 *res);



typedef struct {

   sres_Basic *NBits;
   sres_Chi2 *NRuns;
   long *Count0;
   long *Count1;


} sstring_Res3;



sstring_Res3 * sstring_CreateRes3 (void);



void sstring_DeleteRes3 (sstring_Res3 *res);



typedef struct {

   int L;


   tables_StyleType Style;


   long **Counters;


   double **ZCounters;


   int d;


   long XD [sstring_MAXD + 1][2];


   sres_Basic *Block [sstring_MAXD + 1];


   sres_Basic *Bas;  


} sstring_Res;



sstring_Res * sstring_CreateRes (void);



void sstring_DeleteRes (sstring_Res *res);



void sstring_PeriodsInStrings (unif01_Gen *gen, sres_Chi2 *res, 
                               long N, long n, int r, int s);



void sstring_LongestHeadRun (unif01_Gen *gen, sstring_Res2 *res,
                             long N, long n, int r, int s, long L);



void sstring_HammingWeight (unif01_Gen *gen, sres_Chi2 *res,
                            long N, long n, int r, int s, long L);



void sstring_HammingWeight2 (unif01_Gen *gen, sres_Basic *res,
                             long N, long n, int r, int s, long L);



void sstring_HammingCorr (unif01_Gen *gen, sstring_Res *res,
                          long N, long n, int r, int s, int L);



void sstring_HammingIndep (unif01_Gen *gen, sstring_Res *res,
                           long N, long n, int r, int s, int L, int d);



void sstring_Run (unif01_Gen *gen, sstring_Res3 *res,
                  long N, long n, int r, int s);



void sstring_AutoCor (unif01_Gen *gen, sres_Basic *res,
                      long N, long n, int r, int s, int d);


 
#endif
 

