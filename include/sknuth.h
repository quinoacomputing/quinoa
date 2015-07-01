 
/* sknuth.h  for ANSI C */
#ifndef SKNUTH_H
#define SKNUTH_H
 
#include "gdef.h"
#include "unif01.h"
#include "sres.h"


typedef struct {

   sres_Chi2 *Chi;
   sres_Basic *Bas;

} sknuth_Res1;



sknuth_Res1 * sknuth_CreateRes1 (void);



void sknuth_DeleteRes1 (sknuth_Res1 *res);

typedef struct {

   sres_Poisson *Pois;
   sres_Basic *Bas;

} sknuth_Res2;



sknuth_Res2 * sknuth_CreateRes2 (void);



void sknuth_DeleteRes2 (sknuth_Res2 *res);


void sknuth_Serial (unif01_Gen *gen, sres_Chi2 *res,
                    long N, long n, int r, long d, int t);



void sknuth_SerialSparse (unif01_Gen *gen, sres_Chi2 *res,
                          long N, long n, int r, long d, int t);



void sknuth_Permutation (unif01_Gen *gen, sres_Chi2 *res,
                         long N, long n, int r, int t);



void sknuth_Gap (unif01_Gen *gen, sres_Chi2 *res,
                 long N, long n, int r, double Alpha, double Beta);



void sknuth_SimpPoker (unif01_Gen *gen, sres_Chi2 *res,
                       long N, long n, int r, int d, int k);



void sknuth_CouponCollector (unif01_Gen *gen, sres_Chi2 *res,
                             long N, long n, int r, int d);



void sknuth_Run (unif01_Gen *gen, sres_Chi2 *res,
                 long N, long n, int r, lebool Up);



void sknuth_RunIndep (unif01_Gen *gen, sres_Chi2 *res,
                      long N, long n, int r, lebool Up);



void sknuth_MaxOft (unif01_Gen *gen, sknuth_Res1 *res,
                    long N, long n, int r, int d, int t);



void sknuth_Collision (unif01_Gen *gen, sknuth_Res2 *res,
                       long N, long n, int r, long d, int t);



void sknuth_CollisionPermut (unif01_Gen *gen, sknuth_Res2 *res,
                             long N, long n, int r, int t);

 
#endif
 

