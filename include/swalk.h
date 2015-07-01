 
/* swalk.h  for ANSI C */
#ifndef SWALK_H
#define SWALK_H
 
#include "bitset.h"
#include "unif01.h"
#include "sres.h"


typedef struct {

   long L0, L1, L;


   sres_Chi2 **H;
   sres_Chi2 **M;
   sres_Chi2 **J;
   sres_Chi2 **R;
   sres_Chi2 **C;


   long imax;


   char *name;


   void *work;


} swalk_Res;



swalk_Res * swalk_CreateRes (void);



void swalk_DeleteRes (swalk_Res *res);


void swalk_RandomWalk1 (unif01_Gen *gen, swalk_Res *res, long N, long n,
                        int r, int s, long L0, long L1);



void swalk_RandomWalk1a (unif01_Gen *gen, swalk_Res *res, long N, long n,
                         int r, int s, int t, long L, bitset_BitSet C);



void swalk_VarGeoP (unif01_Gen *gen, sres_Chi2 *res,
                    long N, long n, int r, double Mu);



void swalk_VarGeoN (unif01_Gen *gen, sres_Chi2 *res,
                    long N, long n, int r, double Mu);

  
#endif
 

