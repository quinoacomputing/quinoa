 
/* sentrop.h  for ANSI C */
#ifndef SENTROP_H
#define SENTROP_H
 
#include "statcoll.h"
#include "gofw.h"
#include "unif01.h"
#include "sres.h"


typedef struct {

   long *Count;


   long jmin;
   long jmax;


   sres_Basic *Bas;


} sentrop_Res;



sentrop_Res * sentrop_CreateRes (void);



void sentrop_DeleteRes (sentrop_Res *res);


void sentrop_EntropyDisc (unif01_Gen *gen, sentrop_Res *res,
                          long N, long n, int r, int s, int L);



void sentrop_EntropyDiscOver (unif01_Gen *gen, sentrop_Res *res,
                              long N, long n, int r, int s, int L);



void sentrop_EntropyDiscOver2 (unif01_Gen *gen, sentrop_Res *res,
                               long N, long n, int r, int s, int L);



void sentrop_EntropyDM (unif01_Gen *gen, sres_Basic *res,
                        long N, long n, int r, long m);



void sentrop_EntropyDMCirc (unif01_Gen *gen, sres_Basic *res,
                            long N, long n, int r, long m);

 
#endif
 

