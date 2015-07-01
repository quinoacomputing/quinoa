 
/* sspectral.h  for ANSI C */
#ifndef SSPEC_H
#define SSPEC_H
 
#include "statcoll.h"
#include "gofw.h"
#include "unif01.h"
#include "sres.h"

typedef struct {

   sres_Basic *Bas;


   double *Coef;


   long jmin, jmax;


} sspectral_Res;



sspectral_Res * sspectral_CreateRes (void);



void sspectral_DeleteRes (sspectral_Res *res);


void sspectral_Fourier1 (unif01_Gen *gen, sspectral_Res *res,
                         long N, int k, int r, int s);



void sspectral_Fourier2 (unif01_Gen *gen, sspectral_Res *res,
                         long N, int k, int r, int s);



void sspectral_Fourier3 (unif01_Gen *gen, sspectral_Res *res,
                         long N, int k, int r, int s);


 
#endif
 

