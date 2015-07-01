 
/* scomp.h  for ANSI C */
#ifndef SCOMP_H
#define SCOMP_H
 
#include "unif01.h"
#include "sres.h"

typedef struct {

   sres_Basic *JumpNum;
   sres_Chi2 *JumpSize;
   sres_Chi2 *LinComp;

} scomp_Res;



scomp_Res * scomp_CreateRes (void);



void scomp_DeleteRes (scomp_Res *res);

void scomp_LinearComp (unif01_Gen *gen, scomp_Res *res,
                       long N, long n, int r, int s);



void scomp_LempelZiv (unif01_Gen *gen, sres_Basic *res,
                      long N, int k, int r, int s);

 
#endif
 

