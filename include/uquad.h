
 
/* uquad.h for ANSI C */
#ifndef UQUAD_H
#define UQUAD_H
 
#include "unif01.h"


unif01_Gen * uquad_CreateQuadratic (long m, long a, long b, long c, long s);



unif01_Gen * uquad_CreateQuadratic2 (int e, unsigned long a,
    unsigned long b, unsigned long c, unsigned long s);


void uquad_DeleteGen (unif01_Gen *gen);

 
#endif
 

