
 
/* ucubic.h for ANSI C */
#ifndef UCUBIC_H
#define UCUBIC_H
 
#include "unif01.h"


unif01_Gen * ucubic_CreateCubic (long m, long a, long b, long c, long d,
                                 long s);



unif01_Gen * ucubic_CreateCubicFloat (long m, long a, long b, long c,
                                      long d, long s);



unif01_Gen * ucubic_CreateCubic1 (long m, long a, long s);



unif01_Gen * ucubic_CreateCubic1Float (long m, long a, long s);



unif01_Gen * ucubic_CreateCombCubic2 (long m1, long m2, long a1, long a2, 
                                      long s1, long s2);



unif01_Gen * ucubic_CreateCubicOut (long m, long a, long c, long s);



unif01_Gen * ucubic_CreateCubicOutFloat (long m, long a, long c, long s);


void ucubic_DeleteGen (unif01_Gen *gen);

 
#endif
 

