 
/* wdist.h for ANSI C */
#ifndef WDIST_H
#define WDIST_H
 
#include "fmass.h"


typedef double (*wdist_CFUNC) (double [], double);



typedef double (*wdist_DFUNC) (fmass_INFO, long);


double wdist_Normal (double Par[], double x);



double wdist_ChiSquare (double Par[], double x);



double wdist_Unif (double Par[], double x);

 
#endif
 

