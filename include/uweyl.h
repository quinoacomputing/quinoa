
 
#ifndef UWEYL_H
#define UWEYL_H
/* uweyl.h for ANSI C */
 
#include "unif01.h"


unif01_Gen * uweyl_CreateWeyl (double alpha, long n0);


unif01_Gen * uweyl_CreateNWeyl (double alpha, long n0);



unif01_Gen * uweyl_CreateSNWeyl (long m, double alpha, long n0);


void uweyl_DeleteGen (unif01_Gen *gen);

 
#endif
 

