
 
/* uvaria.h for ANSI C */

#ifndef UVARIA_H
#define UVARIA_H
 
#include "unif01.h"


unif01_Gen * uvaria_CreateACORN (int k, double S[]);



unif01_Gen * uvaria_CreateCSD (long v, long s);



unif01_Gen * uvaria_CreateRanrotB (unsigned int seed);



unif01_Gen * uvaria_CreateRey97 (double a1, double a2, double b2, long n0);



unif01_Gen * uvaria_CreateTindo (long b, long Delta, int s, int k);


void uvaria_DeleteACORN (unif01_Gen *gen);



void uvaria_DeleteRanrotB (unif01_Gen *gen);



void uvaria_DeleteGen (unif01_Gen *gen);

 
#endif
 

