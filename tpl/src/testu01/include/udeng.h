 
/* udeng.h for ANSI C */

#ifndef UDENG_H
#define UDENG_H
 
#include "unif01.h"



unif01_Gen * udeng_CreateDL00a (unsigned long m, unsigned long b, int k,
                                unsigned long S[]);



unif01_Gen * udeng_CreateDX02a (unsigned long m, unsigned long b, int k,
                                unsigned long S[]);



void udeng_DeleteGen (unif01_Gen * gen);

 
#endif
 

