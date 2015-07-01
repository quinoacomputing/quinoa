
 
/* utaus.h for ANSI C */

#ifndef UTAUS_H
#define UTAUS_H
 
#include "gdef.h"
#include "unif01.h"


unif01_Gen * utaus_CreateTaus (unsigned int k, unsigned int q,
                               unsigned int s, unsigned int Y);



unif01_Gen * utaus_CreateTausJ (unsigned int k, unsigned int q,
                                unsigned int s, unsigned int j,
                                unsigned int Y);



#ifdef USE_LONGLONG
   unif01_Gen * utaus_CreateLongTaus (unsigned int k, unsigned int q,
                                      unsigned int s, ulonglong Y1);

#endif


unif01_Gen * utaus_CreateCombTaus2 (
   unsigned int k1, unsigned int k2, unsigned int q1, unsigned int q2,
   unsigned int s1, unsigned int s2, unsigned int Y1, unsigned int Y2);



unif01_Gen * utaus_CreateCombTaus3 (
    unsigned int k1, unsigned int k2, unsigned int k3,
    unsigned int q1, unsigned int q2, unsigned int q3,
    unsigned int s1, unsigned int s2, unsigned int s3,
    unsigned int Y1, unsigned int Y2, unsigned int Y3);



unif01_Gen * utaus_CreateCombTaus3T (
    unsigned int k1, unsigned int k2, unsigned int k3,
    unsigned int q1, unsigned int q2, unsigned int q3,
    unsigned int s1, unsigned int s2, unsigned int s3,
    unsigned int Y1, unsigned int Y2, unsigned int Y3);


void utaus_DeleteGen (unif01_Gen *gen);

 
#endif
 

