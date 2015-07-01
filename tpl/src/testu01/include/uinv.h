
 
/* uinv.h for ANSI C */
#ifndef UINV_H
#define UINV_H
 
#include "unif01.h"


unif01_Gen * uinv_CreateInvImpl (long m, long a1, long a2, long z0);



unif01_Gen * uinv_CreateInvImpl2a (int e, unsigned long a1,
                                   unsigned long a2, unsigned long z0);



unif01_Gen * uinv_CreateInvImpl2b (int e, unsigned long a1,
                                   unsigned long a2, unsigned long z0);



unif01_Gen * uinv_CreateInvExpl (long m, long a, long c);



unif01_Gen * uinv_CreateInvExpl2a (int e, long a, long c);



unif01_Gen * uinv_CreateInvExpl2b (int e, long a, long c);



unif01_Gen * uinv_CreateInvMRG (long m, int k, long A[], long S[]);



unif01_Gen * uinv_CreateInvMRGFloat (long m, int k, long A[], long S[]);


void uinv_DeleteInvMRG (unif01_Gen * gen);



void uinv_DeleteInvMRGFloat (unif01_Gen * gen);



void uinv_DeleteGen (unif01_Gen * gen);

 
#endif
 

