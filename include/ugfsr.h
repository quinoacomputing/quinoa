
 
/* ugfsr.h for ANSI C */

#ifndef UGFSR_H
#define UGFSR_H
 
#include "unif01.h"


unif01_Gen * ugfsr_CreateGFSR3 (unsigned int k, unsigned int r,
                                unsigned int l, unsigned long S[]);



unif01_Gen * ugfsr_CreateToot73 (unsigned long S[]);



unif01_Gen * ugfsr_CreateKirk81 (long s);



unif01_Gen * ugfsr_CreateRipley90 (long s);



unif01_Gen * ugfsr_CreateFushimi (int k, int r, int s);



unif01_Gen * ugfsr_CreateFushimi90 (int s);



unif01_Gen * ugfsr_CreateGFSR5 (unsigned int k, unsigned int r1,
                                unsigned int r2, unsigned int r3,
                                unsigned int l, unsigned long S[]);



unif01_Gen * ugfsr_CreateZiff98 (unsigned long S[]);


unif01_Gen * ugfsr_CreateTGFSR (unsigned int k, unsigned int r,
                                unsigned int l, unsigned long Av,
                                unsigned long S[]);



unif01_Gen * ugfsr_CreateT800 (unsigned long S[]);



unif01_Gen * ugfsr_CreateTGFSR2 (unsigned int k, unsigned int r,
                                 unsigned int l, unsigned int s,
                                 unsigned int t, unsigned long Av, 
                                 unsigned long Bv, unsigned long Cv, 
                                 unsigned long S[]);



unif01_Gen * ugfsr_CreateTT400 (unsigned long S[]);



unif01_Gen * ugfsr_CreateTT403 (unsigned long S[]);



unif01_Gen * ugfsr_CreateTT775 (unsigned long S[]);



unif01_Gen * ugfsr_CreateTT800 (unsigned long S[]);



unif01_Gen * ugfsr_CreateTT800M94 (unsigned long S[]);



unif01_Gen * ugfsr_CreateTT800M96 (unsigned long S[]);



unif01_Gen * ugfsr_CreateMT19937_98 (unsigned long seed);



unif01_Gen * ugfsr_CreateMT19937_02 (unsigned long seed,
                                     unsigned long Key[], int len);


void ugfsr_DeleteGFSR5 (unif01_Gen * gen);



void ugfsr_DeleteGen (unif01_Gen *gen);

 
#endif
 

