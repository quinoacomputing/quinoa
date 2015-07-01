
 
/*  umrg.h  for ANSI C  */
#ifndef UMRG_H
#define UMRG_H
 
#include "gdef.h"
#include "unif01.h"

unif01_Gen * umrg_CreateMRG (long m, int k, long A[], long S[]);



unif01_Gen * umrg_CreateMRGFloat (long m, int k, long A[], long S[]);



#ifdef USE_GMP
   unif01_Gen * umrg_CreateBigMRG (char *m, int k, char *A[], char *S[]);

#endif


unif01_Gen * umrg_CreateLagFibFloat (int k, int r, char Op, int Lux,
                                     unsigned long S[]);



unif01_Gen * umrg_CreateLagFib (int t, int k, int r, char Op, int Lux,
                                unsigned long S[]);


unif01_Gen * umrg_CreateC2MRG (long m1, long m2, int k, long A1[],
                               long A2[], long S1[], long S2[]);



#ifdef USE_GMP
   unif01_Gen * umrg_CreateBigC2MRG (char *m1, char *m2, int k, char *A1[],
                                     char *A2[], char *S1[], char *S2[]);

#endif


void umrg_DeleteMRG    (unif01_Gen * gen);
void umrg_DeleteMRGFloat (unif01_Gen * gen);
void umrg_DeleteLagFib (unif01_Gen * gen);
void umrg_DeleteLagFibFloat (unif01_Gen * gen);
void umrg_DeleteC2MRG  (unif01_Gen * gen);

#ifdef USE_GMP
   void umrg_DeleteBigMRG (unif01_Gen * gen);
   void umrg_DeleteBigC2MRG (unif01_Gen * gen);
#endif

 
#endif
 

