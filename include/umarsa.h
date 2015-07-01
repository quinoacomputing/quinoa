 
/* umarsa.h for ANSI C */
#ifndef UMARSA_H
#define UMARSA_H
 
#include "gdef.h"
#include "unif01.h"


unif01_Gen * umarsa_CreateMarsa90a (int y1, int y2, int y3, int z0,
                                    unsigned int Y0);



unif01_Gen * umarsa_CreateRANMAR (int y1, int y2, int y3, int z0);


#ifdef USE_LONGLONG
   unif01_Gen * umarsa_CreateMother0 (unsigned long x1, unsigned long x2,
      unsigned long x3, unsigned long x4, unsigned long c);
#endif



unif01_Gen * umarsa_CreateCombo (unsigned int x1, unsigned int x2,
                                 unsigned int y1, unsigned int c);



unif01_Gen * umarsa_CreateECG1 (unsigned int x1, unsigned int x2,
                                unsigned int x3);



unif01_Gen * umarsa_CreateECG2 (unsigned int x1, unsigned int x2,
                                unsigned int x3);



unif01_Gen * umarsa_CreateECG3 (unsigned int x1, unsigned int x2,
                                unsigned int x3);



unif01_Gen * umarsa_CreateECG4 (unsigned int x1, unsigned int x2,
                                unsigned int x3);



unif01_Gen * umarsa_CreateMWC97R (unsigned int x0, unsigned int y0);



unif01_Gen * umarsa_CreateULTRA (unsigned int s1, unsigned int s2,
                                 unsigned int s3, unsigned int s4);



unif01_Gen * umarsa_CreateSupDup73 (unsigned int x0, unsigned int y0);



unif01_Gen * umarsa_CreateSupDup96Add (unsigned int x0, unsigned int y0,
                                       unsigned int c);



unif01_Gen * umarsa_CreateSupDup96Xor (unsigned int x0, unsigned int y0,
                                       unsigned int c);



#ifdef USE_LONGLONG
unif01_Gen * umarsa_CreateSupDup64Add (ulonglong x0, ulonglong y0,
                                       ulonglong a, ulonglong c,
                                       int s1, int s2, int s3);
#endif

    

#ifdef USE_LONGLONG
unif01_Gen * umarsa_CreateSupDup64Xor (ulonglong x0, ulonglong y0,
                                       ulonglong a, ulonglong c,
                                       int s1, int s2, int s3);
#endif



unif01_Gen * umarsa_CreateKISS93 (unsigned int x0, unsigned int y0,
                                  unsigned int z0);

 

unif01_Gen * umarsa_CreateKISS96 (unsigned int x0, unsigned int y0,
                                  unsigned int z1, unsigned int z2);



unif01_Gen * umarsa_CreateKISS99 (unsigned int x0, unsigned int y0,
                                  unsigned int z1, unsigned int z2);



unif01_Gen * umarsa_Create4LFIB99 (unsigned int T[256]);



unif01_Gen * umarsa_Create3SHR99 (unsigned int y0);



unif01_Gen * umarsa_CreateSWB99 (unsigned int T[256], int b);


void umarsa_DeleteGen (unif01_Gen *gen);

 
#endif
 

