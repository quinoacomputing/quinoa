
 
/*  ucarry.h  for ANSI C  */

#ifndef UCARRY_H
#define UCARRY_H
 
#include "gdef.h"
#include "unif01.h"


unif01_Gen * ucarry_CreateAWC (unsigned int r, unsigned int s,
                               unsigned long c, unsigned long m,
                               unsigned long S[]);



unif01_Gen * ucarry_CreateSWB (unsigned int r, unsigned int s,
                               unsigned long c, unsigned long m,
                               unsigned long S[]);



unif01_Gen * ucarry_CreateRanlux (unsigned int L, long s);



#ifdef USE_LONGLONG
   unif01_Gen * ucarry_CreateMWC (unsigned int r, unsigned long c,
                                  unsigned int w, unsigned long A[],
                                  unsigned long S[]);
#endif



unif01_Gen * ucarry_CreateMWCFloat (unsigned int r, unsigned long c,
                                    unsigned int w, unsigned long A[],
                                    unsigned long S[]);



#ifdef USE_LONGLONG
   unif01_Gen * ucarry_CreateMWCfix8r4 (unsigned long c, unsigned long S[]);
#endif



#ifdef USE_LONGLONG
   unif01_Gen * ucarry_CreateMWCfix8r8 (unsigned long c, unsigned long S[]);
#endif



unif01_Gen * ucarry_CreateMWCfixCouture (unsigned int c,
                                         unsigned int S[]);



unif01_Gen * ucarry_CreateSWC (unsigned int r, unsigned int h,
                               unsigned int c, unsigned int w,
                               unsigned int A[], unsigned int S[]);



unif01_Gen * ucarry_CreateMWC1616 (unsigned int a, unsigned int b,
                                   unsigned int x, unsigned int y);



void ucarry_DeleteAWC (unif01_Gen *gen);
void ucarry_DeleteSWB (unif01_Gen *gen);
void ucarry_DeleteRanlux (unif01_Gen *gen);
void ucarry_DeleteMWC (unif01_Gen *gen);
void ucarry_DeleteMWCFloat (unif01_Gen *gen);
void ucarry_DeleteMWCfixCouture (unif01_Gen *gen);
void ucarry_DeleteSWC (unif01_Gen *gen);



void ucarry_DeleteGen (unif01_Gen *gen);


 
#endif
 

