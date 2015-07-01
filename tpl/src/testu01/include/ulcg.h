
 
/*  ulcg.h  for ANSI C  */

#ifndef ULCG_H
#define ULCG_H
 
#include "gdef.h"
#include "unif01.h"


unif01_Gen * ulcg_CreateLCG (long m, long a, long c, long s);



unif01_Gen * ulcg_CreateLCGFloat (long m, long a, long c, long s);



#ifdef USE_GMP
   unif01_Gen * ulcg_CreateBigLCG (char *m, char *a, char *c, char *s);

#endif


unif01_Gen * ulcg_CreateLCGWu2 (long m, char o1, unsigned int q, char o2, 
                                unsigned int r, long s);



unif01_Gen * ulcg_CreateLCGPayne (long a, long c, long s);



unif01_Gen * ulcg_CreateLCG2e31m1HD (long a, long s);



unif01_Gen * ulcg_CreateLCG2e31 (long a, long c, long s);



unif01_Gen * ulcg_CreateLCG2e32 (unsigned long a, unsigned long c,
                                 unsigned long s);



unif01_Gen * ulcg_CreatePow2LCG (int e, long a, long c, long s);



#ifdef USE_LONGLONG
unif01_Gen * ulcg_CreateLCG2e48L (ulonglong a, ulonglong c, ulonglong s);



unif01_Gen * ulcg_CreatePow2LCGL (int e, ulonglong a, ulonglong c,
                                  ulonglong s);

#endif



#ifdef USE_GMP
unif01_Gen * ulcg_CreateBigPow2LCG (long e, char *a, char *c, char *s);

#endif


unif01_Gen * ulcg_CreateCombLEC2 (long m1, long m2, long a1, long a2,
                                  long c1, long c2, long s1, long s2);



unif01_Gen * ulcg_CreateCombLEC2Float (long m1, long m2, long a1, long a2,
                                       long c1, long c2, long s1, long s2);



unif01_Gen * ulcg_CreateCombLEC3 (long m1, long m2, long m3, long a1,
                                  long a2, long a3, long c1, long c2,
                                  long c3, long s1, long s2, long s3);



unif01_Gen * ulcg_CreateCombWH2 (long m1, long m2, long a1, long a2,
                                 long c1, long c2, long s1, long s2);



unif01_Gen * ulcg_CreateCombWH2Float (long m1, long m2, long a1, long a2,
                                      long c1, long c2, long s1, long s2);



unif01_Gen * ulcg_CreateCombWH3 (long m1, long m2, long m3, long a1,
                                 long a2, long a3, long c1, long c2,
                                 long c3, long s1, long s2, long s3);



#ifdef USE_GMP
   void ulcg_DeleteBigLCG (unif01_Gen *gen);



   void ulcg_DeleteBigPow2LCG (unif01_Gen *gen);

#endif


void ulcg_DeleteGen (unif01_Gen *gen);

 
#endif
 

