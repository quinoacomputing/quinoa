 
#ifndef UGRANGER_H
#define UGRANGER_H
/* ugranger.h for ANSI C */
 
#include "gdef.h"
#include "unif01.h"


unif01_Gen * ugranger_CreateCombLCGInvExpl (
   long m1, long a1, long c1, long s1, long m2, long a2, long c2);



#ifdef USE_GMP
   unif01_Gen * ugranger_CreateCombBigLCGInvExpl (
      char *m1, char *a1, char *c1, char *s1, long m2, long a2, long c2);

#endif


unif01_Gen * ugranger_CreateCombLCGCub (
   long m1, long a1, long c1, long s1, long m2, long a2, long s2);



#ifdef USE_GMP
   unif01_Gen * ugranger_CreateCombBigLCGCub (
      char *m1, char *a1, char *c1, char *s1, long m2, long a2, long c2);


 
   unif01_Gen * ugranger_CreateCombTausBigLCG (
      unsigned int k1, unsigned int q1, unsigned int s1, unsigned int SS1,
      unsigned int k2, unsigned int q2, unsigned int s2, unsigned int SS2,
      char *m, char *a, char *c, char *SS3);

#endif


unif01_Gen * ugranger_CreateCombTausLCG21xor (
   unsigned int k1, unsigned int q1, unsigned int s1, unsigned int SS1,
   unsigned int k2, unsigned int q2, unsigned int s2, unsigned int SS2,
   long m, long a, long c, long SS3);
 


unif01_Gen * ugranger_CreateCombTausCub21xor (
   unsigned int k1, unsigned int q1, unsigned int s1, unsigned int SS1,
   unsigned int k2, unsigned int q2, unsigned int s2, unsigned int SS2,
   long m, long a, long SS3);



unif01_Gen * ugranger_CreateCombTausInvExpl21xor (
   unsigned int k1, unsigned int q1, unsigned int s1, unsigned int SS1,
   unsigned int k2, unsigned int q2, unsigned int s2, unsigned int SS2,
   long m, long a, long c);


void ugranger_DeleteCombLCGInvExpl (unif01_Gen *gen);
void ugranger_DeleteCombLCGCub (unif01_Gen *gen);
void ugranger_DeleteCombTausLCG21xor (unif01_Gen *gen);
void ugranger_DeleteCombTausCub21xor (unif01_Gen *gen);
void ugranger_DeleteCombTausInvExpl21xor (unif01_Gen *gen);

#ifdef USE_GMP
   void ugranger_DeleteCombBigLCGInvExpl (unif01_Gen *gen);
   void ugranger_DeleteCombBigLCGCub (unif01_Gen *gen);
   void ugranger_DeleteCombTausBigLCG (unif01_Gen *gen);
#endif

 
#endif
 

