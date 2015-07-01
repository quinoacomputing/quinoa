
 
/* usoft.h for ANSI C */

#ifndef USOFT_H
#define USOFT_H
 
#include "gdef.h"
#include "unif01.h"


unif01_Gen * usoft_CreateSPlus (long S1, long S2);



#ifdef HAVE_RANDOM
   unif01_Gen * usoft_CreateUnixRandom (unsigned int s);
#endif



#ifdef USE_LONGLONG
   unif01_Gen * usoft_CreateJava48 (ulonglong s, int jflag);
#endif


unif01_Gen * usoft_CreateExcel97 (double r);



unif01_Gen * usoft_CreateExcel2003 (int x0, int y0, int z0);



unif01_Gen * usoft_CreateVisualBasic (unsigned long s);



#if defined(USE_GMP) && defined(USE_LONGLONG)
   unif01_Gen * usoft_CreateMaple_9 (longlong s);
#endif



#ifdef USE_LONGLONG
   unif01_Gen * usoft_CreateMATLAB (int i, unsigned int j, int bf,
                                    double Z[]);
#endif



#ifdef HAVE_MATHEMATICA
   unif01_Gen * usoft_CreateMathematicaReal (int argc, char * argv[],
                                             long s);
#endif



#ifdef HAVE_MATHEMATICA
   unif01_Gen * usoft_CreateMathematicaInteger (int argc, char * argv[],
                                                long s);
#endif


#ifdef USE_LONGLONG
   void usoft_DeleteMATLAB (unif01_Gen *gen);
#endif



#ifdef HAVE_MATHEMATICA
   void usoft_DeleteMathematicaReal (unif01_Gen *);
   void usoft_DeleteMathematicaInteger (unif01_Gen *);

#endif


#ifdef HAVE_RANDOM
   void usoft_DeleteUnixRandom (unif01_Gen *);

#endif


#if defined(USE_GMP) && defined(USE_LONGLONG)
   void usoft_DeleteMaple_9 (unif01_Gen *gen);

#endif


void usoft_DeleteGen (unif01_Gen *gen);

 
#endif
 

