
 
/*  ulec.h  for ANSI C  */

#ifndef ULEC_H
#define ULEC_H
 
#include "gdef.h"
#include "unif01.h"


unif01_Gen * ulec_CreateCombLec88 (long S1, long S2);



unif01_Gen * ulec_CreateCombLec88Float (long S1, long S2);



unif01_Gen * ulec_CreateCLCG4 (long S1, long S2, long S3, long S4);



unif01_Gen * ulec_CreateCLCG4Float (long S1, long S2, long S3, long S4);



unif01_Gen * ulec_CreateMRG93 (long S1, long S2, long S3, long S4, long S5);



unif01_Gen * ulec_CreateCombMRG96 (long S11, long S12, long S13, 
                                   long S21, long S22, long S23);



unif01_Gen * ulec_CreateCombMRG96Float (long S11, long S12, long S13,
                                        long S21, long S22, long S23);



unif01_Gen * ulec_CreateCombMRG96D (long S11, long S12, long S13,
                                    long S21, long S22, long S23);



unif01_Gen * ulec_CreateCombMRG96FloatD (long S11, long S12, long S13,
                                         long S21, long S22, long S23);



unif01_Gen * ulec_CreateMRG32k3a (double x10, double x11, double x12,
                                  double x20, double x21, double x22);



unif01_Gen * ulec_CreateMRG32k3aL (long x10, long x11, long x12,
                                   long x20, long x21, long x22);



unif01_Gen * ulec_CreateMRG32k3b (double x10, double x11, double x12,
                                  double x20, double x21, double x22);



unif01_Gen * ulec_CreateMRG32k5a (double x10, double x11, double x12,
                                  double x13, double x14, double x20,
                                  double x21, double x22, double x23,
                                  double x24);



unif01_Gen * ulec_CreateMRG32k5b (double x10, double x11, double x12,
                                  double x13, double x14, double x20,
                                  double x21, double x22, double x23,
                                  double x24);



#ifdef USE_LONGLONG
unif01_Gen * ulec_CreateMRG63k3a (longlong s10, longlong s11, longlong s12,
                                  longlong s20, longlong s21, longlong s22);



unif01_Gen * ulec_CreateMRG63k3b (longlong s10, longlong s11, longlong s12,
                                  longlong s20, longlong s21, longlong s22);

#endif


unif01_Gen * ulec_Createlfsr88 (unsigned int s1, unsigned int s2,
                                unsigned int s3);



unif01_Gen * ulec_Createlfsr113 (unsigned int s1, unsigned int s2,
                                 unsigned int s3, unsigned int s4);



#ifdef USE_LONGLONG
unif01_Gen * ulec_Createlfsr258 (ulonglong s1, ulonglong s2, ulonglong s3,
                                 ulonglong s4, ulonglong s5);
#endif



unif01_Gen * ulec_CreateCombTausLCG11 (unsigned int k, unsigned  int q, 
                                       unsigned int s, unsigned S1,
                                       long m, long a, long c, long S2);



unif01_Gen * ulec_CreateCombTausLCG21 (unsigned int k1, unsigned int q1,
                                       unsigned int s1, unsigned int Y1, 
                                       unsigned int k2, unsigned int q2,
                                       unsigned int s2, unsigned int Y2, 
                                       long m, long a, long c, long Y3);



unif01_Gen * ulec_CreateMRG31k3p (long x10, long x11, long x12,
                                  long x20, long x21, long x22);

void ulec_DeleteCombTausLCG11 (unif01_Gen *gen);



void ulec_DeleteCombTausLCG21 (unif01_Gen *gen);



void ulec_DeleteGen (unif01_Gen *gen);
   
 
#endif
 

