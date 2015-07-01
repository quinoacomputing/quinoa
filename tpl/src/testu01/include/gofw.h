 
/* gofw.h for ANSI C */
#ifndef GOFW_H
#define GOFW_H
 
#include "gdef.h"           /* From the library mylib */
#include "bitset.h"         /* From the library mylib */
#include "fdist.h"
#include "wdist.h"
#include <stdio.h>


typedef enum {
   gofw_Gnuplot,
   gofw_Mathematica
   } gofw_GraphType;



extern gofw_GraphType gofw_GraphSoft;



void gofw_GraphFunc (FILE *f, wdist_CFUNC F, double par[], double a,
                     double b, int m, int mono, char Desc[]);



void gofw_GraphDistUnif (FILE *f, double U[], long N, char Desc[]);


extern double gofw_Epsilonp;
extern double gofw_Epsilonp1;



extern double gofw_Suspectp;



double gofw_pDisc (double pL, double pR);



void gofw_Writep0 (double p);



void gofw_Writep1 (double p);



void gofw_Writep2 (double x, double p);



void gofw_WriteKS0 (long N, double DP, double DM, double D);



void gofw_WriteKS1 (double V[], long N, wdist_CFUNC F, double par[]);



void gofw_WriteKSJumpOne0 (long N, double a, double DP);



void gofw_WriteKSJumpOne1 (double V[], long N, 
                           wdist_CFUNC F, double par[], double a);

#if 0

void gofw_KSJumpsMany0 (double DP, double DM, fdist_FUNC_JUMPS *H);



void gofw_KSJumpsMany2 (statcoll_Collector *S, fdist_FUNC_JUMPS *H,
                        int Detail);



#endif


typedef enum {
   gofw_KSP,                      /* Kolmogorov-Smirnov+        */
   gofw_KSM,                      /* Kolmogorov-Smirnov-        */
   gofw_KS,                       /* Kolmogorov-Smirnov         */
   gofw_AD,                       /* Anderson-Darling           */
   gofw_CM,                       /* Cramer-vonMises            */
   gofw_WG,                       /* Watson G                   */
   gofw_WU,                       /* Watson U                   */
   gofw_Mean,                     /* Mean                       */
   gofw_Var,                      /* Variance                   */
   gofw_Cor,                      /* Correlation                */
   gofw_Sum,                      /* Sum                        */
   gofw_NTestTypes                /* Total number of test types */
   } gofw_TestType;



typedef double gofw_TestArray [gofw_NTestTypes];



extern char *gofw_TestNames [gofw_NTestTypes];



extern bitset_BitSet gofw_ActiveTests;



void gofw_InitTestArray (gofw_TestArray A, double x);



void gofw_Tests0 (double U[], long N, gofw_TestArray sVal);



void gofw_Tests1 (double V[], long N, wdist_CFUNC F, double par[],
                  gofw_TestArray sVal);



void gofw_ActiveTests0 (double U[], long N, 
                        gofw_TestArray sVal, gofw_TestArray pVal);



void gofw_ActiveTests1 (double V[], long N, wdist_CFUNC F, double par[],
                        gofw_TestArray sVal, gofw_TestArray pVal);



void gofw_ActiveTests2 (double V[], double U[], long N, wdist_CFUNC F,
                        double par[], gofw_TestArray sVal,
                        gofw_TestArray pVal);



void gofw_WriteActiveTests0 (long N, gofw_TestArray sVal,
                                     gofw_TestArray pVal);



void gofw_WriteActiveTests1 (double V[], long N, 
                             wdist_CFUNC F, double par[]);



void gofw_WriteActiveTests2 (long N, gofw_TestArray sVal,
                             gofw_TestArray pVal, char Desc[]);



void gofw_IterSpacingsTests0 (double U[], long N, int k, 
                              lebool printval, lebool graph, FILE *f);



void gofw_IterPowRatioTests0 (double U[], long N, int k,
                              lebool printval, lebool graph, FILE *f);

 
#endif
 

