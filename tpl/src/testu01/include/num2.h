 
/* num2.h for ANSI C */

#ifndef NUM2_H
#define NUM2_H
 
#include "gdef.h"
#include <math.h>


double num2_Factorial (int n);



double num2_LnFactorial (int n);



double num2_Combination (int n, int s);



#ifdef HAVE_LGAMMA
#define num2_LnGamma lgamma
#else
   double num2_LnGamma (double x);
#endif



double num2_Digamma (double x);



#ifdef HAVE_LOG1P
#define num2_log1p log1p
#else
   double num2_log1p (double x);
#endif



void num2_CalcMatStirling (double *** M, int m, int n);



void num2_FreeMatStirling (double *** M, int m);



double num2_VolumeSphere (double p, int t);



double num2_EvalCheby (const double A[], int N, double x);



double num2_BesselK025 (double x);

 

#endif
 

