 
/* gofs.h for ANSI C */
#ifndef GOFS_H
#define GOFS_H
 
#include "bitset.h"       /* From the library mylib */
#include "fmass.h"
#include "fdist.h"
#include "wdist.h"



extern double gofs_MinExpected;



extern double gofs_EpsilonAD;


void gofs_ContUnifTransform (double V[], long N, wdist_CFUNC F,
                             double par[], double U[]);



void gofs_DiscUnifTransform (double V[], long N, wdist_DFUNC F,
                             fmass_INFO W, double U[]);



void gofs_DiffD (double U[], double D[], long N1, long N2, 
                 double a, double b);



void gofs_DiffL (long U[], long D[], long N1, long N2, long a, long b);

#ifdef USE_LONGLONG
void gofs_DiffLL (longlong U[], longlong D[], long N1, long N2,
                  longlong a, longlong b);
void gofs_DiffULL (ulonglong U[], ulonglong D[], long N1, long N2,
                   ulonglong a, ulonglong b);
#endif



void gofs_IterateSpacings (double V[], double S[], long N);



void gofs_PowerRatios (double U[], long N);



void gofs_MergeClasses (double NbExp[], long Loc[],
                        long *smin, long *smax, long *NbClasses);



void gofs_WriteClasses (double NbExp[], long Loc[], 
                        long smin, long smax, long NbClasses);


double gofs_Chi2 (double NbExp[], long Count[], long smin, long smax);



double gofs_Chi2Equal (double NbExp, long Count[], long smin, long smax);



long gofs_Scan (double U[], long N, double d);



double gofs_CramerMises (double U[], long N);



double gofs_WatsonG (double U[], long N);



double gofs_WatsonU (double U[], long N);



double gofs_AndersonDarling (double U[], long N);



void gofs_KS (double U[], long N, double *DP, double *DM, double *D);



void gofs_KSJumpOne (double U[], long N, double a, double *DP, double *DM);

#if 0
void gofs_KSJumpsMany (double X[], int N, wdist_CFUNC F, double W[],
                       double *DP, double *DM, int Detail);
#endif

 
#endif
 

