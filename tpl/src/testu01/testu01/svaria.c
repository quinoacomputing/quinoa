/*************************************************************************\
 *
 * Package:        TestU01
 * File:           svaria.c
 * Environment:    ANSI C
 *
 * Copyright (c) 2002 Pierre L'Ecuyer, DIRO, Université de Montréal.
 * e-mail: lecuyer@iro.umontreal.ca
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted without a fee for private, research,
 * academic, or other non-commercial purposes.
 * Any use of this software in a commercial environment requires a
 * written licence from the copyright owner.
 *
 * Any changes made to this package must be clearly identified as such.
 *
 * In scientific publications which used this software, a reference to it
 * would be appreciated.
 *
 * Redistributions of source code must retain this copyright notice
 * and the following disclaimer.
 *
 * THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
\*************************************************************************/

#include "gdef.h"
#include "util.h"
#include "tables.h"
#include "chrono.h"
#include "num.h"
#include "num2.h"

#include "svaria.h"
#include "unif01.h"
#include "sres.h"
#include "wdist.h"
#include "swrite.h"
#include "smultin.h"

#include "fmass.h"
#include "gofs.h"
#include "gofw.h"

#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>




lebool svaria_Timer = FALSE;




/*-------------------------------- Constants ------------------------------*/

/* Max string lengths */
#define LEN1 100
#define LEN2 200

/* Sample limit for normal approximation in svaria_SampleMean */
#define SAM_LIM 60

/* Arrays dimension in svaria_AppearanceSpacings */
#define AS_DIM 32

/* The Meschach package for matrix computations */
#undef MESCHACH





/*-------------------------------- Functions ------------------------------*/


static void InitFDistMeans (int n, double Coef[])
/*
 * Initializes the distribution for svaria_SampleMean by computing the 
 * coefficients. We shall keep the value of n in element Coef[SAM_LIM],
 * since we shall need it to get the value of the distribution at x.
 */
{
   int s;
   double z;
   fmass_INFO Q;

   z = num2_Factorial (n);
   /* This uses the binomial formulae, but is not the binomial probability
      distribution since p + q != 1 */
   Q = fmass_CreateBinomial (n, -1.0, 1.0);
   for (s = 0; s <= n; s++)
      Coef[s] = fmass_BinomialTerm2 (Q, s) / z;
   fmass_DeleteBinomial (Q);
   Coef[SAM_LIM] = n;

   if (swrite_Classes) {
      printf ("---------------------------------------\n");
      for (s = 0; s <= n; s++) {
         printf ("   Coeff[%2d] = %14.6g\n", s, Coef[s]);
      }
      printf ("\n");
   }
}

/*-------------------------------------------------------------------------*/

static double FDistMeans (
   double C[],               /* Coefficients and sample size n */
   double x                  /* Argument */
   )
/* 
 * Distribution function of sample mean as in Stephens (1966), p.235.
 * This function is not very precise: the normal approximation is poor
 * for small x, and the computation of the exact function is numerically
 * unstable for large n. This will be used only for n < SAM_LIM. The
 * value of n is in  C[SAM_LIM].
 */
{
   double Sum;
   int M;
   int i;
   double nLR = C[SAM_LIM];
   int n = nLR;

   if (x <= 0.0)
      return 0.0;
   if (x >= n)
      return 1.0;

   M = x;
   Sum = 0.0;
   if (x < n / 2.0) {
      for (i = 0; i <= M; i++) {
         Sum += C[i] * pow (x, nLR);
         x -= 1.0;
      }
   } else {
      x = -x + nLR;
      for (i = n; i >= M + 1; i--) {
         Sum += C[i] * pow (x, nLR);
         x -= 1.0;
      }
      if (!(n & 1))
         Sum = -Sum;
      Sum += 1.0;
   }
   return Sum;
}

/*-------------------------------------------------------------------------*/

void svaria_SampleMean (unif01_Gen * gen, sres_Basic * res,
   long N, long n, int r)
{
   long i;
   long Seq;
   double Sum;
   double Coef[SAM_LIM + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "svaria_SampleMean test";

   Timer = chrono_Create ();
   if (swrite_Basic) {
      swrite_Head (gen, TestName, N, n, r);
      printf ("\n\n");
   }
   util_Assert (n > 1, "svaria_SampleMean:   n < 2");

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateBasic ();
   }
   sres_InitBasic (res, N, "svaria_SampleMean");
   if (n < SAM_LIM)
      InitFDistMeans (n, Coef);

   if (n < SAM_LIM)
      statcoll_SetDesc (res->sVal1, "SampleMean sVal1:   n*<U>");
   else
      statcoll_SetDesc (res->sVal1, "SampleMean sVal1:   standard normal");

   for (Seq = 1; Seq <= N; Seq++) {
      Sum = 0.0;
      for (i = 1; i <= n; i++)
         Sum += unif01_StripD (gen, r);

      if (n < SAM_LIM)
         statcoll_AddObs (res->sVal1, Sum);
      else
         statcoll_AddObs (res->sVal1, sqrt (12.0 / n) * (Sum - 0.5 * n));
   }

   if (n < SAM_LIM) {
      gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, FDistMeans, Coef,
                         res->sVal2, res->pVal2);
   } else {
      /* Normal approximation */
      gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_Normal,
         (double *) NULL, res->sVal2, res->pVal2);
   }
   res->pVal1->NObs = N;

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2,
         "Statistic value                       :");
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeleteBasic (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void svaria_SampleCorr (unif01_Gen * gen, sres_Basic * res,
   long N, long n, int r, int k)
{
   long i;
   long Seq;
   double U;
   double Sum;
   double *Pre;                   /* Previous k generated numbers */
   int pos;                       /* Circular index to element at lag k */
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "svaria_SampleCorr test";

   Timer = chrono_Create ();
   if (swrite_Basic) {
      swrite_Head (gen, TestName, N, n, r);
      printf (",   k = %d\n\n", k);
   }
   util_Assert (n > 2, "svaria_SampleCorr:   n <= 2");

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateBasic ();
   }
   sres_InitBasic (res, N, "svaria_SampleCorr");
   statcoll_SetDesc (res->sVal1,
      "SampleCorr sVal1:   asymptotic standard normal");

   Pre = util_Calloc ((size_t) (k + 1), sizeof (double));

   for (Seq = 1; Seq <= N; Seq++) {
      /* Generate first k numbers U and keep them in Pre */
      for (i = 0; i < k; i++)
         Pre[i] = unif01_StripD (gen, r);

      Sum = 0.0;
      pos = 0;
      /* Element Pre[pos] is at lag k from U */
      for (i = k; i < n; i++) {
         U = unif01_StripD (gen, r);
         Sum += Pre[pos] * U - 0.25;
         Pre[pos] = U;
         pos++;
         pos %= k;
      }
      /* Save standardized correlation */
      statcoll_AddObs (res->sVal1, Sum * sqrt (12.0 / (n - k)));
   }

   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_Normal,
       (double *) NULL, res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetNormalSumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2,
         "Normal statistic                      :");
      swrite_NormalSumTest (N, res);
      swrite_Final (gen, Timer);
   }
   util_Free (Pre);
   if (localRes)
      sres_DeleteBasic (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static double FDistProd (
   double Par[],             /* The parameter t = Par[0] */
   double x                  /* The argument */
   )
/*
 * Compute the distribution function F for the product of t random variables
 * U[0, 1], where 
 *                          t - 1
 *                           __            j
 *  F[u1*u2*...ut <= x] = x \      (-ln(x))
 *                          /__    -----------
 *                          j = 0      j!
 *
 */
{
   double vlog, jterm, vterm, Sum;
   int j, t;

   if (x >= 1.0)
      return 1.0;
   if (x <= 0.0)
      return 0.0;

   vlog = log (x);
   t = Par[0];
   Sum = 1.0;
   vterm = 1.0;
   jterm = 1.0;
   for (j = 1; j < t; j++) {
      vterm *= vlog;
      jterm *= -j;
      Sum += vterm / jterm;
      if (vterm / jterm < DBL_EPSILON)
         break;
   }
   return x * Sum;
}


/*-------------------------------------------------------------------------*/

void svaria_SampleProd (unif01_Gen * gen, sres_Basic * res,
   long N, long n, int r, int t)
{
   long i;
   int j;
   long Seq;
   double *P;
   double temp;
   double Par[1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "svaria_SampleProd test";

   Timer = chrono_Create ();
   if (swrite_Basic) {
      swrite_Head (gen, TestName, N, n, r);
      printf (",   t = %d\n\n", t);
   }

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateBasic ();
   }
   sres_InitBasic (res, N, "svaria_SampleProd");

   P = util_Calloc ((size_t) n + 1, sizeof (double));
   statcoll_SetDesc (res->sVal1, "SampleProd sVal1:   Uniform [0, 1]");
   Par[0] = t;

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 1; i <= n; i++) {
         temp = unif01_StripD (gen, r);
         for (j = 2; j <= t; j++)
            temp *= unif01_StripD (gen, r);
         P[i] = temp;
      }
      gofw_ActiveTests1 (P, n, FDistProd, Par, res->sVal2, res->pVal2);
      statcoll_AddObs (res->sVal1, res->pVal2[gofw_AD]);
   }

   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_Unif,
      (double *) NULL, res->sVal2, res->pVal2);
   res->pVal1->NObs = N;

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2,
         "Anderson-Darling statistic            :");
      swrite_Final (gen, Timer);
   }
   util_Free (P);
   if (localRes)
      sres_DeleteBasic (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void svaria_SumLogs (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r)
{
   const double Eps = DBL_EPSILON / 2.0;      /* To avoid log(0) */
   const double Epsilon = 1.E-100;            /* To avoid underflow */
   long i;
   long Seq;
   double u;
   double Prod;
   double Sum;
   double V[1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "svaria_SumLogs test";
   char chaine[LEN1 + 1] = "";
   char str[LEN2 + 1];

   Timer = chrono_Create ();
   if (swrite_Basic) {
      swrite_Head (gen, TestName, N, n, r);
      printf ("\n\n");
   }
   util_Assert (n < LONG_MAX/2, "2n > LONG_MAX");
   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, -1, "svaria_SumLogs");

   strncpy (chaine, "SumLogs sVal1:   chi2 with ", (size_t) LEN1);
   sprintf (str, "%ld", 2 * n);
   strncat (chaine, str, (size_t) LEN2);
   strncat (chaine, " degrees of freedom", (size_t) LEN1);
   statcoll_SetDesc (res->sVal1, chaine);
   res->degFree = 2 * n;
   if (res->degFree < 1) {
      util_Warning (TRUE, "Chi-square with 0 degree of freedom.");
      if (localRes)
         sres_DeleteChi2 (res);
      return;
   }

   for (Seq = 1; Seq <= N; Seq++) {
      Prod = 1.0;
      Sum = 0.0;
      for (i = 1; i <= n; i++) {
         u = unif01_StripD (gen, r);
         if (u < Eps)
            u = Eps;
         Prod *= u;
         if (Prod < Epsilon) {
            Sum += log (Prod);
            Prod = 1.0;
         }
      }
      statcoll_AddObs (res->sVal1, -2.0 * (Sum + log (Prod)));

   }
   V[0] = 2 * n;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
                      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetChi2SumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      swrite_AddStrChi (str, LEN2, res->degFree);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTest (N, res);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataWeight (unif01_Gen * gen, char *TestName,
   long N, long n, int r, long k, double Alpha, double Beta)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",  k = %1ld,  Alpha = %6.4g,  Beta = %6.4g\n\n",
           k, Alpha, Beta);
}


/*-------------------------------------------------------------------------*/

void svaria_WeightDistrib (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r, long k, double Alpha, double Beta)
{
   long W;
   long j;
   long i;
   long Seq;
   double X;
   double U;
   double p;
   double nLR = n;
   double V[1];
   long NbClasses;
   long *Loc;
   fmass_INFO Q;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "svaria_WeightDistrib test";
   char chaine[LEN1 + 1] = "";
   char str[LEN2 + 1];

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataWeight (gen, TestName, N, n, r, k, Alpha, Beta);

   /*   util_Assert (n >= 3.0 * gofs_MinExpected,
	"svaria_WeightDistrib:   n is too small"); */
   util_Assert (Alpha <= 1.0 && Alpha >= 0.0,
      "svaria_WeightDistrib:    Alpha must be in [0, 1]");
   util_Assert (Beta <= 1.0 && Beta >= 0.0,
      "svaria_WeightDistrib:    Beta must be in [0, 1]");
   p = Beta - Alpha;

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, k, "svaria_WeightDistrib");
   Loc = res->Loc;

   /* Compute binomial probabilities and multiply by n */
   Q = fmass_CreateBinomial (k, p, 1.0 - p);
   for (i = 0; i <= k; i++)
      res->NbExp[i] = nLR * fmass_BinomialTerm2 (Q, i);
   fmass_DeleteBinomial (Q);

   res->jmin = 0;
   res->jmax = k;
   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, Loc, res->jmin, res->jmax, 0);

   /* Merge classes for the chi-square */
   gofs_MergeClasses (res->NbExp, Loc, &res->jmin, &res->jmax, &NbClasses);

   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, Loc, res->jmin, res->jmax, NbClasses);

   strncpy (chaine, "WeightDistrib sVal1:   chi2 with ", (size_t) LEN1);
   sprintf (str, "%ld", NbClasses - 1);
   strncat (chaine, str, (size_t) LEN2);
   strncat (chaine, " degrees of freedom", (size_t) LEN1);
   statcoll_SetDesc (res->sVal1, chaine);
   res->degFree = NbClasses - 1;
   if (res->degFree < 1) {
      if (localRes)
         sres_DeleteChi2 (res);
      return;
   }

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 0; i <= k; i++)
         res->Count[i] = 0;
      for (i = 1; i <= n; i++) {
         W = 0;
         for (j = 1; j <= k; j++) {
            U = unif01_StripD (gen, r);
            if (U >= Alpha && U < Beta)
               ++W;
         }
         if (W > res->jmax)
            ++res->Count[res->jmax];
         else
            ++res->Count[Loc[W]];
      }
      if (swrite_Counters)
         tables_WriteTabL (res->Count, res->jmin, res->jmax, 5, 10,
                           "Observed numbers:");

      X = gofs_Chi2 (res->NbExp, res->Count, res->jmin, res->jmax);
      statcoll_AddObs (res->sVal1, X);
   }

   V[0] = NbClasses - 1;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
                      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetChi2SumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      swrite_AddStrChi (str, LEN2, res->degFree);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTest (N, res);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataArgMax (unif01_Gen * gen, char *TestName,
   long N, long n, int r, long k, long m)
{
   double x;

   swrite_Head (gen, TestName, N, n, r);
   printf (",   k = %1ld,   m = %1ld\n\n", k, m);
   printf ("   Number of balls = n = %1ld\n", n);
   printf ("   Number of urns  = k = %1ld\n", k);

   x = n;
   x = x * x / (2 * k);
   printf ("   Number (approx) of collisions = n^2 / 2k = %g\n\n\n", x);
}

/*-------------------------------------------------------------------------*/

static int svaria_CollisionArgMax_00 (unif01_Gen *gen, sres_Chi2 *res,
   long N, long n, int r, long k, long m)
/*
 * Return 0 if no error, otherwise return != 0.
 */
{
   double X;
   double U;
   double Max;
   long NbColl;
   long Indice = -1;
   long j;
   long i;
   long Rep;
   long Seq;
   long NbClasses;
   long *Loc;
   int *Urne;
   double V[1];
   fmass_INFO Q;
   lebool localRes = FALSE;
   chrono_Chrono *chro, *Timer;
   char *TestName = "svaria_CollisionArgMax test";
   char chaine[LEN1 + 1] = "";
   char str[LEN2 + 1];

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataArgMax (gen, TestName, N, n, r, k, m);

   util_Assert (n <= 4 * k, "svaria_CollisionArgMax:   n > 4k");
   /*   util_Assert (m > 2.0 * gofs_MinExpected,
	"svaria_CollisionArgMax:    m <= 2*gofs_MinExpected"); */

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, n, "svaria_CollisionArgMax");
   Loc = res->Loc;
   Urne = util_Calloc ((size_t) k + 1, sizeof (int));

   if (svaria_Timer) {
      printf ("-----------------------------------------------");
      printf ("\nCPU time to initialize the collision distribution:  ");
      chro = chrono_Create ();
   }
   Q = smultin_CreateCollisions (n, (smultin_CellType) k);
   if (svaria_Timer) {
      chrono_Write (chro, chrono_hms);
      printf ("\n\n");
   }

   /* Compute the expected numbers of collisions: m*P(j) */
   for (j = 0; j <= n; j++)
      res->NbExp[j] = m * smultin_CollisionsTerm (Q, j);
   smultin_DeleteCollisions (Q);

   res->jmin = 0;
   res->jmax = n;
   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, Loc, res->jmin, res->jmax, 0);

   gofs_MergeClasses (res->NbExp, Loc, &res->jmin, &res->jmax, &NbClasses);

   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, Loc, res->jmin, res->jmax, NbClasses);

   strncpy (chaine, "CollisionArgMax sVal1:   chi2 with ", (size_t) LEN1);
   sprintf (str, "%ld", NbClasses - 1);
   strncat (chaine, str, (size_t) LEN2);
   strncat (chaine, " degrees of freedom", (size_t) LEN1);
   statcoll_SetDesc (res->sVal1, chaine);
   res->degFree = NbClasses - 1;
   if (res->degFree < 1) {
      if (localRes)
         sres_DeleteChi2 (res);
      return 1;
   }

   if (svaria_Timer)
      chrono_Init (chro);

   for (Seq = 1; Seq <= N; Seq++) {
      for (j = 0; j <= n; j++)
         res->Count[j] = 0;

      for (Rep = 1; Rep <= m; Rep++) {
         for (j = 0; j <= k; j++)
            Urne[j] = -1;

         NbColl = 0;
         for (j = 1; j <= n; j++) {
            Max = -1.0;
            for (i = 1; i <= k; i++) {
               U = unif01_StripD (gen, r);
               if (U > Max) {
                  Max = U;
                  Indice = i;
               }
            }
            if (Urne[Indice] < 0)
               Urne[Indice] = 1;
            else
               ++NbColl;
         }
         if (NbColl > res->jmax)
            ++res->Count[res->jmax];
         else
            ++res->Count[Loc[NbColl]];
      }
      if (swrite_Counters)
         tables_WriteTabL (res->Count, res->jmin, res->jmax, 5, 10,
                           "Observed numbers:");
      X = gofs_Chi2 (res->NbExp, res->Count, res->jmin, res->jmax);
      statcoll_AddObs (res->sVal1, X);
   }

   if (svaria_Timer) {
      printf ("\n----------------------------------------------\n"
              "CPU time for the test           :  ");
      chrono_Write (chro, chrono_hms);
      printf ("\n\n");
      chrono_Delete (chro);
   }

   V[0] = NbClasses - 1;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
                      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetChi2SumStat (res);
   
   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      swrite_AddStrChi (str, LEN2, res->degFree);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTest (N, res);
      swrite_Final (gen, Timer);
   }
   util_Free (Urne);
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
   return 0;
}


/*-------------------------------------------------------------------------*/

void svaria_CollisionArgMax (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r, long k, long m)
{
   if (m > 1) {
      svaria_CollisionArgMax_00 (gen, res, N, n, r, k, m);

   } else if (m == 1) {
      double ValDelta[] = { -1.0 };
      smultin_Param *par;

      if (swrite_Basic) {
         printf (
          "***********************************************************\n"
          "Test svaria_CollisionArgMax calling smultin_Multinomial\n\n");
      }
      par = smultin_CreateParam (1, ValDelta, smultin_GenerCellMax, -3);
      if (NULL == res) {
         smultin_Multinomial (gen, par, NULL, N, n, r, 0, k, TRUE);
      } else {
         smultin_Res *resm;
         resm = smultin_CreateRes (par);
         smultin_Multinomial (gen, par, resm, N, n, r, 0, k, TRUE);
         sres_InitChi2 (res, N, -1, "svaria_CollisionArgMax");
         statcoll_SetDesc (res->sVal1, "CollisionArgMax sVal1");
         res->sVal1->NObs = resm->Collector[0]->NObs;
         tables_CopyTabD (resm->Collector[0]->V, res->sVal1->V, 1, N);
         tables_CopyTabD (resm->sVal2[0], res->sVal2, 0, gofw_NTestTypes - 1);
         tables_CopyTabD (resm->pVal2[0], res->pVal2, 0, gofw_NTestTypes - 1);
         smultin_DeleteRes (resm);
      }
      smultin_DeleteParam (par);
   } else {
     util_Warning (m <= 0, "svaria_CollisionArgMax:   m <= 0");
   }
}


/*=========================================================================*/

static void WriteDataSumColl (unif01_Gen * gen, char *TestName,
   long N, long n, int r, double g)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   g = %g\n\n", g);
}

/*-------------------------------------------------------------------------*/

static double ProbabiliteG (int jmin, int j, double g)
/* 
 * Returns the probability that the minimum number of random U(0,1) whose
 * sum is larger than g is j+1. g cannot be too large because the
 * calculation here becomes numerically unstable.
 */
{
   int s;
   double temp;
   const double jLR = j;
   double signe;                  /* +1 or -1 */
   double somme;

   signe = 1.0;
   somme = 0.0;
   for (s = 0; s <= jmin; s++) {
      temp = signe * num2_Combination (j + 1, s);
      temp *= pow (g - s, jLR);
      somme += temp;
      signe = -signe;
   }
   somme = (jLR + 1.0 - g) * somme / num2_Factorial (j + 1);
   return somme;
}


/*-------------------------------------------------------------------------*/

void svaria_SumCollector (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r, double g)
{
   const double gmax = 10.0;      /* Maximal value of g */
   const int jmax = 50;           /* Maximal number of classes */
   int j;                         /* Class index */
   long Seq;
   long i;
   double X;
   double Y;
   double Sum;
   long NbClasses;
   long *Loc;
   double V[1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "svaria_SumCollector test";
   char chaine[LEN1 + 1] = "";
   char str[LEN2 + 1];

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataSumColl (gen, TestName, N, n, r, g);

   if (g < 1.0 || g > gmax) {
      util_Error ("svaria_SumCollector:   g < 1.0 or g > 10.0");
   }
   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, jmax, "svaria_SumCollector");
   Loc = res->Loc;

   res->jmin = g;
   res->jmax = jmax;
   Sum = 0.0;
   for (j = res->jmin; j < jmax; j++) {
      res->NbExp[j] = n * ProbabiliteG (res->jmin, j, g);
      Sum += res->NbExp[j];
   }
   res->NbExp[jmax] = util_Max (0.0, n - Sum);

   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, Loc, res->jmin, res->jmax, 0);
   gofs_MergeClasses (res->NbExp, Loc, &res->jmin, &res->jmax, &NbClasses);
   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, Loc, res->jmin, res->jmax, NbClasses);

   strncpy (chaine, "SumCollector sVal1:   chi2 with ", (size_t) LEN1);
   sprintf (str, "%ld", NbClasses - 1);
   strncat (chaine, str, (size_t) LEN2);
   strncat (chaine, " degrees of freedom", (size_t) LEN1);
   statcoll_SetDesc (res->sVal1, chaine);
   res->degFree = NbClasses - 1;
   if (res->degFree < 1) {
      if (localRes)
         sres_DeleteChi2 (res);
      return;
   }

   for (Seq = 1; Seq <= N; Seq++) {
      for (j = 1; j <= jmax; j++)
         res->Count[j] = 0;

      for (i = 1; i <= n; i++) {
         X = 0.0;
         j = 0;
         do {
            X += unif01_StripD (gen, r);
            ++j;
         }
         while (X <= g);
         if (j > res->jmax)
            ++res->Count[res->jmax];
         else
            ++res->Count[Loc[j - 1]];
      }
      if (swrite_Counters)
         tables_WriteTabL (res->Count, res->jmin, res->jmax, 5, 10,
                           "Observed numbers:");
      Y = gofs_Chi2 (res->NbExp, res->Count, res->jmin, res->jmax);
      statcoll_AddObs (res->sVal1, Y);
   }

   V[0] = NbClasses - 1;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
                      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetChi2SumStat (res);
   
   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      swrite_AddStrChi (str, LEN2, res->degFree);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTest (N, res);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void InitAppear (int r, int s, int L, long Q, double E[], double KV[])
{
   util_Assert (r >= 0, "svaria_AppearanceSpacings:   r < 0");
   util_Assert (s > 0, "svaria_AppearanceSpacings:   s <= 0");
   /*   if (L >= s && L % s) {
      util_Error ("svaria_AppearanceSpacings:   L mod s != 0");
      }*/
   if (L < s && s % L) {
      util_Error ("svaria_AppearanceSpacings:   s mod L != 0");
   }
   util_Warning (Q < 10.0 * num_TwoExp[L],
      "svaria_AppearanceSpacings:   Q < 10 * 2^L");

   /* Theoretical mean E and variance KV for different L given by Maurer */
   if (L > 16) {
      /* For L > 16 and near 16, the 6-th decimal of E[L] could be
         erroneous. For KV [L], the 4-th decimal. */
      E[L] = L - 8.32746E-1;
      KV[L] = 3.423715;
   } else {
      E [1]  = 0.73264948;      KV[1]  = 0.68977;
      E [2]  = 1.53743829;      KV[2]  = 1.33774;
      E [3]  = 2.40160681;      KV[3]  = 1.90133;
      E [4]  = 3.31122472;      KV[4]  = 2.35774;
      E [5]  = 4.25342659;      KV[5]  = 2.70455;
      E [6]  = 5.21770525;      KV[6]  = 2.95403;
      E [7]  = 6.19625065;      KV[7]  = 3.12539;
      E [8]  = 7.18366555;      KV[8]  = 3.23866;
      E [9]  = 8.17642476;      KV[9]  = 3.31120;
      E [10] = 9.17232431;      KV[10] = 3.35646;
      E [11] = 10.1700323;      KV[11] = 3.38409;
      E [12] = 11.1687649;      KV[12] = 3.40065;
      E [13] = 12.1680703;      KV[13] = 3.41043;
      E [14] = 13.1676926;      KV[14] = 3.41614;
      E [15] = 14.1674884;      KV[15] = 3.41943;
      E [16] = 15.1673788;      KV[16] = 3.42130;
   }
}

/*-------------------------------------------------------------------------*/

static double CalcSigma (int L, long K, double KV[])
/*
 * Compute the standard deviation for the svaria_AppearanceSpacings test.
 */
{
   double dCor[AS_DIM + 1];        /* Coron-Naccache factor d */
   double eCor[AS_DIM + 1];        /* Coron-Naccache factor e */
   double temp;
   /*   double c; */

#if 0     /* No correction */
   return sqrt (KV[L] / K);

#elif 0   /* Maurer's correction c(L, K) */
   temp = 3.0 / L * num_Log2 ((double) K);
   if (temp >= DBL_MAX_EXP - 1)
      temp = 0.0;
   else
      temp = pow (2.0, -temp);
   c = 0.7 - 0.8/L + (4.0 + 32.0/L) * temp / 15.0;
   if (L < 3 || L > 16)
      c = 1.0;
   return c * sqrt (KV[L] / K);

#else   /* based on Coron and Naccache exact calculation */
   dCor [3]  = 0.2732725;        eCor [3]  = 0.4890883;
   dCor [4]  = 0.3045101;        eCor [4]  = 0.4435381;
   dCor [5]  = 0.3296587;        eCor [5]  = 0.4137196;
   dCor [6]  = 0.3489769;        eCor [6]  = 0.3941338;
   dCor [7]  = 0.3631815;        eCor [7]  = 0.3813210;
   dCor [8]  = 0.3732189;        eCor [8]  = 0.3730195;
   dCor [9]  = 0.3800637;        eCor [9]  = 0.3677118;
   dCor [10] = 0.3845867;        eCor [10] = 0.3643695;
   dCor [11] = 0.3874942;        eCor [11] = 0.3622979;
   dCor [12] = 0.3893189;        eCor [12] = 0.3610336;
   dCor [13] = 0.3904405;        eCor [13] = 0.3602731;
   dCor [14] = 0.3911178;        eCor [14] = 0.3598216;
   dCor [15] = 0.3915202;        eCor [15] = 0.3595571;
   dCor [16] = 0.3917561;        eCor [16] = 0.3594040;
   /* L = infinite */
   dCor [0]  = 0.3920729;        eCor [0]  = 0.3592016;

   if (L < 3)
      return sqrt (KV[L] / K);
   if (L > 16)
      temp = dCor[0] + eCor [0] * num_TwoExp[L] / K;
   else
      temp = dCor[L] + eCor [L] * num_TwoExp[L] / K;
   return sqrt (temp * KV[L] / K);
   
#endif
}

/*-------------------------------------------------------------------------*/

static void WriteDataAppear (unif01_Gen * gen,
   long N, int r, int s, int L, long Q, long K, double n)
{
   printf ("***********************************************************\n");
   printf ("HOST = ");
   if (swrite_Host) {
      gdef_WriteHostName ();
      printf ("\n");
   } else 
      printf ("\n\n");
   unif01_WriteNameGen (gen);
   printf ("\n");
   if (swrite_ExperimentName && strcmp (swrite_ExperimentName, "")) {
      printf ("%s", swrite_ExperimentName);
      printf (":\n\n");
   }

   printf ("svaria_AppearanceSpacings test:\n"
          "-----------------------------------------------\n");

   printf ("   N = %2ld,   Q = %1ld,   K = %1ld,   r = %1d,   s = %1d,"
           "   L = %1d\n\n", N, Q, K, r, s, L);
   printf ("   Sequences of n = (K + Q)L = %12.0f bits\n", n);
   printf ("   Q = %4ld initialization blocks\n", Q);
   printf ("   K = %4ld blocks for the test\n", K);
   printf ("   the blocks have L = %2d bits\n\n\n", L);
}


/*-------------------------------------------------------------------------*/

void svaria_AppearanceSpacings (unif01_Gen * gen, sres_Basic * res,
   long N, long Q, long K, int r, int s, int L)
{
   double E[AS_DIM + 1];          /* Theoretical mean of the log (Base2) of
                                     the most recent occurrence of a block */
   double KV[AS_DIM + 1];         /* K times the theoretical variance of the
                                     same */
   long Seq;
   long block;
   long Nblocks;                  /* 2^L = total number of distinct blocks */
   long K2;
   long Q2;
   long i;
   long sBits;                    /* Numerical value of the s given bits */
   long d;                        /* 2^s */
   const int SdivL = s / L;
   const int LdivS = L / s;
   const int LmodS = L % s;
   long sd;                       /* 2^LmodS */
   long rang;
   double n;                      /* Total number of bits in a sequence */
   double sigma;                  /* Standard deviation = sqrt (Variance) */
   double somme;
   double ARang;                  /* Most recent occurrence of block */
   long *Count;                   /* Index of most recent occurrence of
                                     block */
   double FactMoy;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;

   Timer = chrono_Create ();
   n = ((double) K + (double) Q) * L;
   if (swrite_Basic)
      WriteDataAppear (gen, N, r, s, L, Q, K, n);
   util_Assert (s < 32, "svaria_AppearanceSpacings:   s >= 32");
   InitAppear (r, s, L, Q, E, KV);
   sigma = CalcSigma (L, K, KV);
   d = num_TwoExp[s];
   Nblocks = num_TwoExp[L];
   FactMoy = 1.0 / (num_Ln2 * K);

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateBasic ();
   }
   sres_InitBasic (res, N, "svaria_AppearanceSpacings");
   Count = util_Calloc ((size_t) Nblocks + 2, sizeof (long));

   statcoll_SetDesc (res->sVal1,
      "AppearanceSpacings sVal1:   standard normal");

   if (LdivS > 0) {
      sd = num_TwoExp[LmodS];

      for (Seq = 1; Seq <= N; Seq++) {
         for (i = 0; i < Nblocks; i++)
            Count[i] = 0;

         /* Initialization with Q blocks */
         for (rang = 0; rang < Q; rang++) {
            block = 0;
            for (i = 1; i <= LdivS; i++) {
               sBits = unif01_StripB (gen, r, s);
               block = block * d + sBits;
            }
            if (LmodS > 0) {
               sBits = unif01_StripB (gen, r, LmodS);
               block = block * sd + sBits;
            }
            Count[block] = rang;
         }

         /* Test proper with K blocks */
         somme = 0.0;
         for (rang = Q; rang < Q + K; rang++) {
            block = 0;
            for (i = 1; i <= LdivS; i++) {
               sBits = unif01_StripB (gen, r, s);
               block = block * d + sBits;
            }
            if (LmodS > 0) {
               sBits = unif01_StripB (gen, r, LmodS);
               block = block * sd + sBits;
            }
            ARang = rang - Count[block];
            somme += log (ARang);
            Count[block] = rang;
         }
         statcoll_AddObs (res->sVal1, (somme * FactMoy - E[L]) / sigma);
      }

   } else {                       /* s > L */
      Q2 = Q / SdivL;
      K2 = K / SdivL;
      for (Seq = 1; Seq <= N; Seq++) {
         for (i = 0; i < Nblocks; i++)
            Count[i] = 0;

         /* Initialization: Q blocks */
         for (rang = 0; rang < Q2; rang++) {
            sBits = unif01_StripB (gen, r, s);
            for (i = 0; i < SdivL; i++) {
               block = sBits % Nblocks;
               Count[block] = SdivL * rang + i;
               sBits /= Nblocks;
            }
         }
         /* Test proper with K blocks */
         somme = 0.0;
         for (rang = Q2; rang < Q2 + K2; rang++) {
            sBits = unif01_StripB (gen, r, s);
            for (i = 0; i < SdivL; i++) {
               block = sBits % Nblocks;
               ARang = SdivL * rang + i - Count[block];
               somme += log (ARang);
               Count[block] = SdivL * rang + i;
               sBits /= Nblocks;
            }
         }
         statcoll_AddObs (res->sVal1, (somme * FactMoy - E[L]) / sigma);
      }
   }

   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_Normal,
      (double *) NULL, res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetNormalSumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 12, 4, 3);

   if (swrite_Basic) {
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2,
         "Normal statistic                      :");
      swrite_NormalSumTest (N, res);
      swrite_Final (gen, Timer);
   }
   util_Free (Count);
   if (localRes)
      sres_DeleteBasic (res);
   chrono_Delete (Timer);
}
