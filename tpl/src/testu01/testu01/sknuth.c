/*************************************************************************\
 *
 * Package:        TestU01
 * File:           sknuth.c
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

#include "util.h"
#include "tables.h"
#include "chrono.h"
#include "num2.h"

#include "sknuth.h"
#include "sres.h"
#include "smultin.h"
#include "wdist.h"
#include "swrite.h"
#include "unif01.h"

#include "gofs.h"
#include "gofw.h"

#include <math.h>
#include <stdio.h>


#define LENGTH 200




/*-------------------------------- Functions ------------------------------*/


static void InitRes1 (
   sknuth_Res1 *res,          /* Results holder */
   long N,                    /* Number of replications */
   int d                      /* Max class index for chi2 */
)
/* 
 * Initializes the sknuth_Res structure
 */
{
   sres_InitBasic (res->Bas, N, "sknuth_MaxOft:   Anderson-Darling");
   sres_InitChi2 (res->Chi, N, d, "sknuth_MaxOft:   Chi2");
}


/*-------------------------------------------------------------------------*/

sknuth_Res1 * sknuth_CreateRes1 (void)
{
   sknuth_Res1 *res;
   res = util_Malloc (sizeof (sknuth_Res1));
   res->Bas = sres_CreateBasic ();
   res->Chi = sres_CreateChi2 ();
   return res;
}


/*-------------------------------------------------------------------------*/

void sknuth_DeleteRes1 (sknuth_Res1 *res)
{
   if (res == NULL)
      return;
   sres_DeleteBasic (res->Bas);
   sres_DeleteChi2 (res->Chi);
   util_Free (res);
}


/*=========================================================================*/

static void InitRes2 (
   sknuth_Res2 *res,           /* Results holder */
   long N,                     /* Number of replications */
   double Lambda,              /* Poisson mean */
   char *nam                   /* Test name */
)
/* 
 * Initializes res
 */
{
   sres_InitBasic (res->Bas, N, nam);
   sres_InitPoisson (res->Pois, N, Lambda, nam);
}


/*-------------------------------------------------------------------------*/

sknuth_Res2 * sknuth_CreateRes2 (void)
{
   sknuth_Res2 *res;
   res = util_Malloc (sizeof (sknuth_Res2));
   res->Bas = sres_CreateBasic ();
   res->Pois = sres_CreatePoisson ();
   res->Pois->pLeft = -1.0;
   res->Pois->pRight = -1.0;
   return res;
}


/*-------------------------------------------------------------------------*/

void sknuth_DeleteRes2 (sknuth_Res2 *res)
{
   if (res == NULL)
      return;
   sres_DeleteBasic (res->Bas);
   sres_DeletePoisson (res->Pois);
   util_Free (res);
}


/*=========================================================================*/

void sknuth_Serial (unif01_Gen *gen, sres_Chi2 *res,
                    long N, long n, int r, long d, int t)
{
   double ValDelta[] = { 1.0 };
   smultin_Param *par;

   if (swrite_Basic)
      printf ("***********************************************************\n"
              "Test sknuth_Serial calling smultin_Multinomial\n\n");

   par = smultin_CreateParam (1, ValDelta, smultin_GenerCellSerial, 3);
   if (NULL == res) {
      smultin_Multinomial (gen, par, NULL, N, n, r, d, t, FALSE);
   } else {
      smultin_Res *resm;
      resm = smultin_CreateRes (par);
      smultin_Multinomial (gen, par, resm, N, n, r, d, t, FALSE);
      sres_InitChi2 (res, N, -1, "sknuth_Serial");
      statcoll_SetDesc (res->sVal1, "Serial sVal1");
      res->sVal1->NObs = resm->Collector[0]->NObs;
      tables_CopyTabD (resm->Collector[0]->V, res->sVal1->V, 1, N);
      tables_CopyTabD (resm->sVal2[0], res->sVal2, 0, gofw_NTestTypes - 1);
      tables_CopyTabD (resm->pVal2[0], res->pVal2, 0, gofw_NTestTypes - 1);
      smultin_DeleteRes (resm);
   }
   smultin_DeleteParam (par);
}


/*=========================================================================*/

void sknuth_SerialSparse (unif01_Gen *gen, sres_Chi2 *res,
                          long N, long n, int r, long d, int t)
{
   double ValDelta[] = { 1.0 };
   smultin_Param *par;

   if (swrite_Basic)
      printf ("***********************************************************\n"
              "Test sknuth_SerialSparse calling smultin_Multinomial\n\n");

   par = smultin_CreateParam (1, ValDelta, smultin_GenerCellSerial, 3);
   if (NULL == res) {
      smultin_Multinomial (gen, par, NULL, N, n, r, d, t, TRUE);
   } else {
      smultin_Res *resm;
      resm = smultin_CreateRes (par);
      smultin_Multinomial (gen, par, resm, N, n, r, d, t, TRUE);
      sres_InitChi2 (res, N, -1, "sknuth_SerialSparse");
      statcoll_SetDesc (res->sVal1, "Serial sVal1");
      res->sVal1->NObs = resm->Collector[0]->NObs;
      tables_CopyTabD (resm->Collector[0]->V, res->sVal1->V, 1, N);
      tables_CopyTabD (resm->sVal2[0], res->sVal2, 0, gofw_NTestTypes - 1);
      tables_CopyTabD (resm->pVal2[0], res->pVal2, 0, gofw_NTestTypes - 1);
      smultin_DeleteRes (resm);
   }
   smultin_DeleteParam (par);
}


/*=========================================================================*/

static void WriteDataGap (unif01_Gen *gen, char *TestName,
   long N, long n, int r, double Alpha, double Beta)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   Alpha = %8.6g,   Beta  = %8.6g\n\n", Alpha, Beta);
}


/*-------------------------------------------------------------------------*/

void sknuth_Gap (unif01_Gen *gen, sres_Chi2 *res,
                 long N, long n, int r, double Alpha, double Beta)
{
   int len;
   int t;
   long m;                        /* Number of observed Gaps */
   long Seq;                      /* Current replication number */
   double p;                      /* Probability of U01 in (Alpha, Beta) */
   double X2;
   double U;
   double Mult;
   double V[1];                   /* Number of degrees of freedom for Chi2 */
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sknuth_Gap test";

   Timer = chrono_Create ();
   p = Beta - Alpha;
   t = log (gofs_MinExpected / n) / num2_log1p (-p);
   len = 1 + log (gofs_MinExpected / (n*p)) / num2_log1p (-p);
   t = util_Min(t, len);
   t = util_Max(t, 0);

   Mult = p * n;
   if (swrite_Basic)
      WriteDataGap (gen, TestName, N, n, r, Alpha, Beta);

   util_Assert (Alpha >= 0.0 && Alpha <= 1.0,
                "sknuth_Gap:   Alpha outside interval [0..1]");
   util_Assert (Beta <= 1.0 && Beta > Alpha,
                "sknuth_Gap:   Beta outside interval (Alpha..1]");

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, t, "sknuth_Gap");

   sprintf (str, "The N statistic values (a ChiSquare with %1d degrees"
                 " of freedom):", t);
   statcoll_SetDesc (res->sVal1, str);
   res->degFree = t;
   if (res->degFree < 1) {
      util_Warning (TRUE, "Chi-square with 0 degree of freedom.");
      if (localRes)
         sres_DeleteChi2 (res);
      chrono_Delete (Timer);
      return;
   }

   /* Compute the probabilities for each gap length */
   res->NbExp[0] = Mult;
   res->Loc[0] = 0;
   for (len = 1; len < t; len++) {
      Mult *= 1.0 - p;
      res->NbExp[len] = Mult;
      res->Loc[len] = len;
   }
   res->NbExp[t] = Mult * (1.0 - p) / p;
   res->Loc[t] = t;
   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, res->Count, 0, t, 0);

   /* Beginning of test */
   for (Seq = 1; Seq <= N; Seq++) {
      for (len = 0; len <= t; len++)
         res->Count[len] = 0;
      for (m = 1; m <= n; m++) {
         /* Process one gap */
         len = 0;
         U = unif01_StripD (gen, r);
         while ((U < Alpha || U >= Beta) && len < n) {
            ++len;
            U = unif01_StripD (gen, r);
         }
         if (len >= n) {
            util_Warning (TRUE,
   "sknuth_Gap:   one gap of length > n\n*********  Interrupting the test\n");
            printf ("\n\n");
            res->pVal2[gofw_Mean] = res->pVal2[gofw_AD]
                   = res->pVal2[gofw_KSM] = res->pVal2[gofw_KSP] = 0.0;
            if (localRes)
               sres_DeleteChi2 (res);
            chrono_Delete (Timer);
            return;
         }
         if (len >= t)
            ++res->Count[t];
         else
            ++res->Count[len];
      }
      if (swrite_Counters)
         tables_WriteTabL (res->Count, 0, t, 5, 10, "Observed numbers:");

      X2 = gofs_Chi2 (res->NbExp, res->Count, 0, t);
      statcoll_AddObs (res->sVal1, X2);
   }

   V[0] = t;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
                      res->sVal2, res->pVal2);
   sres_GetChi2SumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      swrite_AddStrChi (str, LENGTH, res->degFree);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTest (N, res);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataPoker (unif01_Gen * gen, char *TestName,
   long N, long n, int r, int d, int k)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   d = %4d,   k = %4d\n\n", d, k);
}


/*-------------------------------------------------------------------------*/

#define Maxkd 127

void sknuth_SimpPoker (unif01_Gen *gen, sres_Chi2 *res,
                       long N, long n, int r, int d, int k)
{
   long Seq;                      /* Replication number */
   long NbGroups;                 /* Number of classes */
   long jhigh;
   long jlow;
   long Groupe;
   long L;
   int Minkd;
   int s, j;
   double X2;
   double Mult;
   double *NbExp;
   long *Loca;
   long *Nb;
   lebool Occurs[1 + Maxkd];
   double **M;
   double V[1];                   /* Number degrees of freedom for Chi2 */
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sknuth_SimpPoker test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataPoker (gen, TestName, N, n, r, d, k);

   util_Assert (d <= Maxkd, "sknuth_SimpPoker:   d > 127");
   util_Assert (k <= Maxkd, "sknuth_SimpPoker:   k > 127");
   util_Assert (d > 1, "sknuth_SimpPoker:   d < 2");
   util_Assert (k > 1, "sknuth_SimpPoker:   k < 2");
   if (k < d)
      Minkd = k;
   else
      Minkd = d;

   num2_CalcMatStirling (&M, Minkd, k);

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, Minkd, "sknuth_SimpPoker");
   NbExp = res->NbExp;
   Nb = res->Count;
   Loca = res->Loc;

   /* NbExp[s] = n * d * (d-1) * ... * (d-s+1) * M [s,k] / d^k.  */
   Mult = n * pow ((double) d, -(double) k);
   for (s = 1; s <= Minkd; s++) {
      Mult *= d - s + 1;
      NbExp[s] = Mult * M[s][k];
   }
   jlow = 1;
   jhigh = Minkd;
   if (swrite_Classes)
      gofs_WriteClasses (NbExp, Loca, jlow, jhigh, 0);
   gofs_MergeClasses (NbExp, Loca, &jlow, &jhigh, &NbGroups);
   if (swrite_Classes)
      gofs_WriteClasses (NbExp, Loca, jlow, jhigh, NbGroups);
   res->jmin = jlow;
   res->jmax = jhigh;
   res->degFree = NbGroups - 1;
   if (res->degFree < 1) {
      if (localRes)
         sres_DeleteChi2 (res);
      return;
   }
   sprintf (str, "The N statistic values (a ChiSquare with %1ld degrees"
                 " of freedom):", NbGroups - 1);
   statcoll_SetDesc (res->sVal1, str);

   for (Seq = 1; Seq <= N; Seq++) {
      for (s = 1; s <= Minkd; s++)
         Nb[s] = 0;
      for (Groupe = 1; Groupe <= n; Groupe++) {
         /* Draw one poker hand */
         for (j = 0; j < d; j++)
            Occurs[j] = FALSE;
         s = 0;                   /* s = number of different values */
         for (j = 1; j <= k; j++) {
            L = unif01_StripL (gen, r, d);
            if (!Occurs[L]) {
               Occurs[L] = TRUE;
               ++s;
            }
         }
         ++Nb[Loca[s]];
      }
      if (swrite_Counters)
         tables_WriteTabL (Nb, jlow, jhigh, 5, 10, "Observed numbers:");

      X2 = gofs_Chi2 (NbExp, Nb, jlow, jhigh);
      statcoll_AddObs (res->sVal1, X2);
   }

   V[0] = NbGroups - 1;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetChi2SumStat (res);

   if (swrite_Collectors) {
      statcoll_Write (res->sVal1, 5, 14, 4, 3);
   }
   if (swrite_Basic) {
      swrite_AddStrChi (str, LENGTH, res->degFree);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTest (N, res);
      swrite_Final (gen, Timer);
   }
   num2_FreeMatStirling (&M, Minkd);
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

#define MAXT 62

static void WriteDataCoupCol (unif01_Gen *gen, char *TestName,
   long N, long n, int r, int d)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   d = %4d\n\n", d);
}


/*-------------------------------------------------------------------------*/

static long NRepet (
   unif01_Gen *gen, 
   int dInt,                  /* d */
   int r,
   lebool Occurs[]
   )
/*
 * Used by CouponCollector. Counts the number of values generated before
 * each possible value of d appears at least once.
 */
{
   int u, j;
   int s = 0;

   for (j = 1; j <= dInt; j++) {
      do {
         ++s;
         if (s >= MAXT)
            return MAXT;
         u = unif01_StripL (gen, r, dInt);
      } while (Occurs[u]);
      Occurs[u] = TRUE;
   }
   /* j is the number of different values observed up to now */
   return s;
}


/*-------------------------------------------------------------------------*/

void sknuth_CouponCollector (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r, int d)
{
   long Seq;                      /* Replication number */
   long Segm;
   const int t = MAXT;
   long tt = t;
   int dInt = d;
   long dd = d;
   int s, k;
   long NbGroups;
   double Moydes;
   double Mult;
   double dReal = d;
   double **M;
   double *NbExp;
   long *Loca;
   long *Nb;
   lebool Occurs[1 + MAXT];
   double X2;
   double V[1];                   /* Number degrees of freedom for Chi2 */
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sknuth_CouponCollector test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataCoupCol (gen, TestName, N, n, r, d);

   util_Assert (d < MAXT, "sknuth_CouponCollector:  d >= 62");
   util_Assert (d > 1, "sknuth_CouponCollector:  d < 2");

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, MAXT, "sknuth_CouponCollector");
   NbExp = res->NbExp;
   Nb = res->Count;
   Loca = res->Loc;

   /* Compute the expected number of segments of each length */
   /* NbExp [s] = n * d! * Stirling (d-1, s-1) / d^s for d <= s <= t - 1 */
   /* NbExp [t] = n * (1 - d! * Stirling (d, t-1) / d^{t-1}) */
   dInt = d;
   num2_CalcMatStirling (&M, d, t - 1);
   Mult = n;
   for (s = 1; s <= d; s++) {
      Mult *= s / dReal;
   }
   NbExp[d] = Mult;
   Moydes = d * Mult;
   for (s = d + 1; s < t; s++) {
      Mult /= dReal;
      NbExp[s] = Mult * M[d - 1][s - 1];
      Moydes += s * NbExp[s];
   }
   NbExp[t] = n - Mult * M[d][t - 1];
   Moydes += t * NbExp[t];
   Moydes /= n;
 /* 
   if (swrite_Basic) {
       printf ("   Expected value of s = ");
	   num_WriteD (Moydes, 10, 2, 2);
      printf ("\n\n");
   }
 */
   if (swrite_Classes)
      gofs_WriteClasses (NbExp, Loca, d, t, 0);
   gofs_MergeClasses (NbExp, Loca, &dd, &tt, &NbGroups);
   if (swrite_Classes)
      gofs_WriteClasses (NbExp, Loca, dd, tt, NbGroups);
   res->jmin = dd;
   res->jmax = tt;
   res->degFree = NbGroups - 1;
   if (res->degFree < 1) {
      if (localRes)
         sres_DeleteChi2 (res);
      return;
   }

   sprintf (str, "The N statistic values (a ChiSquare with %1ld degrees"
                 " of freedom):", NbGroups - 1);
   statcoll_SetDesc (res->sVal1, str);

   /* Beginning of test */
   for (Seq = 1; Seq <= N; Seq++) {
      for (s = dInt; s <= MAXT; s++)
         Nb[s] = 0;
      for (Segm = 1; Segm <= n; Segm++) {
         /* One collection of values. */
         for (k = 0; k < dInt; k++)
            Occurs[k] = FALSE;
         ++Nb[Loca[NRepet (gen, dInt, r, Occurs)]];
      }
      if (swrite_Counters)
         tables_WriteTabL (Nb, dd, tt, 5, 10, "Observed numbers:");

      X2 = gofs_Chi2 (NbExp, Nb, dd, tt);
      statcoll_AddObs (res->sVal1, X2);
   }

   V[0] = NbGroups - 1;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetChi2SumStat (res);

   if (swrite_Collectors) {
      statcoll_Write (res->sVal1, 5, 14, 4, 3);
   }
   if (swrite_Basic) {
      swrite_AddStrChi (str, LENGTH, res->degFree);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTest (N, res);
      swrite_Final (gen, Timer);
   }
   num2_FreeMatStirling (&M, d);
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void sknuth_Permutation (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r, int t)
{
   double ValDelta[] = { 1.0 };
   smultin_Param *par;

   if (swrite_Basic)
      printf ("***********************************************************\n"
              "Test sknuth_Permutation calling smultin_Multinomial\n\n");

   par = smultin_CreateParam (1, ValDelta, smultin_GenerCellPermut, 3);
   if (NULL == res) {
      smultin_Multinomial (gen, par, NULL, N, n, r, 1, t, FALSE);
   } else {
      smultin_Res *resm;
      resm = smultin_CreateRes (par);
      smultin_Multinomial (gen, par, resm, N, n, r, 1, t, FALSE);
      sres_InitChi2 (res, N, -1, "sknuth_Permutation");
      statcoll_SetDesc (res->sVal1, "Serial sVal1");
      res->sVal1->NObs = resm->Collector[0]->NObs;
      tables_CopyTabD (resm->Collector[0]->V, res->sVal1->V, 1, N);
      tables_CopyTabD (resm->sVal2[0], res->sVal2, 0, gofw_NTestTypes - 1);
      tables_CopyTabD (resm->pVal2[0], res->pVal2, 0, gofw_NTestTypes - 1);
      smultin_DeleteRes (resm);
   }
   smultin_DeleteParam (par);
}


/*=========================================================================*/

static void WriteDataRun (unif01_Gen * gen, char *TestName,
   long N, long n, int r, lebool Up)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   Up = %5s\n\n", Up ? "TRUE" : "FALSE");
}


/*-------------------------------------------------------------------------*/

void sknuth_Run (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r, lebool Up)
{
   long Seq;                      /* Replication number */
   double U;
   double UPrec;                  /* Preceding value of U */
   double nReal = n;
   double A[6][6];
   double B[6];
   double *NbExp;
   long k;
   int j, i;
   long Longueur;                 /* Current length of the sequence */
   double Khi;
   long *Count;
   char str[LENGTH + 1];
   double V[1];                   /* Number degrees of freedom for Chi2 */
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sknuth_Run test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataRun (gen, TestName, N, n, r, Up);

   if (n < 600)
      return;
   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, 6, "sknuth_Run");
   NbExp = res->NbExp;
   Count = res->Count;
   res->jmin = 1;
   res->jmax = 6;

   A[0][0] =   4529.35365;
   A[0][1] =   9044.90208;
   A[0][2] =  13567.9452;
   A[0][3] =  18091.2672;
   A[0][4] =  22614.7139;
   A[0][5] =  27892.1588;
   A[1][1] =  18097.0254;
   A[1][2] =  27139.4552;
   A[1][3] =  36186.6493;
   A[1][4] =  45233.8198;
   A[1][5] =  55788.8311;
   A[2][2] =  40721.3320;
   A[2][3] =  54281.2656;
   A[2][4] =  67852.0446;
   A[2][5] =  83684.5705;
   A[3][3] =  72413.6082;
   A[3][4] =  90470.0789;
   A[3][5] = 111580.110;
   A[4][4] = 113261.815;
   A[4][5] = 139475.555;
   A[5][5] = 172860.170;

   for (i = 2; i <= 6; i++) {
      for (j = 1; j < i; j++)
         A[i - 1][j - 1] = A[j - 1][i - 1];
   }

   B[0] = 1.0 / 6.0;
   B[1] = 5.0 / 24.0;
   B[2] = 11.0 / 120.0;
   B[3] = 19.0 / 720.0;
   B[4] = 29.0 / 5040.0;
   B[5] = 1.0 / 840.0;
   for (i = 1; i <= 6; i++) {
      NbExp[i] = nReal * B[i - 1];
      res->Loc[i] = i;
   }

   if (swrite_Classes)
      /* gofs_Classes (NbExp, NULL, 1, 6, 0); */
      tables_WriteTabD (NbExp, 1, 6, 1, 20, 2, 1, "Expected numbers:");

   statcoll_SetDesc (res->sVal1,
      "The N statistic values (a ChiSquare with 6 degrees of freedom):");
   res->degFree = 6;

   /* Beginning of test */
   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 1; i <= 6; i++)
         Count[i] = 0;
      Longueur = 1;
      UPrec = unif01_StripD (gen, r);
      /* Generate n numbers */
      for (k = 1; k < n; k++) {
         U = unif01_StripD (gen, r);
         if ((Up && U < UPrec) || (!Up && U > UPrec)) {
            /* The end of a "Run" */
            ++Count[Longueur];
            Longueur = 1;
         } else if (Longueur < 6)
            ++Longueur;
         UPrec = U;
      }
      ++Count[Longueur];

      if (swrite_Counters)
         tables_WriteTabL (Count, 1, 6, 5, 10, "Observed numbers:");

      /* Compute modified Chi2 for a sequence */
      Khi = 0.0;
      for (i = 1; i <= 6; i++) {
	 for (j = 1; j <= 6; j++) {
	    Khi += A[i-1][j-1]*(Count[i] - NbExp[i])*(Count[j] - NbExp[j]);
	 }
      }
      statcoll_AddObs (res->sVal1, Khi / (nReal - 6.0));
   }

   V[0] = 6;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetChi2SumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      swrite_AddStrChi (str, LENGTH, res->degFree);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTest (N, res);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void sknuth_RunIndep (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r, lebool Up)
{
   long Seq;                      /* Replication number */
   double U;
   double UPrec;                  /* Preceding value of U */
   double X2;
   long Nb;
   long k;
   int i;
   long Longueur;                 /* Current length of the sequence */
   long *Count;
   double *NbExp;
   double Prob[7];
   char str[LENGTH + 1];
   double V[1];                   /* Number degrees of freedom for Chi2 */
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sknuth_RunIndep test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataRun (gen, TestName, N, n, r, Up);

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, 6, "sknuth_RunIndep");
   NbExp = res->NbExp;
   Count = res->Count;
   res->jmin = 1;
   res->jmax = 6;
   sprintf (str, "NumExpected[6] < %.1f", gofs_MinExpected);

   for (i = 1; i <= 5; i++) {
      Prob[i] = 1.0 / num2_Factorial (i) - 1.0 / num2_Factorial (i + 1);
   }
   Prob[6] = 1.0 / num2_Factorial (6);

   statcoll_SetDesc (res->sVal1,
      "The N statistic values (a ChiSquare with 5 degrees of freedom):");
   res->degFree = 5;

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 1; i <= 6; i++)
         Count[i] = 0;
      Longueur = 1;
      UPrec = unif01_StripD (gen, r);
      for (k = 1; k <= n; k++) {
         U = unif01_StripD (gen, r);
         if ((Up && U < UPrec) || (!Up && U > UPrec)) {
            /* The end of a "Run" */
            ++Count[Longueur];
            Longueur = 1;
            U = unif01_StripD (gen, r);
         } else if (Longueur < 6)
            ++Longueur;
         UPrec = U;
      }
      ++Count[Longueur];

      Nb = 0;
      for (i = 1; i <= 6; i++)
         Nb += Count[i];
      for (i = 1; i <= 6; i++)
         NbExp[i] = Nb * Prob[i];

      if (swrite_Counters) {
         tables_WriteTabD (NbExp, 1, 6, 1, 20, 2, 1, "Expected numbers:");
         tables_WriteTabL (Count, 1, 6, 1, 17, "Observed numbers:");
      }
      /*     util_Warning (NbExp[6] < gofs_MinExpected, str); */

      X2 = gofs_Chi2 (NbExp, Count, 1, 6);
      statcoll_AddObs (res->sVal1, X2);
   }

   V[0] = 5;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetChi2SumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      swrite_AddStrChi (str, LENGTH, res->degFree);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTest (N, res);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataMaxOft (unif01_Gen * gen, char *TestName,
   long N, long n, int r, int d, int t, double NbExp)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   d = %4d,   t = %2d\n\n", d, t);
   printf ("      Number of categories = %d\n", d);
   printf ("      Expected number per category  = %.2f\n\n", NbExp);
}


/*-------------------------------------------------------------------------*/

static double FDistMax (
   double Par[],             /* The parameter t = Par[0] */
   double x                  /* The argument */
   )
/*
 * Distribution function for the maximum of t random variables U01 = x^t
 */
{
  /*   double Prod;
   int j;
   const int t = Par[0] + 0.5;
  */
   if (x >= 1.0)
      return 1.0;
   if (x <= 0.0)
      return 0.0;
   return pow (x, Par[0]);
   /*
   Prod = x;
   for (j = 1; j < t; j++)
      Prod *= x;
   return Prod;
   */
}


/*-------------------------------------------------------------------------*/

void sknuth_MaxOft (unif01_Gen * gen, sknuth_Res1 * res,
   long N, long n, int r, int d, int t)
{
   long Seq;                      /* Replication number */
   double tReal = t;
   double dReal = d;
   double NbExp;                  /* Expected number in each class */
   double MaxU;
   double U;
   long Groupe;
   int j, Indice;
   double *P;
   double Par[1];
   double X2;
   double V[1];                   /* Number degrees of freedom for Chi2 */
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sknuth_MaxOft test";
   sres_Basic *Bas;
   sres_Chi2 *Chi;

   Timer = chrono_Create ();
   Par[0] = t;

   NbExp = n / dReal;
   if (swrite_Basic)
      WriteDataMaxOft (gen, TestName, N, n, r, d, t, NbExp);
   util_Assert (NbExp >= gofs_MinExpected,
      "MaxOft:   NbExp < gofs_MinExpected");
   if (res == NULL) {
      localRes = TRUE;
      res = sknuth_CreateRes1 ();
   }
   InitRes1 (res, N, d);
   Bas = res->Bas;
   Chi = res->Chi;
   Chi->jmin = 0;
   Chi->jmax = d - 1;
   for (j = 0; j < d; j++) {
      Chi->Loc[j] = j;
      Chi->NbExp[j] = NbExp;
   }

   sprintf (str, "The N statistic values (a ChiSquare with %1d degrees"
                 " of freedom):", d - 1);
   statcoll_SetDesc (Chi->sVal1, str);
   Chi->degFree = d - 1;
   statcoll_SetDesc (Bas->sVal1,
      "The N statistic values (the Anderson-Darling p-values):");
   P = util_Calloc ((size_t) n + 1, sizeof (double));

   for (Seq = 1; Seq <= N; Seq++) {
      for (Indice = 0; Indice < d; Indice++)
         Chi->Count[Indice] = 0;
      for (Groupe = 1; Groupe <= n; Groupe++) {
         /* Generate a vector and find the max value */
         MaxU = unif01_StripD (gen, r);
         for (j = 1; j < t; j++) {
            U = unif01_StripD (gen, r);
            if (U > MaxU)
               MaxU = U;
         }
         /* For the chi2 */
         Indice = pow (MaxU, tReal) * dReal;
         ++Chi->Count[Indice];

         /* For the Anderson-Darling */
         P[Groupe] = MaxU;
      }
      if (swrite_Counters)
         tables_WriteTabL (Chi->Count, 0, d - 1, 5, 10, "Observed numbers:");

      /* Value of the chi2 statistic */
      X2 = gofs_Chi2Equal (NbExp, Chi->Count, 0, d - 1);
      statcoll_AddObs (Chi->sVal1, X2);

      /* Value of the Anderson-Darling statistic */
      gofw_ActiveTests1 (P, n, FDistMax, Par, Bas->sVal2, Bas->pVal2);
      statcoll_AddObs (Bas->sVal1, Bas->pVal2[gofw_AD]);
   }
   util_Free (P);

   V[0] = d - 1;
   gofw_ActiveTests2 (Chi->sVal1->V, Chi->pVal1->V, N, wdist_ChiSquare, V,
      Chi->sVal2, Chi->pVal2);
   Chi->pVal1->NObs = N;
   sres_GetChi2SumStat (Chi);

   gofw_ActiveTests2 (Bas->sVal1->V, Bas->pVal1->V, N, wdist_Unif,
      (double *) NULL, Bas->sVal2, Bas->pVal2);
   Bas->pVal1->NObs = N;

   if (swrite_Collectors) {
      statcoll_Write (Chi->sVal1, 5, 14, 4, 3);
      statcoll_Write (Bas->sVal1, 5, 14, 4, 3);
   }
   if (swrite_Basic) {
      if (N == 1) {
         swrite_AddStrChi (str, LENGTH, Chi->degFree);
         gofw_WriteActiveTests2 (N, Chi->sVal2, Chi->pVal2, str);
      } else {
         printf ("\n-----------------------------------------------\n");
         printf ("Test results for chi2 with %2ld degrees of freedom:\n",
                 Chi->degFree);
         gofw_WriteActiveTests0 (N, Chi->sVal2, Chi->pVal2);
         swrite_Chi2SumTest (N, Chi);
      }

      if (N == 1) {
         gofw_WriteActiveTests2 (N, Bas->sVal2, Bas->pVal2,
            "Anderson-Darling statistic            :");
      } else {
         printf ("\n-----------------------------------------------\n");
         printf ("Test results for Anderson-Darling:\n");
         gofw_WriteActiveTests0 (N, Bas->sVal2, Bas->pVal2);
      }
      printf ("\n");
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sknuth_DeleteRes1 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void sknuth_Collision (unif01_Gen * gen, sknuth_Res2 * res,
   long N, long n, int r, long d, int t)
{
   double ValDelta[] = { -1.0 };
   smultin_Param *par;

   if (swrite_Basic)
      printf ("***********************************************************\n"
              "Test sknuth_Collision calling smultin_Multinomial\n\n");

   par = smultin_CreateParam (1, ValDelta, smultin_GenerCellSerial, -3);
   if (NULL == res) {
      smultin_Multinomial (gen, par, NULL, N, n, r, d, t, TRUE);
   } else {
      smultin_Res *resm;
      resm = smultin_CreateRes (par);
      smultin_Multinomial (gen, par, resm, N, n, r, d, t, TRUE);
      InitRes2 (res, N, resm->Mu[0], "sknuth_Collision");
      statcoll_SetDesc (res->Bas->sVal1, "Collision sVal1");
      statcoll_SetDesc (res->Pois->sVal1, "Collision sVal1");
      res->Pois->sVal1->NObs = resm->Collector[0]->NObs;
      res->Bas->sVal1->NObs = resm->Collector[0]->NObs;
      res->Pois->pLeft = resm->pCollLeft;
      res->Pois->pRight = resm->pCollRight;
      tables_CopyTabD (resm->Collector[0]->V, res->Bas->sVal1->V, 1, N);
      tables_CopyTabD (resm->Collector[0]->V, res->Pois->sVal1->V, 1, N);
      res->Pois->pVal2 = resm->pColl;
      res->Pois->sVal2 = resm->NbCollisions;
      tables_CopyTabD (resm->sVal2[0], res->Bas->sVal2, 0,
         gofw_NTestTypes - 1);
      tables_CopyTabD (resm->pVal2[0], res->Bas->pVal2, 0,
         gofw_NTestTypes - 1);
      smultin_DeleteRes (resm);
   }
   smultin_DeleteParam (par);
}


/*=========================================================================*/

void sknuth_CollisionPermut (unif01_Gen * gen, sknuth_Res2 * res,
   long N, long n, int r, int t)
{
   double ValDelta[] = { -1.0 };
   smultin_Param *par;

   if (swrite_Basic)
      printf ("***********************************************************\n"
         "Test sknuth_CollisionPermut calling smultin_Multinomial\n\n");

   par = smultin_CreateParam (1, ValDelta, smultin_GenerCellPermut, -3);
   if (NULL == res) {
      smultin_Multinomial (gen, par, NULL, N, n, r, 0, t, TRUE);
   } else {
      smultin_Res *resm;
      resm = smultin_CreateRes (par);
      smultin_Multinomial (gen, par, resm, N, n, r, 0, t, TRUE);
      InitRes2 (res, N, resm->Mu[0], "sknuth_CollisionPermut");
      statcoll_SetDesc (res->Bas->sVal1, "CollisionPermut sVal1");
      statcoll_SetDesc (res->Pois->sVal1, "CollisionPermut sVal1");
      res->Pois->pLeft = resm->pCollLeft;
      res->Pois->pRight = resm->pCollRight;
      res->Pois->sVal1->NObs = resm->Collector[0]->NObs;
      res->Bas->sVal1->NObs = resm->Collector[0]->NObs;
      tables_CopyTabD (resm->Collector[0]->V, res->Bas->sVal1->V, 1, N);
      tables_CopyTabD (resm->Collector[0]->V, res->Pois->sVal1->V, 1, N);
      res->Pois->pVal2 = resm->pColl;
      res->Pois->sVal2 = resm->NbCollisions;
      tables_CopyTabD (resm->sVal2[0], res->Bas->sVal2, 0,
         gofw_NTestTypes - 1);
      tables_CopyTabD (resm->pVal2[0], res->Bas->pVal2, 0,
         gofw_NTestTypes - 1);
      smultin_DeleteRes (resm);
   }
   smultin_DeleteParam (par);
}
