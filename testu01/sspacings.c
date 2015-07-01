/*************************************************************************\
 *
 * Package:        TestU01
 * File:           sspacings.c
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
#include "chrono.h"
#include "num.h"
#include "tables.h"

#include "sspacings.h"
#include "unif01.h"
#include "wdist.h"
#include "swrite.h"

#include "statcoll.h"
#include "gofw.h"
#include "fbar.h"

#include <math.h>
#include <string.h>
#include <stdio.h>






/*------------------------------- Constants -------------------------------*/

/* String length */
#define LEN 50

/* Euler's constant */
static const double EULER = 0.577215664901533;

/* Necessary so that the multiplication of a large number of U01 will not
   underflow; the precise value is not so important: it could be bigger. */
static const double Epsilon = 1.0E-200;


#define MAXM 201


/*--------------------------------- Types -------------------------------*/

typedef enum
{
   LOG_EXACT_CIRC,              /* Log of spacings, exact, circular */
   LOG_EXACT_LIN,               /* Log of spacings, exact, linear */
   LOG_ASYMP_CIRC,              /* Log of spacings, asymptotic, circular */
   LOG_ASYMP_LIN,               /* Log of spacings, asymptotic, linear */
   SQUA_EXACT_CIRC,             /* Square of spacings, exact, circular */
   SQUA_EXACT_LIN,              /* Square of spacings, exact, linear */
   SQUA_ASYMP_CIRC,             /* Square of spacings, asymptotic, circular */
   SQUA_ASYMP_LIN,              /* Square of spacings, asymptotic, linear */
   Stats_N                      /* Number of statistics in this enum */
} Stats;


typedef struct
{
   int NbColl;                    /* Number of collectors used for each m */
   int Nbm;                       /* Number of values of m to consider */
   int Loc[MAXM];                 /* The statistic for m is in Loc [m] */
   double Mu[MAXM][Stats_N];      /* Mean */
   double Sig[MAXM][Stats_N];     /* Standard deviation */
   double HMu[MAXM][Stats_N];     /* Empirical mean */
   double HSig[MAXM][Stats_N];    /* Empirical standard deviation */
} Param;



/*------------------------------- Functions -------------------------------*/


sspacings_Res * sspacings_CreateRes (void)
{
   sspacings_Res *res;
   res = util_Malloc (sizeof (sspacings_Res));
   memset (res, 0, sizeof (sspacings_Res));
   res->name = util_Calloc (1, sizeof (char));
   res->smax = -1;
   return res;
}


/*-------------------------------------------------------------------------*/

void sspacings_DeleteRes (sspacings_Res *res)
{
   int j;
   if (res == NULL)
      return;
   for (j = 0; j <= res->smax; j += 2)
      res->Collectors[j] = statcoll_Delete (res->Collectors[j]);
   util_Free (res->Collectors);

   for (j = 0; j <= res->imax; j++) {
      sres_DeleteBasic (res->LogCAMu[j]);
      sres_DeleteBasic (res->LogCEMu[j]);
      sres_DeleteBasic (res->SquareCAMu[j]);
      sres_DeleteBasic (res->SquareCEMu[j]);
   }

   util_Free (res->LogCEMu);
   util_Free (res->LogCAMu);
   util_Free (res->SquareCEMu);
   util_Free (res->SquareCAMu);

   util_Free (res->LogCESig_sVal);
   util_Free (res->LogCESig_pVal);
   util_Free (res->LogCASig_sVal);
   util_Free (res->LogCASig_pVal);
   util_Free (res->SquareCESig_sVal);
   util_Free (res->SquareCESig_pVal);
   util_Free (res->SquareCASig_sVal);
   util_Free (res->SquareCASig_pVal);

   util_Free (res->name);
   util_Free (res);
}


/*-------------------------------------------------------------------------*/

static void InitRes (
   sspacings_Res *res,          /* Results holder */
   long N,                    /* Number of replications */
   int Nbm,                   /* Number of values of m to consider */
   char *nam                  /* Test name */
)
/* 
 * Initializes res
 */
{
   char nom[LEN + 1];
   char spindex[LEN + 1];
   int j;

   if (res->smax < 0) {
      res->Collectors = util_Calloc ((size_t) 8 * Nbm,
         sizeof (statcoll_Collector *));
      for (j = 0; j < 8 * Nbm; j += 2)
         res->Collectors[j] = statcoll_Create (N, "");

      res->LogCAMu = util_Calloc ((size_t) Nbm, sizeof (sres_Basic *));
      res->LogCEMu = util_Calloc ((size_t) Nbm, sizeof (sres_Basic *));
      res->SquareCAMu = util_Calloc ((size_t) Nbm, sizeof (sres_Basic *));
      res->SquareCEMu = util_Calloc ((size_t) Nbm, sizeof (sres_Basic *));
      for (j = 0; j < Nbm; j++) {
         res->LogCAMu[j] = sres_CreateBasic ();
         res->LogCEMu[j] = sres_CreateBasic ();
         res->SquareCAMu[j] = sres_CreateBasic ();
         res->SquareCEMu[j] = sres_CreateBasic ();
      }

      res->LogCESig_sVal = util_Calloc ((size_t) Nbm, sizeof (double));
      res->LogCESig_pVal = util_Calloc ((size_t) Nbm, sizeof (double));
      res->LogCASig_sVal = util_Calloc ((size_t) Nbm, sizeof (double));
      res->LogCASig_pVal = util_Calloc ((size_t) Nbm, sizeof (double));
      res->SquareCESig_sVal = util_Calloc ((size_t) Nbm, sizeof (double));
      res->SquareCESig_pVal = util_Calloc ((size_t) Nbm, sizeof (double));
      res->SquareCASig_sVal = util_Calloc ((size_t) Nbm, sizeof (double));
      res->SquareCASig_pVal = util_Calloc ((size_t) Nbm, sizeof (double));

   } else {
      for (j = 8 * Nbm; j <= res->smax; j += 2)
         res->Collectors[j] = statcoll_Delete (res->Collectors[j]);
      res->Collectors = util_Realloc (res->Collectors,
         8 * Nbm * sizeof (statcoll_Collector *));
      for (j = res->smax + 2; j < 8 * Nbm; j += 2)
         res->Collectors[j] = statcoll_Create (N, "");

      for (j = Nbm; j <= res->imax; j++) {
         sres_DeleteBasic (res->LogCAMu[j]);
         sres_DeleteBasic (res->LogCEMu[j]);
         sres_DeleteBasic (res->SquareCAMu[j]);
         sres_DeleteBasic (res->SquareCEMu[j]);
      }
      res->LogCAMu = util_Realloc (res->LogCAMu,
         Nbm * sizeof (sres_Basic *));
      res->LogCEMu = util_Realloc (res->LogCEMu,
         Nbm * sizeof (sres_Basic *));
      res->SquareCAMu = util_Realloc (res->SquareCAMu,
         Nbm * sizeof (sres_Basic *));
      res->SquareCEMu = util_Realloc (res->SquareCEMu,
         Nbm * sizeof (sres_Basic *));

      for (j = res->imax + 1; j < Nbm; j++) {
         res->LogCAMu[j] = sres_CreateBasic ();
         res->LogCEMu[j] = sres_CreateBasic ();
         res->SquareCAMu[j] = sres_CreateBasic ();
         res->SquareCEMu[j] = sres_CreateBasic ();
      }

      res->LogCESig_sVal =
         util_Realloc (res->LogCESig_sVal, Nbm * sizeof (double));
      res->LogCESig_pVal =
         util_Realloc (res->LogCESig_pVal, Nbm * sizeof (double));
      res->LogCASig_sVal =
         util_Realloc (res->LogCASig_sVal, Nbm * sizeof (double));
      res->LogCASig_pVal =
         util_Realloc (res->LogCASig_pVal, Nbm * sizeof (double));
      res->SquareCESig_sVal =
         util_Realloc (res->SquareCESig_sVal, Nbm * sizeof (double));
      res->SquareCESig_pVal =
         util_Realloc (res->SquareCESig_pVal, Nbm * sizeof (double));
      res->SquareCASig_sVal =
         util_Realloc (res->SquareCASig_sVal, Nbm * sizeof (double));
      res->SquareCASig_pVal =
         util_Realloc (res->SquareCASig_pVal, Nbm * sizeof (double));
   }

   for (j = 0; j < 8 * Nbm; j += 2)
      statcoll_Init (res->Collectors[j], N);
   res->smax = 8 * Nbm - 2;

   for (j = 0; j < Nbm; j++) {
      strncpy (nom, "LogCEMu[", (size_t) LEN);
      sprintf (spindex, "%1d", j);
      strncat (nom, spindex, (size_t) 5);
      strncat (nom, "]", (size_t) 2);
      sres_InitBasic (res->LogCEMu[j], N, nom);

      strncpy (nom, "LogCAMu[", (size_t) LEN);
      sprintf (spindex, "%1d", j);
      strncat (nom, spindex, (size_t) 5);
      strncat (nom, "]", (size_t) 2);
      sres_InitBasic (res->LogCAMu[j], N, nom);

      strncpy (nom, "SquareCEMu[", (size_t) LEN);
      sprintf (spindex, "%1d", j);
      strncat (nom, spindex, (size_t) 5);
      strncat (nom, "]", (size_t) 2);
      sres_InitBasic (res->SquareCEMu[j], N, nom);

      strncpy (nom, "SquareCAMu[", (size_t) LEN);
      sprintf (spindex, "%1d", j);
      strncat (nom, spindex, (size_t) 5);
      strncat (nom, "]", (size_t) 2);
      sres_InitBasic (res->SquareCAMu[j], N, nom);

      res->imax = j;
   }
   res->name = util_Realloc (res->name, 1 + strlen (nam) * sizeof (char));
   strcpy (res->name, nam);

   memset (res->LogCESig_sVal, 0, sizeof (res->LogCESig_sVal));
   memset (res->LogCESig_pVal, 0, sizeof (res->LogCESig_pVal));
   memset (res->LogCASig_sVal, 0, sizeof (res->LogCASig_sVal));
   memset (res->LogCASig_pVal, 0, sizeof (res->LogCASig_pVal));
   memset (res->SquareCESig_sVal, 0, sizeof (res->SquareCESig_sVal));
   memset (res->SquareCESig_pVal, 0, sizeof (res->SquareCESig_pVal));
   memset (res->SquareCASig_sVal, 0, sizeof (res->SquareCASig_sVal));
   memset (res->SquareCASig_pVal, 0, sizeof (res->SquareCASig_pVal));
}


/*=========================================================================*/


void sspacings_SumLogsSpacings (unif01_Gen * gen, sspacings_Res * res,
   long N, long n, int r, int m)
{
   util_Error ("sspacings_SumLogsSpacings is implemented inside "
               "sspacings_AllSpacings2");
}


/*=========================================================================*/

void sspacings_SumSquaresSpacings (unif01_Gen * gen, sspacings_Res * res,
   long N, long n, int r, int m)
{
   util_Error ("sspacings_SumSquaresSpacings is implemented inside "
               "sspacings_AllSpacings2");
}


/*=========================================================================*/

void sspacings_ScanSpacings (unif01_Gen * gen, sspacings_Res * res,
   long N, long n, int r, double d)
{
   util_Error ("sspacings_ScanSpacings not implemented");
}


/*=========================================================================*/

static void WrRes (char S1[], long N, Param * par, int m, int i,
   statcoll_Collector * Collectors[], gofw_TestArray sv, gofw_TestArray pv)
{
   printf ("%s", S1);
   printf ("\n   Mu    = ");
   num_WriteD (par->Mu[m][i], 12, 8, 7);
   printf ("\n   Sigma = ");
   num_WriteD (par->Sig[m][i], 12, 8, 7);
   printf ("\n\nEmpirical mean of standardized values :");
   num_WriteD (par->HMu[m][i] / N, 12, 8, 7);
   printf ("\n");
   gofw_Writep1 (fbar_Normal1 (par->HMu[m][i] / N));
   printf ("Second empirical moment of standardized values:");
   num_WriteD (par->HSig[m][i] / N, 12, 8, 7);
   printf ("\n");
   gofw_Writep1 (fbar_ChiSquare2 (N, 12, par->HSig[m][i]));
   i += 8 * par->Loc[m];
   if (N > 1)
      gofw_WriteActiveTests0 (N, sv, pv);

   if (swrite_Collectors) {
      statcoll_Write (Collectors[i], 5, 14, 4, 3);
      printf ("\n");
   }
   printf ("\n");
}

/*-------------------------------------------------------------------------*/

static void UpdateStat (Param * par, int m, int i, double v,
   statcoll_Collector * Collectors[])
{
   double x;

   x = (v - par->Mu[m][i]) / par->Sig[m][i];
   par->HMu[m][i] += x;
   par->HSig[m][i] += x * x;
   i += 8 * par->Loc[m];
   statcoll_AddObs (Collectors[i], x);
}

/*-------------------------------------------------------------------------*/

static void InitAllSpacings (unif01_Gen * gen, char *TestName, Param * par,
   long N, long n0, int r, int M0, int M1, int D, int LgEps)
{
   double Rmn[MAXM];              /* R(m,n) and Q(m,n) from Cressie (1976) */
   double Qmn[MAXM];
   double Rmm1[MAXM];
   double Qmm1[MAXM];
   double Rm1;
   double Qm1;                    /* R(1,m-1) and Q(1,m-1) from Cressie */
   double Fact;                   /* n^2 / (n+2) */
   double Mu0;
   double x, y;
   double I;
   double R;
   double Q;
   double nLR;
   double mLR;
   int m;
   int i;

   if (swrite_Basic) {
      swrite_Head (gen, TestName, N, n0, r);
      printf (",   M0 = %1d,   M1 = %1d,   D  = %1d\n", M0, M1, D);
      printf ("   LgEps = %1d\n\n\n", LgEps);
   }
   util_Assert (M1 < MAXM, "InitAllSpacings:   M1 is too large");

   par->Nbm = 1 + (M1 - M0) / D;
   for (i = 0; i < par->Nbm; i++)
      par->Loc[M0 + i * D] = i;
   if (M0 == 0)
      par->Loc[1] = 0;

   /* Initialize constants and arrays */
   nLR = n0;
   Fact = nLR * nLR / (nLR + 2.0);
   R = 0.0;
   Q = 0.0;

   for (i = n0; i >= M1; i--) {
      I = 1.0 / i;
      R += I;
      Q += I * I;
   }

   Rmn[M1] = R;
   Qmn[M1] = Q;
   for (m = M1 - 1; m >= 1; m--) {
      I = 1.0 / m;
      Rmn[m] = Rmn[m + 1] + I;
      Qmn[m] = Qmn[m + 1] + I * I;
   }

   m = 1;
   Rmm1[1] = 0.0;
   Qmm1[1] = 0.0;
   for (m = 2; m <= M1; m++) {
      I = 1.0 / (m - 1);
      Rmm1[m] = Rmm1[m - 1] + I;
      Qmm1[m] = Qmm1[m - 1] + I * I;
   }

   /* Compute theoretical means and variances */
   if (M0 == 0)
      m = 1;
   else
      m = M0;

   while (m <= M1) {
      mLR = m;
      Rm1 = Rmm1[m];
      Qm1 = Qmm1[m];
      par->Mu[m][LOG_EXACT_CIRC] = -(nLR + 1.0) * Rmn[m];
      par->Mu[m][LOG_EXACT_LIN] =
         par->Mu[m][LOG_EXACT_CIRC] * (nLR + 2.0 - mLR) / (nLR + 1.0);
      par->Mu[m][LOG_ASYMP_CIRC] =
         -(nLR + 1.0) * (log (nLR + 1.0) + EULER - Rm1);
      par->Mu[m][LOG_ASYMP_LIN] =
         par->Mu[m][LOG_ASYMP_CIRC] * (nLR + 2.0 - mLR) / (nLR + 1.0);
      Mu0 = mLR * (mLR + 1.0);
      par->Mu[m][SQUA_EXACT_CIRC] = Mu0 * Fact;
      par->Mu[m][SQUA_EXACT_LIN] =
         par->Mu[m][SQUA_EXACT_CIRC] * (nLR - mLR + 2.0) / (nLR + 1.0);
      par->Mu[m][SQUA_ASYMP_CIRC] = Mu0 * (nLR + 1.0);
      par->Mu[m][SQUA_ASYMP_LIN] = Mu0 * (nLR - mLR + 2.0);

      x = (2 * m * (m - 1) + 1) * ((num_Pi * num_Pi) / 6.0 - Qm1) +
         (-2 * m + 1);
      util_Assert (x > 0.0, "Negative Sig [m, 2]");
      par->Sig[m][LOG_ASYMP_CIRC] = sqrt (nLR * x);
      par->Sig[m][LOG_ASYMP_LIN] = par->Sig[m][LOG_ASYMP_CIRC]; /* See Holst 
                                                                   1979 */

      x = Qmn[m] + nLR * Qmn[1] - 2.0 * (mLR - 1.0) * (mLR * Qm1 + 1.0) +
         (2.0 * mLR * (mLR - 1.0) - nLR) * num_Pi * num_Pi / 6.0;
      util_Assert (x > 0.0, "Negative Sig [m, 0] ...");
      par->Sig[m][LOG_EXACT_CIRC] = sqrt ((nLR + 1.0) * x);

      y = 2 * m * (m + 1) * (2 * m + 1) / 3.0;
      par->Sig[m][SQUA_ASYMP_CIRC] = sqrt (nLR * y);
      par->Sig[m][SQUA_ASYMP_LIN] = par->Sig[m][SQUA_ASYMP_CIRC];

      x = 2.0 * mLR * (1.0 + mLR) * (2.0 + mLR * (1.0 - 3.0 * mLR) +
         nLR * (1.0 + 2.0 * mLR)) / 3.0;
      x = x / ((nLR + 3.0) * (nLR + 4.0));
      util_Assert (x > 0.0, "Negative Sig [m, 4]");
      par->Sig[m][SQUA_EXACT_CIRC] = sqrt (x) * Fact;

      x = 20.0 + mLR * (-54.0 + mLR * (6.0 + mLR * (58.0 - 30.0 * mLR)));
      x += nLR * (34.0 + mLR * (-37.0 + mLR * (-27.0 + mLR * (48.0 -
                  12.0 * mLR))));
      x += nLR * nLR * (16.0 + mLR * (3.0 + mLR * (-15.0 + 8.0 * mLR)));
      x += nLR * nLR * nLR * 2.0 * (1.0 + 2.0 * mLR);
      x = x * mLR * (1.0 + mLR) / 3.0;
      x = x / ((nLR + 3.0) * (nLR + 4.0));
      util_Assert (x > 0.0, "Negative Sig [m, 5]");
      par->Sig[m][SQUA_EXACT_LIN] = sqrt (x) * Fact / (nLR + 1.0);

      /* Initalization of tables for the tests */
      for (i = 0; i < Stats_N; i++) {
         par->HMu[m][i] = 0.0;
         par->HSig[m][i] = 0.0;
      }
      if (M0 == 0 && m == 1)
         m = D;
      else
         m += D;
   }
   /*
   if (swrite_Basic)
      printf ("   Number of collectors initialized to N:   %1d\n\n",
         par->NbColl * par->Nbm);
   */
}


/*=========================================================================*/

static void CopyResults (sspacings_Res * res, Param * par, long N,
   int M0, int M1, int D, int flag)
{
   int m, s, i, k;

   if (M0 == 0)
      m = 1;
   else
      m = M0;
   s = 0;

   while (m <= M1) {
      i = LOG_EXACT_CIRC;
      k = i + 8 * par->Loc[m];
      tables_CopyTabD (res->Collectors[k]->V, res->LogCEMu[s]->sVal1->V, 1,
         N);
      res->LogCEMu[s]->sVal1->NObs = N;
      gofw_ActiveTests2 (res->Collectors[k]->V, res->LogCEMu[s]->pVal1->V,
         N, wdist_Normal, (double *) NULL, res->LogCEMu[s]->sVal2,
         res->LogCEMu[s]->pVal2);

      res->LogCEMu[s]->sVal2[gofw_Mean] = par->HMu[m][i] / N;
      res->LogCEMu[s]->pVal2[gofw_Mean] = fbar_Normal1 (par->HMu[m][i] / N);
      res->LogCESig_sVal[s] = par->HSig[m][i] / N;
      res->LogCESig_pVal[s] = fbar_ChiSquare2 (N, 12, par->HSig[m][i]);

      i = LOG_ASYMP_CIRC;
      k = i + 8 * par->Loc[m];
      tables_CopyTabD (res->Collectors[k]->V, res->LogCAMu[s]->sVal1->V, 1,
         N);
      res->LogCAMu[s]->sVal1->NObs = N;
      if (flag) {
         gofw_ActiveTests2 (res->Collectors[k]->V, res->LogCAMu[s]->pVal1->V,
            N, wdist_Normal, (double *) NULL, res->LogCAMu[s]->sVal2,
            res->LogCAMu[s]->pVal2);

         res->LogCAMu[s]->sVal2[gofw_Mean] = par->HMu[m][i] / N;
         res->LogCAMu[s]->pVal2[gofw_Mean] =
            fbar_Normal1 (par->HMu[m][i] / N);
         res->LogCASig_sVal[s] = par->HSig[m][i] / N;
         res->LogCASig_pVal[s] = fbar_ChiSquare2 (N, 12, par->HSig[m][i]);
      } else {
         res->LogCAMu[s]->sVal2[gofw_Mean] = 0.0;
         res->LogCAMu[s]->pVal2[gofw_Mean] = 0.0;
         res->LogCASig_sVal[s] = 0.0;
         res->LogCASig_pVal[s] = 0.0;
         memset (res->LogCAMu[s]->sVal2, 0, sizeof (res->LogCAMu[s]->sVal2));
         memset (res->LogCAMu[s]->pVal2, 0, sizeof (res->LogCAMu[s]->pVal2));
      }

      i = SQUA_EXACT_CIRC;
      k = i + 8 * par->Loc[m];
      tables_CopyTabD (res->Collectors[k]->V, res->SquareCEMu[s]->sVal1->V,
         1, N);
      res->SquareCEMu[s]->sVal1->NObs = N;
      gofw_ActiveTests2 (res->Collectors[k]->V, res->SquareCEMu[s]->pVal1->V,
         N, wdist_Normal, (double *) NULL, res->SquareCEMu[s]->sVal2,
         res->SquareCEMu[s]->pVal2);

      res->SquareCEMu[s]->sVal2[gofw_Mean] = par->HMu[m][i] / N;
      res->SquareCEMu[s]->pVal2[gofw_Mean] =
         fbar_Normal1 (par->HMu[m][i] / N);
      res->SquareCESig_sVal[s] = par->HSig[m][i] / N;
      res->SquareCESig_pVal[s] = fbar_ChiSquare2 (N, 12, par->HSig[m][i]);

      i = SQUA_ASYMP_CIRC;
      k = i + 8 * par->Loc[m];
      tables_CopyTabD (res->Collectors[k]->V, res->SquareCAMu[s]->sVal1->V,
         1, N);
      res->SquareCAMu[s]->sVal1->NObs = N;
      if (flag) {
         gofw_ActiveTests2 (res->Collectors[k]->V,
            res->SquareCAMu[s]->pVal1->V, N, wdist_Normal, (double *) NULL,
            res->SquareCAMu[s]->sVal2, res->SquareCAMu[s]->pVal2);
         res->SquareCAMu[s]->sVal2[gofw_Mean] = par->HMu[m][i] / N;
         res->SquareCAMu[s]->pVal2[gofw_Mean] =
            fbar_Normal1 (par->HMu[m][i] / N);
         res->SquareCASig_sVal[s] = par->HSig[m][i] / N;
         res->SquareCASig_pVal[s] = fbar_ChiSquare2 (N, 12, par->HSig[m][i]);

      } else {
         res->SquareCAMu[s]->sVal2[gofw_Mean] = 0.0;
         res->SquareCAMu[s]->pVal2[gofw_Mean] = 0.0;
         res->SquareCASig_sVal[s] = 0.0;
         res->SquareCASig_pVal[s] = 0.0;
         memset (res->SquareCAMu[s]->sVal2, 0,
            sizeof (res->SquareCAMu[s]->sVal2));
         memset (res->SquareCAMu[s]->pVal2, 0,
            sizeof (res->SquareCAMu[s]->pVal2));
      }

      if (M0 == 0 && m == 1)
         m = D;
      else
         m += D;
      s++;
   }
}


/*=========================================================================*/

void sspacings_AllSpacings (unif01_Gen * gen, sspacings_Res * res,
   long N, long n0, int r, int M0, int M1, int D, int LgEps)
{
   long i;
   int m, s;
   long Seq;
   double Eps;                    /* Minimal spacing for SumLogsSpacings */
   double *U;
   double Prod, x;
   double LnProd;
   double SumSq;
   int NbMinus[MAXM];             /* Number of spacings < Eps */
   Param par;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sspacings_AllSpacings test";

   Timer = chrono_Create ();
   memset (&par, 0, sizeof (Param));
   par.NbColl = 4;
   InitAllSpacings (gen, TestName, &par, N, n0, r, M0, M1, D, LgEps);
   Eps = 1.0 / num_TwoExp[LgEps];

   if (res == NULL) {
      localRes = TRUE;
      res = sspacings_CreateRes ();
   }
   InitRes (res, N, par.Nbm, "sspacings_AllSpacings");
   res->step = 2;

   U = util_Calloc ((size_t) n0 + M1 + 2, sizeof (double));
   U[0] = 0.0;

   /* Beginning of tests */
   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 1; i <= n0; i++)
         U[i] = unif01_StripD (gen, r);
      tables_QuickSortD (U, 1, n0);
      util_Assert (U[1] >= 0.0, "sspacings_AllSpacings:   U[1] < 0.0");
      util_Assert (U[n0] <= 1.0, "sspacings_AllSpacings:   U[n] > 1.0");

      for (i = 1; i <= M1; i++)
         U[n0 + i] = 1.0 + U[i - 1];
      if (M0 == 0)
         m = 1;
      else
         m = M0;

      while (m <= M1) {
         NbMinus[m] = 0;
         Prod = 1.0;
         LnProd = 0.0;
         SumSq = 0.0;

         for (i = 0; i <= n0; i++) {
            x = U[i + m] - U[i];
            SumSq += x * x;
            /* In case a spacing is zero */
            if (x < Eps) {
               x = Eps;
               ++NbMinus[m];
            }
            /* Compute log of product instead of sum of logs; it is faster */
            Prod *= x;
            if (Prod < Epsilon) {
               /* Take log once in a while to avoid underflow */
               LnProd += log (Prod);
               Prod = 1.0;
            }
         }
         LnProd += log (Prod);

         UpdateStat (&par, m, LOG_EXACT_CIRC, LnProd, res->Collectors);
         UpdateStat (&par, m, LOG_ASYMP_CIRC, LnProd, res->Collectors);
         UpdateStat (&par, m, SQUA_EXACT_CIRC, SumSq * n0 * n0,
            res->Collectors);
         UpdateStat (&par, m, SQUA_ASYMP_CIRC, SumSq * n0 * n0,
            res->Collectors);
         if (M0 == 0 && m == 1)
            m = D;
         else
            m += D;
      }
   }

   CopyResults (res, &par, N, M0, M1, D, 1);
   if (swrite_Basic) {
      printf ("\nResults:");
      if (M0 == 0)
         m = 1;
      else
         m = M0;
      s = 0;
      while (m <= M1) {
         printf ("\n----------------------------------------------------\n");
         printf ("m = %1d\n\n", m);
         if (NbMinus[m] > 0)
            printf ("%1d spacings < 1 / 2^%1d\n\n", NbMinus[m], LgEps);

         printf ("Logs of spacings:\n-----------------\n\n");
         WrRes ("Exact mean and standard deviation, circular:",
            N, &par, m, LOG_EXACT_CIRC, res->Collectors,
            res->LogCEMu[s]->sVal2, res->LogCEMu[s]->pVal2);
         WrRes ("Asymptotic mean and standard deviation, circular:",
            N, &par, m, LOG_ASYMP_CIRC, res->Collectors,
            res->LogCAMu[s]->sVal2, res->LogCAMu[s]->pVal2);

         printf ("\nSquares of spacings:\n--------------------\n\n");
         WrRes ("Exact mean and standard deviation, circular:",
            N, &par, m, SQUA_EXACT_CIRC, res->Collectors,
            res->SquareCEMu[s]->sVal2, res->SquareCEMu[s]->pVal2);
         WrRes ("Asymptotic mean and standard deviation, circular:",
            N, &par, m, SQUA_ASYMP_CIRC, res->Collectors,
            res->SquareCAMu[s]->sVal2, res->SquareCAMu[s]->pVal2);

         if (M0 == 0 && m == 1)
            m = D;
         else
            m += D;
         s++;
      }
      printf ("\n");
      swrite_Final (gen, Timer);
   }

   U = util_Free (U);
   if (localRes)
      sspacings_DeleteRes (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void sspacings_AllSpacings2 (unif01_Gen * gen, sspacings_Res * res,
   long N, long n0, int r, int M0, int M1, int D, int LgEps)
{
   long i;
   int m, j, s;
   long Seq;
   double Eps;                    /* Minimal spacing for SumLogsSpacings */
   double *U;
   double Prod, x;
   double LnProd;
   double SumSq;
   int NbMinus[MAXM];             /* Number of spacings < Eps */
   Param par;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sspacings_AllSpacings2 test";

   Timer = chrono_Create ();
   memset (&par, 0, sizeof (Param));
   par.NbColl = 2;
   InitAllSpacings (gen, TestName, &par, N, n0, r, M0, M1, D, LgEps);
   Eps = 1.0 / num_TwoExp[LgEps];

   if (res == NULL) {
      localRes = TRUE;
      res = sspacings_CreateRes ();
   }
   InitRes (res, N, par.Nbm, "sspacings_AllSpacings2");
   res->step = 4;

   U = util_Calloc ((size_t) n0 + M1 + 2, sizeof (double));
   U[0] = 0.0;

   /* Beginning of tests */
   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 1; i <= n0; i++)
         U[i] = unif01_StripD (gen, r);
      tables_QuickSortD (U, 1, n0);
      for (j = 1; j <= M1; j++)
         U[n0 + j] = 1.0 + U[j - 1];
      if (M0 == 0)
         m = 1;
      else
         m = M0;

      while (m <= M1) {
         NbMinus[m] = 0;
         Prod = 1.0;
         LnProd = 0.0;
         SumSq = 0.0;

         for (i = 0; i <= n0; i++) {
            x = U[i + m] - U[i];
            SumSq += x * x;
            /* In case a spacing is zero */
            if (x < Eps) {
               x = Eps;
               ++NbMinus[m];
            }
            /* Compute log of product instead of sum of logs; it is faster */
            Prod *= x;
            if (Prod < Epsilon) {
               /* Take log once in a while to avoid underflow */
               LnProd += log (Prod);
               Prod = 1.0;
            }
         }
         LnProd += log (Prod);

         UpdateStat (&par, m, LOG_EXACT_CIRC, LnProd, res->Collectors);
         UpdateStat (&par, m, SQUA_EXACT_CIRC, SumSq * n0 * n0,
            res->Collectors);
         if (M0 == 0 && m == 1)
            m = D;
         else
            m += D;
      }
   }

   CopyResults (res, &par, N, M0, M1, D, 0);
   if (swrite_Basic) {
      printf ("\nResults:");
      if (M0 == 0)
         m = 1;
      else
         m = M0;
      s = 0;

      while (m <= M1) {
         printf ("\n----------------------------------------------------\n");
         printf ("m = %1d\n\n", m);
         if (NbMinus[m] > 0)
            printf ("%1d spacings < 1 / 2^%1d\n\n", NbMinus[m], LgEps);
         printf ("Logs of spacings:\n-----------------\n\n");
         WrRes ("Exact mean and standard deviation, circular:",
            N, &par, m, LOG_EXACT_CIRC, res->Collectors,
            res->LogCEMu[s]->sVal2, res->LogCEMu[s]->pVal2);
         printf ("\nSquares of spacings:\n--------------------\n\n");
         WrRes ("Exact mean and standard deviation, circular:",
            N, &par, m, SQUA_EXACT_CIRC, res->Collectors,
            res->SquareCEMu[s]->sVal2, res->SquareCEMu[s]->pVal2);
         if (M0 == 0 && m == 1)
            m = D;
         else
            m += D;
         s++;
      }
      printf ("\n");
      swrite_Final (gen, Timer);
   }

   U = util_Free (U);
   if (localRes)
      sspacings_DeleteRes (res);
   chrono_Delete (Timer);
}
