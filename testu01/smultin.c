/*************************************************************************\
 *
 * Package:        TestU01
 * File:           smultin.c
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
#include "num2.h"
#include "tables.h"

#include "smultin.h"
#include "wdist.h"
#include "swrite.h"
#include "unif01.h"

#include "statcoll.h"
#include "gofw.h"
#include "fmass.h"
#include "fdist.h"
#include "fbar.h"

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <float.h>






/*============================= constants ===============================*/

/* Max string length */
#define LENGTH 100

/* Max dimension t */
#define MAX_DIM 64

/* Used for secondary hashing */
#define HACHE2 41

/* Upper Limit of precomputed tables of 2nI, when SPARSE = TRUE */
static const long LIM_SPARSE = 64;

/* LIM_DENSE*n/k = upper limit of the precomputed tables for 2nI, when */
/* SPARSE = FALSE and 6n/k < k */
static const double LIM_DENSE = 6.0;

/* Precision with which we measure specific values of ValDelta */
static const double EPS_LAM = 1.0E-14;

#ifdef USE_LONGLONG
#define MAXK 9223372036854775808.0  /* 2^63 */
#else
#define MAXK 9007199254740992.0     /* 2^53 */
#endif

/* Our gamma distribution is not good for parameters larger than this. */
#define EMPTYLIM 200000000000000.0

#define MASK64  0x8000000000000000ULL  /* 2^63: set bit 64 to 1 */

#define TRACE(N)  printf ("*********   "#N"%13ld\n", N);



/*=============================== Types =================================*/

/* Index for the different values of ValDelta */
typedef int DeltaIndex;




/*========================= Extern variables ============================*/

smultin_Envir smultin_env = {
   MAXK,
   1024 * 1024,                   /* SeuilHash */
   0.75,                          /* HashLoad */
   5.0E+6,                        /* SeuilEColl */

   12.0,                          /* SeuilCOverDense */
   5.0,                           /* SeuilCOverNorSup */
   1.999,                         /* SeuilCOverNorInf */
   1.0001                         /* SeuilCOverSparse */
};


/* The stable parameters values used by default */
smultin_Param smultin_ParamDefault = {
   2,                             /* NbDelta */
   {-1, 1},                       /* ValDelta */
   smultin_GenerCellSerial,       /* GenerCell */
   -3                             /* bmax */
};



/*============================== Functions ==============================*/


smultin_Param *smultin_CreateParam (int NbDelta, double ValDelta[],
   smultin_GenerCellType GenerCell, int bmax)
{
   smultin_Param *par;
   int j;

   par = util_Malloc (sizeof (smultin_Param));
   par->NbDelta = NbDelta;
   for (j = 0; j < NbDelta; j++) {
      util_Assert (ValDelta[j] >= -1.0,
         "smultin_CreateParam:   ValDelta[j] < -1");
      par->ValDelta[j] = ValDelta[j];
   }
   util_Assert (bmax <= smultin_MAXB,
      "smultin_CreateParam:   bmax > smultin_MAXB");
   par->bmax = bmax;
   par->GenerCell = GenerCell;
   return par;
}

/*-------------------------------------------------------------------------*/

void smultin_DeleteParam (smultin_Param *par)
{
   if (par == NULL)
      return;
   util_Free (par);
}


/*=========================================================================*/

static void CleanPD (smultin_Res *res)
{
   DeltaIndex s;

   if (res == NULL)
      return;

   for (s = 0; s < res->NbDeltaOld; s++) {
      res->TabFj[s] = util_Free (res->TabFj[s]);
   }

   res->Count = util_Free (res->Count);
   res->Count1 = util_Free (res->Count1);
   res->Cell = util_Free (res->Cell);
   res->Cell1 = util_Free (res->Cell1);
   res->Nb = util_Free (res->Nb);
   res->Nb1 = util_Free (res->Nb1);
}


/*=========================================================================*/

static void InitRes (
   smultin_Param *par, 
   smultin_Res *res,          /* Results holder */
   long N                     /* Number of replications */
)
/* 
 * Initializes the smultin_Res structure. The old NbDelta is in res, the
 * new NbDelta is in par. Delete the unused collectors if old NbDelta > new
 * NbDelta, and create the needed collectors if old NbDelta < new NbDelta.
 */
{
   DeltaIndex s;

   if (par == NULL)
      par = &smultin_ParamDefault;
   CleanPD (res);

   for (s = par->NbDelta; s < res->NbDeltaOld; s++)
      res->Collector[s] = statcoll_Delete (res->Collector[s]);

   for (s = res->NbDeltaOld; s < par->NbDelta; s++)
      res->Collector[s] = statcoll_Create (N, "");

   for (s = 0; s < par->NbDelta; s++) {
      statcoll_Init (res->Collector[s], N);
      gofw_InitTestArray (res->sVal2[s], -1.0);
      gofw_InitTestArray (res->pVal2[s], -1.0);
   }

   res->NbDeltaOld = par->NbDelta;
   res->flagTab = FALSE;
   res->nLimit = 1;
   res->pColl = res->pEmpty = -1.0;
   res->pCollLeft = -1.0;
   res->pCollRight = -1.0;
}


/*-------------------------------------------------------------------------*/

smultin_Res * smultin_CreateRes (smultin_Param *par)
{
   smultin_Res *res;
   DeltaIndex s;

   res = util_Malloc (sizeof (smultin_Res));
   memset (res, 0, sizeof (smultin_Res));

   if (par == NULL)
      par = &smultin_ParamDefault;

   for (s = 0; s < par->NbDelta; s++) {
      res->Collector[s] = statcoll_Create (1, "");
      res->TabFj[s] = NULL;
   }

   res->Count = NULL;
   res->Count1 = NULL;
   res->Cell = NULL;
   res->Cell1 = NULL;
   res->NbDeltaOld = 0;
   res->Nb = NULL;
   res->Nb1 = NULL;

   return res;
}


/*-------------------------------------------------------------------------*/

void smultin_DeleteRes (smultin_Res * res)
{
   DeltaIndex s;

   if (res == NULL)
      return;

   for (s = 0; s < res->NbDeltaOld; s++)
      res->Collector[s] = statcoll_Delete (res->Collector[s]);

   CleanPD (res);
   util_Free (res);
}


/*=========================================================================*/

double smultin_MNTermeKhi2 (double junk, double NbEsp, long j)
/*
 * One term of Chi2 = 2nI1; there are j balls in this urn.
 */
{
   double Diff;
   Diff = j - NbEsp;
   return Diff * Diff / NbEsp;
}


/*=======================================================================*/

double smultin_MNTermePowDiv (double Delta, double NbEsp, long j)
/*
 * One term of Power Divergence = 2nI; there are j balls in this urn.
 */
{
   double y;
   if (j == 0)
      return 0.0;
   y = pow (j / NbEsp, Delta) - 1.0;
   return (2.0 * j * y) / (Delta * (Delta + 1.0));
}


/*=======================================================================*/

double smultin_MNTermeLogLikhood (double junk, double NbEsp, long j)
/*
 * One term of loglikelihood ratio = 2nI0; there are j balls in this urn
 */
{
   if (j == 0)
      return 0.0;
   return 2.0 * j * log (j / NbEsp);
}


/*=======================================================================*/

double smultin_MNTermeColl (double junk1, double junk2, long j)
/*
 * Number of collisions when there are j balls in this urn
 */
{
   if (j <= 1)
      return 0.0;
   return (double) (j - 1);
}


/*=======================================================================*/

static void MNCalcMuSigma (
   double V[],           /* Contains the terms of the statistic */
   long nlim,            /* Limit on the non negligible terms */
   long n,               /* Number of balls */
   double k,             /* Number of urns */
   double *Mu,           /* Mean */
   double *Sigma         /* Standard deviation */
   )
/*
 * Compute the mean Mu and the standard deviation Sigma. For some values
 * of n and k, some of the terms to be subtracted will be huge and loss
 * of precision may give negative variance. We shall then stop the tests
 * and exit in WriteDataPowDiv.
 */
{
   const double Epsilon = 1.0E-100;    /* To avoid division by 0 */
   const double Eps = 1.0E-18;
   long i, j, Mid;
   double cond2, cond1, temp;
   double Sum3, Sum2, Sum1;
   double muk, Var;
   double km;
   double Terme2, TermeMid, Terme1;
   double x, nr = n;

   util_Assert (nlim <= n, "MNCalcMuSigma;  nlim > n");
   Mid = n / k;
   while (Mid < nlim && fabs (V[Mid]) < Eps)
      ++Mid;
   util_Assert (Mid <= nlim, "MNCalcMuSigma;  Mid > nlim");

   /******* Compute Mean x = *Mu */
   TermeMid = fmass_BinomialTerm3 (n, 1.0 / k, Mid) * k;
   x = V[Mid] * TermeMid;
   Terme1 = x;
   i = Mid;
   km = k - 1;
   while (i < n && fabs (Terme1 / x) > Eps) {
      util_Assert (i < nlim, "MNCalcMuSigma: nlim too small --> Espion001");
      Terme1 *= V[i + 1] * (n - i) / (V[i] * (i + 1) * km);
      x += Terme1;
      ++i;
   }
   Terme1 = V[Mid] * TermeMid;
   i = Mid;
   while (i > 0 && fabs (Terme1) / x > Eps) {
      Terme1 *= V[i - 1] * i * km / (V[i] * (n - i + 1));
      x += Terme1;
      --i;
   }
   /* The cases when |Terme[i]| < Eps, but |Terme[i-1]| > Eps */
   /* these terms must be included also. */
   --i;
   if (i >= 0) {
      Terme1 = V[i] * fmass_BinomialTerm3 (n, 1.0 / k, i) * k;
      x += Terme1;
      while (i > 0 && fabs (Terme1 / x) > Eps) {
         Terme1 = Terme1 * V[i - 1] * i * km / (V[i] * (n - i + 1));
         x += Terme1;
         --i;
      }
   }

   /****** Calculate variance: first series of terms */
   muk = x / k;
   Terme1 = (V[Mid] - muk) * (V[Mid] - muk) * TermeMid;
   cond2 = TermeMid * V[Mid] * V[Mid];
   Sum1 = Terme1;
   Terme2 = Terme1;
   i = Mid;
   km = k - 1;
   while (i < n && fabs (cond2 / Sum1) + fabs (Terme2 / Sum1) > Eps) {
      util_Assert (i < nlim, "MNCalcMuSigma: nlim too small --> Espion002");
      temp = V[i + 1] / (V[i] + Epsilon);
      cond2 *= temp * temp * (n - i) / ((i + 1) * km);
      temp = (V[i + 1] - muk) / (V[i] - muk);
      Terme2 *= temp * temp * (n - i) / ((i + 1) * km);
      Sum1 += Terme2;
      ++i;
   }
   Terme2 = Terme1;
   cond2 = TermeMid * V[Mid] * V[Mid];
   i = Mid;
   while (i > 0 && fabs (cond2 / Sum1) + fabs (Terme2 / Sum1) > Eps) {
      temp = V[i - 1] / (V[i] + Epsilon);
      cond2 *= temp * temp * i * km / (n - i + 1);
      temp = (V[i - 1] - muk) / (V[i] - muk);
      Terme2 *= temp * temp * i * km / (n - i + 1);
      Sum1 += Terme2;
      --i;
   }

   /****** Calculate variance: terms i = j */
   temp = fmass_BinomialTerm4 (n - Mid, 1.0 / k, 2.0 / k, Mid) *
      fmass_BinomialTerm4 (n, 1.0 / k, 0.0, Mid) * k * (k - 1);
   TermeMid = temp * (V[Mid] - muk) * (V[Mid] - muk);
   cond1 = TermeMid * V[Mid] * V[Mid] 
             / ((V[Mid] - muk) * (V[Mid] - muk) + Epsilon);
   Terme1 = TermeMid;
   Sum2 = Terme1;
   i = Mid;
   km = k - 2;
   while ((i < n / 2)
      && fabs (cond1 / Sum2) + fabs (Terme1 / Sum2) > Eps) {
      util_Assert (i < nlim, "MNCalcMuSigma: nlim too small --> Espion003");
      temp = V[i + 1] / ((V[i] + Epsilon) * (i + 1) * km);
      cond1 *= (nr - 2 * i) * (nr - 2 * i - 1) * temp * temp;
      temp = (V[i + 1] - muk) / ((V[i] - muk) * (i + 1) * km);
      Terme1 *= (nr - 2 * i) * (nr - 2 * i - 1) * temp * temp;
      Sum2 += Terme1;
      ++i;
   }
   i = Mid;
   Terme1 = TermeMid;
   cond1 = TermeMid * V[Mid] * V[Mid] / ((V[Mid] - muk) * (V[Mid] - muk));
   while (i > 0 && fabs (cond1 / Sum2) + fabs (Terme1 / Sum2) > Eps) {
      temp = (V[i - 1] * i * km) / (V[i] + Epsilon);
      cond1 *= temp * temp / ((nr - 2 * i + 2) * (nr - 2 * i + 1) + Epsilon);
      temp = (V[i - 1] - muk) * i * km / (V[i] - muk);
      Terme1 *= temp * temp / ((nr - 2 * i + 2) * (nr - 2 * i + 1) + Epsilon);
      Sum2 += Terme1;
      --i;
   }
   /****** Calculate variance: terms i <> j */
   i = Mid + 1;
   Sum3 = 1.0E-40;
   cond1 = 1.0;
   Terme1 = 1.0;
   while (i <= n && fabs (cond1 / Sum3) + fabs (Terme1 / Sum3) > Eps) {
      util_Assert (i <= nlim, "MNCalcMuSigma: nlim too small --> Espion004");
      j = Mid;
      if (j > n - i)
         j = n - i;
      temp = fmass_BinomialTerm4 (n, 1.0 / k, 0.0, i) *
         fmass_BinomialTerm4 (n - i, 1.0 / k, 2.0 / k, j);
      Terme1 = temp * (V[i] - muk) * (V[j] - muk) * k * (k - 1);
      cond1 = Terme1 * V[i] * V[j] / ((V[i] - muk) * (V[j] - muk));
      Sum3 += Terme1;
      Terme2 = Terme1;
      cond2 = cond1;
      while (j > 0 && fabs (cond2 / Sum3) + fabs (Terme2 / Sum3) > Eps) {
         cond2 *= V[j - 1] * j * km / (V[j] * (n - i - j + 1) + Epsilon);
         Terme2 *= (V[j - 1] - muk) * j * km /
                    ((V[j] - muk) * (n - i - j + 1) + Epsilon);
         Sum3 += Terme2;
         --j;
      }
      Terme2 = Terme1;
      cond2 = cond1;
      j = Mid;
      while ((j < i - 1 && i + j < n)
         && fabs (cond2 / Sum3) + fabs (Terme2 / Sum3) > Eps) {
         util_Assert (j < nlim, "MNCalcMuSigma: nlim too small --> Espion005");
         cond2 *= V[j + 1] * (n - i - j) / (V[j] * (j + 1) * km + Epsilon);
         Terme2 *= (V[j + 1] - muk) * (n - i - j)
                     / ((V[j] - muk) * (j + 1) * km + Epsilon);
         Sum3 += Terme2;
         ++j;
      }
      ++i;
   }
   i = Mid;
   cond1 = Sum3 + 1.0;
   Terme1 = 1.0;
   while (i > 0 && fabs (cond1 / Sum3) + fabs (Terme1 / Sum3) > Eps) {
      j = i - 1;
      temp = fmass_BinomialTerm4 (n, 1.0 / k, 0.0, i) *
         fmass_BinomialTerm4 (n - i, 1.0 / k, 2.0 / k, j);
      Terme1 = temp * (V[i] - muk) * (V[j] - muk) * k * (k - 1.0);
      cond1 = (Terme1 * V[i] * V[j]) / ((V[i] - muk) * (V[j] - muk));
      Sum3 += Terme1;
      Terme2 = Terme1;
      cond2 = cond1;
      while (j > 0 && fabs (cond2 / Sum3) + fabs (Terme2 / Sum3) > Eps) {
         cond2 *= V[j - 1] * j * km / (V[j] * (n - i - j + 1) + Epsilon);
         Terme2 *= (V[j - 1] - muk) * j * km /
                       ((V[j] - muk) * (n - i - j + 1) + Epsilon);
         Sum3 += Terme2;
         --j;
      }
      --i;
   }
   Var = Sum1 + Sum2 + 2.0 * Sum3;
   util_Warning (Var < 0.0, "MNCalcMuSigma:   negative variance");
   if (Var >= 0.0)
      *Sigma = sqrt (Var);
   else
      *Sigma = -1.0;
   *Mu = x;
}


/*=======================================================================*/

void smultin_MultinomMuSigma (
   long n,                    /* Number of balls */
   double k,                  /* Number of urns */
   double Theta1,             /* First parameter of the term F */
   double Theta2,             /* Second parameter of the term F */
   smultin_MNTermeType F,     /* One term of the statistic */
   double *Mu,                /* Mean */
   double *Sigma              /* Standard deviation */
   )
/*
 * Compute the mean Mu and the standard deviation Sigma
 */
{

   /* For densities n/k < 8, only the ~ 25 first terms will contribute */
   /* significantly to the normal approximation; thus we precompute only */
   /* elements [0..LIM_SPARSE] of the tables in the case Sparse = TRUE */

   long nlim;
   long j;
   double densite;
   double *PV;

   /* We may choose n >>> nlim because the probabilities will be concen- */
   /* trated near j = 0 for low densites ( < 8). Large values of j will  */
   /* practically never occur. It is not necessary to compute all the    */
   /* PV[0..n]. For high densities n/k, nlim will have to be increased.  */

   densite = n / (double) k;
   nlim = 8 * densite;
   if (nlim < LIM_SPARSE)
      nlim = LIM_SPARSE;          /* Sparse = TRUE */
   if (nlim > n)
      nlim = n;                   /* Sparse = FALSE */
   PV = util_Calloc ((size_t) nlim + 2, sizeof (double));
   for (j = 0; j <= nlim; j++)
      PV[j] = F (Theta1, Theta2, j);
   MNCalcMuSigma (PV, nlim, n, k, Mu, Sigma);
   util_Free (PV);
}


/*=======================================================================*/

static void CalcTabFj (
   smultin_Param *par,
   smultin_Res *res,
   lebool Sparse,
   double k,                  /* Number of cells or urns */
   double NbExp               /* Expected number per cell */
   )
/*
 * May pre-calculate all terms. Will then calculate all non negligible
 * terms for all values of s, and keep them in tables TabFj[s][].
 */
{
   long i;
   DeltaIndex s;
   double delta;
   double c;
   double temp;
   double *F;

   if (!Sparse && LIM_DENSE * NbExp > k) {
      /* Do not precompute tables when Sparse = FALSE and we have a very */
      /* small number k of cells */
      res->flagTab = FALSE;
      return;
   }

   /* Precompute the values and keep them in arrays */
   res->flagTab = TRUE;
   if (Sparse)
      res->nLimit = LIM_SPARSE;
   else {
      res->nLimit = LIM_DENSE * NbExp;
      if (res->nLimit < 1)
         res->nLimit = 2;
   }

   for (s = 0; s < par->NbDelta; s++) {
      res->TabFj[s] = util_Calloc (2 + (size_t) res->nLimit, sizeof (double));
      delta = par->ValDelta[s];
      util_Assert (delta >= -1.0 - EPS_LAM,
         "CalcTabFj:   par->ValDelta[s] < -1");
      F = res->TabFj[s];
      F[0] = 0.0;

      if (fabs (delta - 1.0) < EPS_LAM) {
         /* ChiSquare */
         for (i = 0; i <= res->nLimit; i++) {
            temp = i - NbExp;
            F[i] = temp * temp / NbExp;
         }

      } else if (fabs (delta) < EPS_LAM) {
         /* LogLikelyhood */
         for (i = 1; i <= res->nLimit; i++) {
            temp = i;
            F[i] = 2.0 * temp * log (temp / NbExp);
         }

      } else if (fabs (delta + 1.0) < EPS_LAM) {
         /* Collision */
         for (i = 1; i <= res->nLimit; i++) {
            F[i] = i - 1;
         }

      } else {
         /* PowerDivergence, delta > -1 */
         c = 2.0 / (delta * (delta + 1.0));
         for (i = 1; i <= res->nLimit; i++) {
            temp = i;
            F[i] = c * temp * (pow (temp / NbExp, delta) - 1.0);
         }
      }
   }
}


/*=======================================================================*/

static void ReCalcTabFj (
   smultin_Param *par,
   smultin_Res *res,
   double NbExp              /* Expected number per cell */
   )
/*
 * Update tables TabFj when one of the Count becomes larger than res->nLimit
 */
{
   long i;
   DeltaIndex s;
   double delta;
   double c;
   double temp;
   double *F;
   long i0 = res->nLimit;
   res->nLimit *= 2;

   for (s = 0; s < par->NbDelta; s++) {
      delta = par->ValDelta[s];
      res->TabFj[s] = util_Realloc (res->TabFj[s],
         (res->nLimit + 1) * sizeof (double));
      F = res->TabFj[s];

      if (fabs (delta - 1.0) < EPS_LAM) {
         /* ChiSquare */
         for (i = i0 + 1; i <= res->nLimit; i++) {
            temp = i - NbExp;
            F[i] = temp * temp / NbExp;
         }

      } else if (fabs (delta) < EPS_LAM) {
         /* LogLikelyhood */
         for (i = i0 + 1; i <= res->nLimit; i++) {
            temp = i;
            F[i] = 2.0 * temp * log (temp / NbExp);
         }

      } else if (fabs (delta + 1.0) < EPS_LAM) {
         /* Collision Test */
         for (i = i0 + 1; i <= res->nLimit; i++) {
            F[i] = i - 1;
         }

      } else {
         c = 2.0 / (delta * (delta + 1.0));
         for (i = i0 + 1; i <= res->nLimit; i++) {
            temp = i;
            F[i] = c * temp * (pow (temp / NbExp, delta) - 1.0);
         }
      }
   }
}


/*=======================================================================*/

void smultin_PowDivMomCorChi (
   double Delta,
   long n,                    /* Number of balls */
   double k,                  /* Number of urns */
   double *MuC,               /* Corrected mean */
   double *SigmaC             /* Corrected standard deviation */
   )
/*
 * Compute the corrected mean and standard deviation in the dense case
 *  for the ChiSquare approximation. (See Read and Cressie)
 */
{
   double t = k * k;
   double temp;
   if (Delta < EPS_LAM - 1.0) {
      *MuC = -1.0;
      *SigmaC = -1.0;
      return;
   }
   temp = (8.0 - 12.0 * k - 2.0 * k * k + 6.0 * t +
      (Delta - 1.0) * (4.0 - 6.0 * k - 3.0 * k * k + 5.0 * t) / 3.0
      + 2.0 * (Delta - 2.0) * (1.0 - 2.0 * k + t));
   *SigmaC = 2.0 - 2.0 * k - (double) k * k + t + (Delta - 1.0) * temp;
   *SigmaC = sqrt (1.0 + *SigmaC / (2.0 * n * (k - 1.0)));
   temp = (2.0 - 3.0 * k + t) / 3.0 +
      (Delta - 2.0) * (1.0 - 2.0 * k + t) / 4.0;
   *MuC = (k - 1.0) * (1.0 - *SigmaC) + (Delta - 1.0) * temp / n;
}


/*=======================================================================*/

void smultin_PowDivMom (
   double Delta,              /* Which Power Divergence */
   long n,                    /* Number of balls */
   double k,                  /* Number of urns */
   double NbExp,              /* Expected number per urn */
   double *Mu,                /* Mean */
   double *Sigma              /* Standard deviation */
   )
/*
 * Compute the mean and standard deviation in the sparse case
 */
{

   if ((double) n / k > 8.0) {
      printf ("*************  Call of smultin_PowDivMom with n/k > 8\n");
      *Mu = -1.0;
      *Sigma = -1.0;
      return;
   }
   if (k <= 2) {
      printf ("*************  Call of smultin_PowDivMom with k <= 2\n");
      *Mu = -1.0;
      *Sigma = -1.0;
      return;
   }

   util_Assert ((double) n / k <= 8.0,
      "smultin: Call of PowDivMom with n/k > 8");
   util_Assert (k > 2, "smultin: Call of PowDivMom with k <= 2");

   if (fabs (Delta - 1.0) < EPS_LAM) {
      /* ChiSquare test */
      *Mu = k - 1;
      *Sigma = sqrt (2.0 * (k - 1) * (n - 1.0) / n);

   } else if (fabs (Delta + 1.0) < EPS_LAM) {
      /* Collision test */
      smultin_MultinomMuSigma (n, k, 0.0, 0.0, smultin_MNTermeColl, Mu, Sigma);

   } else if (fabs (Delta) < EPS_LAM) {
      /* Delta = 0, LogLikelyhood */
      smultin_MultinomMuSigma (n, k, 0.0, NbExp, smultin_MNTermeLogLikhood,
                               Mu, Sigma);

   } else if (Delta > -1.0) {
      smultin_MultinomMuSigma (n, k, Delta, NbExp, smultin_MNTermePowDiv,
                               Mu, Sigma);

   } else
      util_Error ("smultin_PowDivMom:   Delta < -1.0");
}


/*=======================================================================*/
#if 0

static void CalcPowDiv (double Delta, double NbExp[], long Count[],
   long smin, long smax, double *X)
/*
 * Compute the statistic $2n I^\delta$ defined in (\ref{powdiv}),
 * for $\delta = {\tt Delta}$, and return its value in $X$.
 * We assume that the expected values $n p_i$ in cell $i$  are in
 * {\tt NbExp[smin..smax]}, the observed values $X_i$ are in
 * {\tt Count[smin..smax]}, {\tt smin} and  {\tt smax} are the indices 
 * of the first and the last cell,
 * and the number of cells is $k = {\tt smax} - {\tt smin} + 1$.
 * The $p_i$ are not necessarily equal.
 */
{
   double temp;
   long s;
   *X = 0.0;

   if (fabs (Delta - 1.0) < EPS_LAM) {
      /* ChiSquare */
      for (s = smin; s <= smax; s++) {
         temp = Count[s] - NbExp[s];
         *X += (temp * temp) / NbExp[s];
      }

   } else if (fabs (Delta) < EPS_LAM) {
      /* Loglikelihood */
      for (s = smin; s <= smax; s++) {
         if (Count[s] > 0) {
            temp = Count[s];
            *X += temp * log (temp / NbExp[s]);
         }
      }
      *X *= 2.0;

   } else if (Delta <= EPS_LAM - 1.0) {
      util_Error ("smultin_CalcPowDiv:   Delta <= -1.0");
      /* We do the collisions test only when probabilities are equal. See
         smultin_CalcPowDivEqual. */

   } else {
      /* Other values of Delta.  */
      for (s = smin; s <= smax; s++) {
         if (Count[s] > 0) {
            temp = Count[s];
            *X += temp * (pow (temp / NbExp[s], Delta) - 1.0);
         }
      }
      *X = 2.0 * *X / (Delta * (Delta + 1.0));
   }
}

#endif

/*=======================================================================*/

static void CalcPowDivEqual (
   smultin_Param *par,
   smultin_Res *res,
   DeltaIndex s,
   double NbExp,                  /* Expected number per cell */
   long Count[],                  /* Counters */
   long jmin,                     /* First cell */
   long jmax,                     /* Last cell */
   lebool flagTab,                /* TRUE: use precomputed table */
   double *X                      /* Computed statistic */
   )
/*
 * This function is called only when we do not use hashing.
 *
 * As in CalcPowDiv, except that the values of np_i are all equal to
 * np = NbExp. The lebool flagTab indicates if the values of
 * $2 i ln (i/np)$ and
 * $ {2\over \delta(1+\delta)}
 *   i \left[\left(i/np\right)^\delta -1\right]$ have already been
 * pre-computed and kept in the arrays
 * F (it is so if {\tt flagTab = TRUE}).
 */
{
   double temp;
   double *F = res->TabFj[s];     /* The precomputed table */
   double Delta = par->ValDelta[s];
   long j;
   *X = 0.0;

   if (flagTab) {
      /* For low densities, we use precomputed tables since the observed */
      /* values will all be very small. We shall thus need only a few */
      /* terms (expensive to compute) of the statistics. */
      util_Assert (res->nLimit > 0,
         "smultin_CalcPowDivEqual BUG: res->nLimit <= 0");

      for (j = jmin; j <= jmax; j++) {
         /* A larger than expected counter needs terms that have not been */
         /* precomputed: Recompute missing terms. */
         while (Count[j] > res->nLimit) {
            ReCalcTabFj (par, res, NbExp);
            F = res->TabFj[s];
         }
         *X += F[Count[j]];
      }
      return;
   }

   /* High densities: no precomputed tables */
   if (fabs (Delta - 1.0) < EPS_LAM) {
      /* ChiSquare: Delta = 1 */
      for (j = jmin; j <= jmax; j++) {
         temp = Count[j] - NbExp;
         *X += temp * temp;
      }
      *X /= NbExp;

   } else if (fabs (Delta) < EPS_LAM) {
      /* Loglikelihood ratio */
      for (j = jmin; j <= jmax; j++) {
         if (Count[j] > 0) {
            temp = Count[j];
            *X += temp * log (temp / NbExp);
         }
      }
      *X *= 2.0;

   } else if (fabs (Delta + 1.0) < EPS_LAM) {
      /* Collision test */
      for (j = jmin; j <= jmax; j++) {
         if (Count[j] > 1)
            *X += Count[j] - 1;
      }

   } else if (Delta > (-1.0)) {
      for (j = jmin; j <= jmax; j++) {
         if (Count[j] > 0) {
            temp = Count[j];
            *X += temp * (pow (temp / NbExp, Delta) - 1.0);
         }
      }
      *X = (2.0 * *X) / (Delta * (Delta + 1.0));

   } else
      util_Error ("smultin_CalcPowDivEqual: Delta < -1");
}


/*=======================================================================*/

static void CalcPoDiEqHache (
   smultin_Param *par,
   smultin_Res *res, 
   DeltaIndex i,
   double NbExp,              /* Expected number per cell */
   smultin_CellType Nb[],     /* Number of cells with s balls */
   long CountMax,             /* Max number of balls in any cell */
   lebool flagTab,           /* TRUE if use precomputed table */
   double *X                  /* Computed statistic */
   )
/*
 * Compute the Power Divergence statistic or the number of collisions in 
 * the sparse case. We use a hashing table of smultin_CellType.
 */
{
   double temp;
   double *F = res->TabFj[i];     /* The precomputed table */
   double Delta = par->ValDelta[i];
   long s;
   *X = 0.0;

   if (flagTab) {
      /* For low densities, we use precomputed tables since the observed */
      /* values will all be very small. We shall thus need only a few */
      /* terms (expensive to compute) of the statistics. */
      util_Assert (res->nLimit > 0, "CalcPoDiEqHache BUG: res->nLimit <= 0");

      /* A larger than expected counter needs terms that have not been */
      /* precomputed: Recompute missing terms. */
      while (CountMax > res->nLimit) {
         ReCalcTabFj (par, res, NbExp);
         F = res->TabFj[i];
      }

      for (s = 0; s <= CountMax; s++) {
         *X += F[s] * Nb[s];
      }
      return;
   }

   /* High densities: no precomputed tables */
   if (fabs (Delta - 1.0) < EPS_LAM) {
      /* ChiSquare: Delta = 1 */
      for (s = 1; s <= CountMax; s++) {
         temp = s - NbExp;
         *X += temp * temp * Nb[s];
      }
      *X = *X / NbExp + NbExp * Nb[0];

   } else if (fabs (Delta) < EPS_LAM) {
      /* Delta = 0: Loglikelihood ratio */
      for (s = 1; s <= CountMax; s++) {
         temp = s;
         *X += temp * log (temp / NbExp) * Nb[s];
      }
      *X *= 2.0;

   } else if (fabs (Delta + 1.0) < EPS_LAM) {
      /* Collision test */
      for (s = 2; s <= CountMax; s++) {
         *X += (s - 1.0) * Nb[s];
      }

   } else if (Delta > -1.0) {
      for (s = 1; s <= CountMax; s++) {
         temp = s;
         *X += temp * (pow (temp / NbExp, Delta) - 1.0) * Nb[s];
      }
      *X = 2.0 * *X / (Delta * (Delta + 1.0));

   } else
      util_Error ("CalcPoDiEqHache: Delta < -1");
}


/*=======================================================================*/

static void CalcNbCells (
   smultin_Param *par,
   smultin_Res *res,
   long jmin,                 /* First cell */
   long jmax,                 /* Last cell */
   long CoMax                 /* Maximum number of balls in any cell */
   )
/*
 * Compute the number of cells containing j balls or more.
 */
{

   long j;
   smultin_CellType wb[smultin_MAXB + 1];
   long *Count = res->Count;      /* Counters */
   smultin_CellType *Nb = res->Nb; /* Nb[j] = number of cells with j balls */

   util_Assert (par->bmax <= smultin_MAXB,
      "CalcNbCells:   smultin_MAXB is too small");

   for (j = 0; j <= smultin_MAXB; j++)
      wb[j] = 0;

   if (res->Hashing) {
      for (j = smultin_MAXB; j <= CoMax; j++)
         wb[smultin_MAXB] += Nb[j];
      for (j = smultin_MAXB - 1; j >= 0; j--)
         wb[j] = wb[j + 1] + Nb[j];

   } else {
      Nb[0] = 0;
      for (j = jmin; j <= jmax; j++) {
         if (Count[j] > smultin_MAXB) {
            wb[smultin_MAXB] += 1;
         } else
            Nb[Count[j]] += 1;
      }
      wb[smultin_MAXB] += Nb[smultin_MAXB];
      for (j = smultin_MAXB - 1; j >= 0; j--)
         wb[j] = wb[j + 1] + Nb[j];
   }

   /* the local array wb is necessary in the case N > 1 and Poisson since */
   /* then, the statistic used is the sum of the N  Poisson statistics;   */
   /* we are here summing the numbers for the N replications of the test. */

   for (j = 0; j <= smultin_MAXB; j++) {
      res->WbCells[j] += wb[j];
      res->NbCells[j] += Nb[j];
   }
}


/*=======================================================================*/

smultin_CellType smultin_GenerCellSerial (unif01_Gen *gen,
    int r, int t, long d)
{
   int j;
   smultin_CellType dr = d;
   smultin_CellType Cell;

   Cell = unif01_StripL (gen, r, d);
   for (j = 2; j <= t; j++)
      Cell = Cell * dr + unif01_StripL (gen, r, d);
   return Cell;
}


/*=======================================================================*/

smultin_CellType smultin_GenerCellSerial2 (unif01_Gen *gen,
   int r, int t, long d)
{
   int j;
   smultin_CellType dr = d;
   smultin_CellType Cell;

   Cell = unif01_StripL (gen, r, d);
   for (j = 2; j <= t; j++) {
      Cell += dr * unif01_StripL (gen, r, d);
      dr *= d;
   }
   return Cell;
}


/*=======================================================================*/

smultin_CellType smultin_GenerCellPermut (unif01_Gen *gen,
   int r, int t, long junk)
{
   int s, i, j;
   smultin_CellType Cell = 0;
   double U[64];

   for (j = 1; j <= t; j++)
      U[j] = unif01_StripD (gen, r);

   for (i = t; i >= 2; i--) {
      /* Find the U[s] = max (U[1],...,U[i]) */
      s = 1;
      for (j = 2; j <= i; j++) {
         if (U[j] > U[s])
            s = j;
      }
      Cell = Cell * i + (s - 1);
      U[s] = U[i];
   }
   return Cell;
}


/*=======================================================================*/

smultin_CellType smultin_GenerCellMax (unif01_Gen *gen,
   int r, int t, long junk)
{
   int i, MaxI;
   double U, MaxU = -1.0;

   /* Don't forget that cells are numbered from 0 to k - 1 */
   for (i = 0; i < t; i++) {
      U = unif01_StripD (gen, r);
      if (U > MaxU) {
         MaxU = U;
         MaxI = i;
      }
   }
   return (smultin_CellType) MaxI;
}


/*=======================================================================*/

smultin_CellType smultin_GenerCellSerialBits (unif01_Gen * gen,
   int r, int s, long L)
{
   const int t = L / s;
   const smultin_CellType dr = num_TwoExp[s];
   smultin_CellType Cell;
   int j;

   Cell = unif01_StripB (gen, r, s);
   for (j = 2; j <= t; j++)
      Cell = Cell * dr + unif01_StripB (gen, r, s);
   return Cell;
}


/*=======================================================================*/

fmass_INFO smultin_CreateCollisions (long n, smultin_CellType k)
{
   const long nLim = 100000;
   const int MaxIter = 32;
   const double Epsilon = DBL_EPSILON;
   const double DensityLim = 1.0001;
   long J1, J0, j, i, Dim;
   double terme, v, u, mu, sigma, x;
   double kinv = 1.0 / k;
   double *A;
   fmass_INFO W;

   util_Assert (k > 0, "smultin_CreateCollisions:  k <= 0");
   util_Assert (n > 0, "smultin_CreateCollisions:  n <= 0");

   /* Poisson Approximation */
   if ((n > nLim) && ((double) n / k <= DensityLim)) {
      if ((double) n / k <= 0.1) {
         int jj;
         /* To avoid loss of precision when n/k --> 0, we expand the */
         /* formula below in a MacLaurin series */

         jj = 3;
         u = n - 1;
         v = 2.0;
         terme = (n * u) / (2.0 * k * k);
         mu = terme;
         while (fabs (terme / mu) > Epsilon && jj < MaxIter) {
            u -= 1.0;
            v += 1.0;
            terme = -terme * u / (k * v);
            mu += terme;
            ++jj;
         }
         util_Assert (jj < MaxIter,
                      "smultin_CreateCollisions: limit MaxIter hit");

      } else if (n <= 100) {
         mu = ((double) n / k - 1.0) + pow (1.0 - 1.0 / k, (double) n);

      } else {
         const int ITER = 10;
         int i;
         terme = kinv;
         mu = terme;

         /* Compute the log of pow(1 - 1/k, n) by Maclaurin series */
         for (i = 2; i < ITER; i++) {
            terme *= kinv;
            mu += terme / i;
         }
         mu = ((double) n / k - 1.0) + exp (-n * mu);
      }

      mu *= k;
      W = fmass_CreatePoisson (mu);
      /* W->paramR[0] now contains the Poisson parameter mu */
      W->paramR = util_Realloc (W->paramR, 3 * sizeof (double));
      W->paramR[1] = n;
      W->paramR[2] = k;
      W->paramI = util_Malloc (sizeof (long));
      W->paramI[0] = smultin_CollPoissonSparse;
      return W;
   }

   W = util_Malloc (sizeof (struct fmass_INFO_T));
   W->paramI = util_Malloc (sizeof (long));
   W->paramR = util_Calloc (5, sizeof (double));
   W->paramR[1] = n;
   W->paramR[2] = k;


   /* Normal Approximation */
   if (n > nLim) {
      smultin_MultinomMuSigma (n, (double) k, 0.0, 0.0, smultin_MNTermeColl,
                               &mu, &sigma);
      W->paramR[3] = mu;
      W->paramR[4] = sigma;
      W->paramI[0] = smultin_CollNormal;
      W->pdf = NULL;
      W->cdf = NULL;
      W->smin = -1;
      W->smax = -1;
      return W;
   }


   /* Exact Distribution */
   A = util_Calloc ((size_t) n + 2, sizeof (double));
   for (j = 0; j <= n; j++)
      A[j] = 0.0;
   A[1] = 1.0;
   J1 = J0 = 1;
   for (j = 1; j <= n - 1; j++) {
      ++J1;
      i = J1;
      while (i >= J0) {
         x = i * kinv;
         A[i] = x * A[i] + (1.0 + kinv - x) * A[i - 1];
         if (A[i] <= Epsilon) {
            A[i] = 0.0;
            if (i == J1)
               --J1;
            else if (i == J0)
               ++J0;
         }
         --i;
      }
   }
   Dim = n - J0 + 1;
   W->pdf = util_Calloc ((size_t) Dim + 1, sizeof (double));
   W->cdf = util_Calloc ((size_t) Dim + 1, sizeof (double));

   W->pdf[0] = A[n];
   W->cdf[0] = A[n];
   j = 0;
   while (j < Dim && W->cdf[j] < 1.0) {
      ++j;
      W->pdf[j] = A[n - j];
      W->cdf[j] = W->pdf[j] + W->cdf[j - 1];
   }
   while (j <= Dim) {
      W->pdf[j] = A[n - j];
      W->cdf[j] = 1.0;
      ++j;
   }
   util_Free (A);
   W->paramI[0] = smultin_CollExact;
   W->smin = 0;
   W->smax = Dim;
   return W;
}


/*-------------------------------------------------------------------------*/

void smultin_DeleteCollisions (fmass_INFO W)
{
   if (W == NULL)
      return;
   util_Free (W->paramI);
   util_Free (W->paramR);
   util_Free (W->pdf);
   util_Free (W->cdf);
   util_Free (W);
}


/*-------------------------------------------------------------------------*/

double smultin_CollisionsTerm (fmass_INFO W, long s)
{
   int par;
   double z;
   double Mu;
   double Sigma;

   util_Assert (W != NULL,
      "smultin_CollisionsTerm:   fmass_INFO is NULL pointer");
   if (s < 0)
      return 0.0;
   par = W->paramI[0];

   switch (par) {
   case smultin_CollPoissonSparse:
      return fmass_PoissonTerm2 (W, s);
   case smultin_CollNormal:
      Mu = W->paramR[3];
      Sigma = W->paramR[4];
      z =  fdist_Normal2 ((s - Mu) / Sigma) -
           fdist_Normal2 ((s - 1 - Mu) / Sigma);
      return z;
   case smultin_CollExact:
      if (s > W->smax)
         return 0.0;
      return W->pdf[s];
   default:
      util_Error ("smultin_CollisionsTerm:  Not initialized");
      return 0.0;
   }
}


/*-------------------------------------------------------------------------*/

double smultin_FDistCollisions (fmass_INFO W, long s)
{
   int par;

   util_Assert (W != NULL,
      "smultin_FDistCollisions: fmass_INFO is NULL pointer");
   if (s < 0)
      return 0.0;
   par = W->paramI[0];

   switch (par) {
   case smultin_CollPoissonSparse:
      return fdist_Poisson2 (W, s);
   case smultin_CollNormal:
      /* W->paramR[3] = Mu, W->paramR[4] = Sigma */
      return fdist_Normal2 ((s - W->paramR[3]) / W->paramR[4]);
   case smultin_CollExact:
      if (s > W->smax)
         return 1.0;
      return W->cdf[s];
   default:
      util_Error ("smultin_FDistCollisions:  Not initialized");
      return 0.0;
   }
}


/*-------------------------------------------------------------------------*/

double smultin_FBarCollisions (fmass_INFO W, long s)
{
   return 1.0 - smultin_FDistCollisions (W, s - 1);
}


/*=======================================================================*/

static void InitPowDiv (
   smultin_Param *par,
   smultin_Res *res,
   long N,                    /* Number of replications */
   lebool Sparse,
   long n,                    /* Number of balls */
   smultin_CellType z         /* Number of urns (not quite for PowDivOver) */
   )
/*
 * Initialize the multinomial tests
 */
{

   DeltaIndex s;
   long j;
   double NbExp;
   char chaine[LENGTH + 1];
   char Str[LENGTH + 1];
   double Mu;                     /* Mean */
   double Sigma;                  /* Standard Deviation */

   NbExp = (double) n / z;
   if (z >= smultin_env.SeuilHash && NbExp < 1.0)
      res->Hashing = TRUE;
   else
      res->Hashing = FALSE;

   res->EsCells[0] = N * (double) z * exp (-NbExp);
   res->EsEmpty = res->EsCells[0];
   res->NbCells[0] = 0;
   res->WbCells[0] = 0;

   util_Assert (par->NbDelta <= smultin_MAX_DELTA,
      "par->NbDelta > smultin_MAX_DELTA");
   for (s = 0; s < par->NbDelta; s++) {

      if (Sparse) {
         smultin_PowDivMom (par->ValDelta[s], n, (double) z, (double) n / z,
            &Mu, &Sigma);

      } else if (fabs (par->ValDelta[s] + 1.0) > EPS_LAM) {
         /* Non collision tests */
         smultin_PowDivMomCorChi (par->ValDelta[s], n, (double) z, &Mu,
            &Sigma);

      } else {
         /* Meaningless values as flags */
         Mu = -1.0;
         Sigma = -1.0;
      }
      res->Mu[s] = Mu;
      res->Sigma[s] = Sigma;

      if (fabs (par->ValDelta[s] + 1.0) < EPS_LAM) {
         /* Collision test */
         strncpy (Str, "The N statistic values for Collision:", (size_t) 64);
         res->NbCollisions = 0.0;
         for (j = 1; j <= smultin_MAXB; j++) {
            res->NbCells[j] = 0;
            res->WbCells[j] = 0;
         }
         /* The exact expected numbers (from Knuth)
             temp = n * log ((z - 1.0)/z);
             res->EsCells[0] = z * exp (temp);
             res->EsCells[1] = n * res->EsCells[0] / (z - 1.0);
             res->EsCells[2] = (n - 1.0) * res->EsCells[1] / (2.0*(z - 1.0));
             res->EsCells[>=3] = z - res->EsCells[0] - res->EsCells[1] -
                                  res->EsCells[2];
         */

         /* Expected numbers of urns with exactly j balls in the Poisson */
         /* approximation */
         for (j = 1; j <= smultin_MAXB; j++) {
            res->EsCells[j] = (res->EsEmpty * pow (NbExp, (double) j))
               / num2_Factorial (j);
         }
         /* Expected numbers of urns with >= j balls */
         for (j = smultin_MAXB - 1; j >= 0; j--) {
            res->EsCells[j] += res->EsCells[j + 1];
         }

      } else {
         /* Non Collision tests */
         strncpy (Str, "The N statistic values for Delta = ", (size_t) 64);
         sprintf (chaine, "%4.2f:", par->ValDelta[s]);
         strncat (Str, chaine, (size_t) 10);
      }

      statcoll_SetDesc (res->Collector[s], Str);
   }
}


/*=======================================================================*/

static void WriteDataPowDiv (
   unif01_Gen *gen,
   smultin_Param *par,
   smultin_Res *res, 
   char *TestName,
   long N,                    /* Number of replications */
   long n,                    /* Number of balls */
   int r,                     /* Drop r bits from each random number */
   long d,                    /* Number of segments on 1-dimensional line */
   int t,                     /* Dimension */
   lebool Sparse,
   smultin_CellType k         /* Number of urns */
)
/*
 * Write the parameters of the test
 */
{
   double EC;
   double NbExp;
   DeltaIndex s;

   swrite_Head (gen, TestName, N, n, r);

   if (par->GenerCell == smultin_GenerCellSerial) {
      printf (",   d = %4ld,   t = %2d,\n       Sparse = ", d, t);
      util_WriteBool (Sparse, 6);
      printf ("\n\n");
      printf ("       GenerCell = smultin_GenerCellSerial\n");
      printf ("       Number of cells = d^t = ");
   } else if (par->GenerCell == smultin_GenerCellSerial2) {
      printf (",   d = %4ld,   t = %2d,\n       Sparse = ", d, t);
      util_WriteBool (Sparse, 6);
      printf ("\n\n");
      printf ("       GenerCell = smultin_GenerCellSerial2\n");
      printf ("       Number of cells = d^t = ");
   } else if (par->GenerCell == smultin_GenerCellPermut) {
      printf (",   t = %2d,\n       Sparse = ", t);
      util_WriteBool (Sparse, 6);
      printf ("\n\n");
      printf ("       GenerCell = smultin_GenerCellPermut\n");
      util_Assert (!res->Over,
         "MultinomialOver: non implemented for smultin_GenerCasePermut");
      printf ("       Number of cells = t! = ");
   } else if (par->GenerCell == smultin_GenerCellMax) {
      printf (",   k = %2d,\n       Sparse = ", t);
      util_WriteBool (Sparse, 6);
      printf ("\n\n");
      printf ("       GenerCell = smultin_GenerCellMax\n");
      printf ("       Number of cells = k = ");
   }

#ifdef USE_LONGLONG
   printf ("%18" PRIuLEAST64 "\n", k);
#else
   printf ("%18.0f\n", k);
#endif

   util_Assert (k <= smultin_env.Maxk, "Multinomial:  k is too large");
   printf ("       Expected number per cell =  ");
   NbExp = (double) n / k;
   if (NbExp < 1.0)
      printf ("1 / %10.8g\n", 1.0 / NbExp);
   else
      printf ("%10.8g\n", NbExp);

   EC = (double) n * n / (2.0 * k);
   if (Sparse)
      printf ("       EColl = n^2 / (2k) = %12.10g\n", EC);
   printf ("       Hashing = ");
   util_WriteBool (res->Hashing, 6);
   printf ("\n\n");
   if (par->NbDelta == 1 && par->ValDelta[0] == -1)
      ;
   else {
      if (Sparse) {
         printf ("   For Delta > -1, we use the normal approximation\n");
         printf ("   Mean and standard deviation: \n");
      } else {
         printf ("   For Delta > -1, we use the ChiSquare approximation\n");
         printf ("   Correction factor of the ChiSquare: \n");
      }
   }

   for (s = 0; s < par->NbDelta; s++) {
      if (fabs (par->ValDelta[s] + 1.0) < EPS_LAM) {
         if ((Sparse == FALSE) && (res->Over == FALSE)) {
            /* The collision test is meaningless when Sparse = FALSE */
            continue;
         }
         printf ("       Collision test");
      } else {
         printf ("       Delta = %5.2g", par->ValDelta[s]);
      }
      if (!res->Over || fabs (par->ValDelta[s] + 1.0) > EPS_LAM) {
         printf (",    Mu = %14.8g", res->Mu[s]);
         printf (",    Sigma = %10.6g", res->Sigma[s]);
         util_Assert (res->Sigma[s] > 0.0, "Negative Variance");
      }
      printf ("\n");
   }
   printf ("\n");
}


/*=======================================================================*/

static void WriteDataMNBits (
   unif01_Gen *gen,
   smultin_Param *par,
   smultin_Res *res, 
   char *TestName,
   long N,                    /* Number of replications */
   long n,                    /* Number of balls */
   int r,                     /* Drop r bits from each random number */
   long L,                    /* Number of bits for a cell */
   int s,                     /* Number of bits taken from each rand. num. */
   lebool Sparse,
   smultin_CellType k,        /* Number of cells = 2^L */
   lebool Over               /* Overlapping case = TRUE */
)
/*
 * Write the parameters of the test
 */
{

   double EC;
   double NbExp;
   DeltaIndex j;

   swrite_Head (gen, TestName, N, n, r);

   printf (",   s = %2d,   L = %4ld,\n       Sparse = ", s, L);
   util_WriteBool (Sparse, 6);
   if (Over)
      printf ("\n\n       Number of bits = n = %1ld\n", n);
   else
      printf ("\n\n       Number of bits = n*L = %1ld\n", L * n);

   /* printf (" GenerCell = smultin_GenerCellSerialBits\n"); */

#ifdef USE_LONGLONG
   printf ("       Number of cells = 2^L = %18" PRIuLEAST64 "\n", k);
#else
   printf ("       Number of cells = 2^L = %18.0f\n", k);
#endif
   util_Assert (k <= smultin_env.Maxk, "Multinom:  k is too large");

   printf ("       Expected number per cell =  ");
   NbExp = (double) n / k;
   if (NbExp < 1.0)
      printf ("1 / %10.8g\n", 1.0 / NbExp);
   else
      printf ("%10.8g\n", NbExp);

   EC = (double) n * n / (2.0 * k);
   if (Sparse)
      printf ("       EColl = n^2 / (2k) = %12.10g\n", EC);
   printf ("       Hashing = ");
   util_WriteBool (res->Hashing, 6);
   printf ("\n\n");
   if (par->NbDelta == 1 && par->ValDelta[0] == -1)
      ;
   else {
      if (Sparse) {
         printf ("   For Delta > -1, we use the normal approximation\n");
         printf ("   Mean and standard deviation: \n");
      } else {
         printf ("   For Delta > -1, we use the ChiSquare approximation\n");
         printf ("   Correction factor of the ChiSquare: \n");
      }
   }

   for (j = 0; j < par->NbDelta; j++) {
      if (fabs (par->ValDelta[j] + 1.0) < EPS_LAM) {
         if ((Sparse == FALSE) && (res->Over == FALSE)) {
            /* The collision test is meaningless when Sparse = FALSE */
            continue;
         }
         printf ("       Collision test");
      } else {
         printf ("       Delta = %5.2g", par->ValDelta[j]);
      }
      if (!res->Over || fabs (par->ValDelta[j] + 1.0) > EPS_LAM) {
         printf (",    Mu = %14.8g", res->Mu[j]);
         printf (",    Sigma = %10.6g\n", res->Sigma[j]);
         util_Assert (res->Sigma[j] > 0.0, "Negative Variance");
      }
   }
   printf ("\n");
}


/*=======================================================================*/

static void CalcResultsPowDiv (
   smultin_Param *par,
   smultin_Res *res, 
   DeltaIndex s,                  /* Which statistic */
   long n,                        /* Number of balls */
   lebool Sparse,
   smultin_CellType DegreLib,     /* Number of degrees of freedom */
   double Mu,                     /* Mean */
   double SumX[],
   double SumX2[]
   )
{

   double pR, pL;
   double pCollLeft;              /* Left p-value of Collision test */
   double pCollRight;             /* Right p-value of Collision test */
   int j;
   statcoll_Collector *SC = res->Collector[s];
   fmass_INFO Mass1, Mass2;
   double V[1];
   long N = SC->NObs;
   double racN;

   if (fabs (par->ValDelta[s] + 1.0) < EPS_LAM) {
      /* Collision Test */
      if (!Sparse) {
         res->pColl = -1.0;
         res->pCollLeft = -1.0;
         res->pCollRight = -1.0;
         return;
      }
      /* The total number of collisions of N replications is NbCollisions */
      Mass1 = smultin_CreateCollisions (n, DegreLib + 1);
      if (N == 1) {
         pCollLeft = smultin_FDistCollisions (Mass1, (long) res->NbCollisions);
         pCollRight = smultin_FBarCollisions (Mass1, (long) res->NbCollisions);
         res->pCollLeft = pCollLeft;
         res->pCollRight = pCollRight;
         res->sVal2[s][gofw_Mean] = SumX[s];
         res->pVal2[s][gofw_Mean] = fbar_Normal1 (res->sVal2[s][gofw_Mean]);
      } else {
         if (Mu < smultin_env.SeuilEColl) {
            Mass2 = fmass_CreatePoisson (N * Mu);
            pCollLeft = fdist_Poisson2 (Mass2, (long) res->NbCollisions);
            pCollRight = fbar_Poisson2 (Mass2, (long) res->NbCollisions);
            fmass_DeletePoisson (Mass2);
         }
      }
      smultin_DeleteCollisions (Mass1);
      res->pColl = gofw_pDisc (pCollLeft, pCollRight);

      /* Total number of empty urns of the N replications: res->NbCells[0]
      Mass2 = fmass_CreatePoisson (res->EsEmpty);
      pL = fdist_Poisson2 (Mass2, res->NbCells[0]);
      pR = fbar_Poisson2 (Mass2, res->NbCells[0]);
      fmass_DeletePoisson (Mass2);
      res->pEmpty = gofw_pDisc (pL, pR); */

      /* Since we can have very large values of res->EsEmpty (2^63), we
         compute the Poisson pL, pR by calling the Gamma distribution, since
         our Poisson takes a long argument */

      if (res->NbCells[0] <= EMPTYLIM && res->EsEmpty <= EMPTYLIM) {
	 pL = fbar_Gamma (res->NbCells[0] + 1.0, 12, res->EsEmpty);
	 if ((res->NbCells[0] <= 0)  || (res->NbCells[0] > res->NbCellsTotal))
	    pR = 1.0;
	 else
	    pR = fdist_Gamma ((double) (res->NbCells[0]), 12, res->EsEmpty);
	 res->pEmpty = gofw_pDisc (pL, pR);
      }
      /* The total number of urns containing >= j balls */
      for (j = 2; j <= par->bmax; j++) {
         Mass2 = fmass_CreatePoisson ((double) (res->EsCells[j]));
         pL = fdist_Poisson2 (Mass2, (long) res->WbCells[j]);
         pR = fbar_Poisson2 (Mass2, (long) res->WbCells[j]);
         fmass_DeletePoisson (Mass2);
         res->pWb[j] = gofw_pDisc (pL, pR);
      }

   } else if (Sparse) {
      /* Tests other than Collision test */
      gofw_ActiveTests1 (SC->V, N, wdist_Normal, (double *) NULL,
                         res->sVal2[s], res->pVal2[s]);

   } else {
      V[0] = DegreLib;
      gofw_ActiveTests1 (SC->V, N, wdist_ChiSquare, V,
                         res->sVal2[s], res->pVal2[s]);
   }

   /* Now compute the mean and the correlation with their p-values. */
   if (N > 1) {
      racN = sqrt ((double) N);
      res->sVal2[s][gofw_Mean] = SumX[s] / racN;
      res->pVal2[s][gofw_Mean] = fbar_Normal1 (res->sVal2[s][gofw_Mean]);
      res->sVal2[s][gofw_Cor] = racN * SumX2[s] / (N - 1);
      res->pVal2[s][gofw_Cor] = fbar_Normal1 (res->sVal2[s][gofw_Cor]);
   }
}


/*=======================================================================*/

static void WriteResultsPowDiv (
   smultin_Param *par,
   smultin_Res *res, 
   DeltaIndex s,
   long N,
   double EColl,              /* Approximate expected number of collisions */
   smultin_CellType DegreLib, /* Number of degrees of freedom */
   lebool Sparse,
   double Mu                  /* Exact expected mean */
   )
{
   long j;
   printf ("-----------------------------------------------\n");
   printf ("Test Results for ");

   if (fabs (par->ValDelta[s] + 1.0) < EPS_LAM)
      printf ("Collisions\n\n");
   else {
      printf ("Delta = %8.4f\n\n", par->ValDelta[s]);
      if (N == 1) {
         if (!Sparse) {
#ifdef USE_LONGLONG
            printf ("Number of degrees of freedom          : %4" PRIuLEAST64
               "\n", DegreLib);
#else
            printf ("Number of degrees of freedom          : %4.0f\n",
               DegreLib);
#endif

         }
         printf ("Value of the statistic                :");
         gofw_Writep2 (res->sVal2[s][gofw_Mean], res->pVal2[s][gofw_Mean]);

      } else {
         gofw_WriteActiveTests0 (N, res->sVal2[s], res->pVal2[s]);
         printf ("For the sum of the N observations, we use\n");
         printf ("      the Normal approximation:\n");
         printf ("Standardized empirical mean           :");
         gofw_Writep2 (res->sVal2[s][gofw_Mean], res->pVal2[s][gofw_Mean]);
         printf ("Standardized empirical correlation    :");
         gofw_Writep2 (res->sVal2[s][gofw_Cor], res->pVal2[s][gofw_Cor]);
      }
   }

   if (swrite_Collectors) {
      if (fabs (par->ValDelta[s] + 1.0) < EPS_LAM)
         statcoll_Write (res->Collector[s], 5, 14, 0, 0);
      else
         statcoll_Write (res->Collector[s], 5, 14, 4, 3);
   }

   if (fabs (par->ValDelta[s] + 1.0) < EPS_LAM) {
      if (N > 1 && Mu < smultin_env.SeuilEColl) {
         printf ("For the total number of collisions, we use\n"
            "      the Poisson approximation:\n");
         /* "Value of N * EColl : "); num_WriteD(N * EColl, 11, 2, 2); */
         printf ("Expected number of collisions = N*Mu  : ");
         num_WriteD (N * Mu, 11, 2, 2);
         printf ("\nObserved number of collisions         : %8ld\n",
            (long) res->NbCollisions);
         gofw_Writep1 (res->pColl);
         printf ("\n");
      }
      if (N == 1) {
	/*     printf ("Value of EColl                        : ");
	       num_WriteD(EColl, 11, 2, 2);*/
         printf ("Expected number of collisions = Mu    : ");
         num_WriteD (Mu, 11, 2, 2);
         printf ("\nObserved number of collisions         : %8ld\n",
            (long) res->NbCollisions);
         gofw_Writep1 (res->pColl);
      }
      printf ("-----------------------------\n");
      printf ("Total number of cells containing j balls\n\n");
      for (j = 0; j <= smultin_MAXB / 2; j++) {
         printf ("  j = %2ld", j);

#ifdef USE_LONGLONG
         printf ("                              : %16" PRIuLEAST64 "\n",
            res->NbCells[j]);
#else
         printf ("                              : %16.0f\n", res->NbCells[j]);
#endif
      }

      if (par->bmax >= 0 && res->NbCells[0] <= EMPTYLIM &&
            res->EsEmpty <= EMPTYLIM) {
         printf ("\n-----------------------------\n"
            "Results for the number of empty cells\n\n"
            "Expected number                       : ");
         num_WriteD (res->EsEmpty, 19, 2, 2);
         printf ("\nObserved number                       :");
#ifdef USE_LONGLONG
         printf (" %16" PRIuLEAST64 "\n", res->NbCells[0]);
#else
         printf (" %16.0f\n", res->NbCells[0]);
#endif
         gofw_Writep1 (res->pEmpty);
      }
      if (par->bmax >= 1) {
         printf ("\n-----------------------------\n");
         printf ("Results for the number of cells containing at least"
            " j balls\n\n");
         for (j = 2; j <= par->bmax; j++) {
            printf ("  j = %2ld\n", j);
            printf ("Expected number                       : %11.2f\n",
               (double) (res->EsCells[j]));
            printf ("Observed number                       : %8.0f\n",
               (double) (res->WbCells[j]));
            gofw_Writep1 (res->pWb[j]);
         }
      }
   }
   printf ("\n");
}


/*=======================================================================*/

static void UpdateCountHash (
   smultin_Res *res, 
   smultin_CellType Ind,
   long Hache,
   double UnSurHache,
   long *CoMax,
   lebool DimFlag          /* TRUE for t-1 dimension, FALSE for t dim. */
   )
/*
 * We use hashing. A ball falls in cell Ind: update counters
 * Speed is essential here.
 */
{
   long *Count;
   smultin_CellType *Cell;
   smultin_CellType *Nb;
   long Decal, Pos, Tem;

   if (DimFlag == FALSE) {
      Count = res->Count;         /* Counters in t dimensions */
      Cell = res->Cell;           /* Cell numbers in t dimensions */
      Nb = res->Nb;
   } else {
      Count = res->Count1;        /* Counters in t - 1 dimensions */
      Cell = res->Cell1;          /* Cell numbers in t - 1 dimensions */
      Nb = res->Nb1;
   }

#ifdef USE_LONGLONG
   Pos = Ind % Hache;
#else
   Tem = Ind * UnSurHache;
   Pos = Ind - (double) Hache * Tem;
#endif

   Decal = HACHE2 + Pos % HACHE2;

   /* Insert in hashing table; if sign bit is 1, cell is empty. */
   for (;;) {
#ifdef USE_LONGLONG
      if (Cell[Pos] & MASK64) {
#else
      if (Cell[Pos] < 0.0) {
#endif
         Cell[Pos] = Ind;
         break;
      }
      if (Cell[Pos] == Ind)
         break;
      Pos = (Pos + Decal) % Hache;
   }

   Nb[Count[Pos]] -= 1;
   ++(Count[Pos]);
   if (Count[Pos] > *CoMax)
      ++(*CoMax);
   if (DimFlag == FALSE) {
      if (*CoMax > res->NbSize) {
         int i;
         res->NbSize *= 2;
         res->Nb = util_Realloc (res->Nb,
            (res->NbSize + 1) * sizeof (smultin_CellType));
         Nb = res->Nb;
         for (i = res->NbSize / 2 + 1; i <= res->NbSize; i++)
            Nb[i] = 0;
      }
   } else {
      if (*CoMax > res->Nb1Size) {
         int i;
         res->Nb1Size *= 2;
         res->Nb1 = util_Realloc (res->Nb1,
            (res->Nb1Size + 1) * sizeof (smultin_CellType));
         Nb = res->Nb1;
         for (i = res->Nb1Size / 2 + 1; i <= res->Nb1Size; i++)
            Nb[i] = 0;
      }
   }
   Nb[Count[Pos]] += 1;
}


/*=======================================================================*/

static void GenerAllPointsHash (unif01_Gen * gen, smultin_Param * par,
   smultin_Res * res, long n, int r, long d, int t, long *pCoMax,
   long Hache, double UnSurHache)
/*
 * Generate all n points in hashing case
 */
{
   smultin_CellType Indice;       /* Cell number */
   long i;

   for (i = 0; i <= Hache; i++)
#ifdef USE_LONGLONG
      res->Cell[i] = MASK64;      /* Empty cells */
#else
      res->Cell[i] = -1.0;        /* Empty cells */
#endif
   *pCoMax = 0;
   for (i = 1; i <= n; i++) {
      Indice = par->GenerCell (gen, r, t, d);
      UpdateCountHash (res, Indice, Hache, UnSurHache, pCoMax, FALSE);
   }
}


/*=======================================================================*/

static void GenerAllPoints2 (unif01_Gen * gen, smultin_Param * par,
   smultin_Res * res, long n, int r, long d, int t)
/*
 * Generate all n points; no hashing
 */
{
   smultin_CellType Indice;       /* Cell number */
   long i;
   for (i = 1; i <= n; i++) {
      Indice = par->GenerCell (gen, r, t, d);
      ++res->Count[(long) Indice];
   }
}


/*=======================================================================*/

static void GenerAllPointsHashBits (unif01_Gen *gen, smultin_Res *res,
   long n, int r, long L, int s, long *pCoMax, long Hache,
   double UnSurHache)
/*
 * Generate all n points of L bits each in hashing case
 */
{
   smultin_CellType Indice;       /* Cell number */
   long i;
   int j;
   unsigned long Z;
   const int t = s / L;           /* Number of points in a U01 */
   const long Last = n % t;
   const unsigned long MASK = num_TwoExp[L] - 1.0;

   for (i = 0; i <= Hache; i++)
#ifdef USE_LONGLONG
      res->Cell[i] = MASK64;      /* Empty cells */
#else
      res->Cell[i] = -1.0;        /* Empty cells */
#endif
   *pCoMax = 0;

   for (i = 1; i <= n / t; i++) {
      Z = unif01_StripB (gen, r, s);
      for (j = 1; j <= t; j++) {
         Indice = Z & MASK;
         UpdateCountHash (res, Indice, Hache, UnSurHache, pCoMax, FALSE);
         Z >>= L;
      }
   }
   /* The last points */
   if (Last > 0) {
      Z = unif01_StripB (gen, r, s);
      /* The most significant bits make the points */
      for (j = 1; j <= t - Last; j++)
         Z >>= L;
      for (j = 1; j <= Last; j++) {
         Indice = Z & MASK;
         UpdateCountHash (res, Indice, Hache, UnSurHache, pCoMax, FALSE);
         Z >>= L;
      }
   }
}


/*=======================================================================*/

static void GenerAllPoints2Bits (unif01_Gen * gen, smultin_Res * res,
   long n, int r, long L, int s)
/*
 * Generate all n points of L bits each; no hashing
 */
{
   long i;
   int j;
   unsigned long Z;
   const int t = s / L;
   const long Last = n % t;
   const unsigned long MASK = num_TwoExp[L] - 1.0;

   for (i = 1; i <= n / t; i++) {
      Z = unif01_StripB (gen, r, s);
      for (j = 1; j <= t; j++) {
         ++res->Count[Z & MASK];
         Z >>= L;
      }
   }
   /* The last points */
   if (Last > 0) {
      Z = unif01_StripB (gen, r, s);
      /* The most significant bits make the points */
      for (j = 1; j <= t - Last; j++)
         Z >>= L;
      for (j = 1; j <= Last; j++) {
         ++res->Count[Z & MASK];
         Z >>= L;
      }
   }
}


/*=======================================================================*/

static void Multinom (unif01_Gen * gen, smultin_Param * par,
   smultin_Res * res, long N, long n, int r, long d, int t, lebool Sparse,
   smultin_CellType k, char *TestName, chrono_Chrono * Timer, lebool BitFlag)
/* 
 * If BitFlag = TRUE, this procedure was called from smultin_MultinomialBits,
 * otherwise from smultin_Multinomial. 
 * In the case BitFlag = TRUE, t stand for s, d for L. Otherwise, all
 * parameters are as in smultin_Multinomial.
 *
 * Sparse:   normal approximation for Delta != -1.
 * Non sparse:  chi-square approximation.
 * Collisions test meaningfull only in Sparse case.
 */
{
   long Seq;                      /* Replication number */
   double NbExp;                  /* Expected number per cell */
   double EColl;                  /* Approx. expected number of collisions */
   long Hache;                    /* Hashing module */
   double UnSurHache;
   double HacheLR;                /* Dimension of hashing table */
   long i;
   long CoMax;                    /* Maximum number of balls in any cell */
   double X0, X;                  /* Statistics */
   DeltaIndex j;                  /* Which power divergence case */
   double SumX2[smultin_MAX_DELTA];
   double SumX[smultin_MAX_DELTA];
   double X0Pre[smultin_MAX_DELTA]; /* For empirical mean and correlation */
   lebool localRes = FALSE;

   NbExp = (double) n / k;
   EColl = (double) n / (2.0 * k) * n;

   if (par == NULL)
      par = &smultin_ParamDefault;
   if (res == NULL) {
      localRes = TRUE;
      res = smultin_CreateRes (par);
   } else
      /* Clean memory from a previous call */
      CleanPD (res);

   InitRes (par, res, N);
   res->NbCellsTotal = k;
   res->Over = FALSE;
   InitPowDiv (par, res, N, Sparse, n, k);

   if (swrite_Basic) {
      if (BitFlag)
         /* Here t stand for s, d for L */
         WriteDataMNBits (gen, par, res, TestName, N, n, r, d, t, Sparse, k,
                          FALSE);
      else
         WriteDataPowDiv (gen, par, res, TestName, N, n, r, d, t, Sparse, k);
   }
   /* Initialize the hashing constants and tables */
   CalcTabFj (par, res, Sparse, (double) k, NbExp);
   for (j = 0; j < par->NbDelta; j++) {
      SumX[j] = 0.0;
      SumX2[j] = 0.0;
      X0Pre[j] = 0.0;
   }
   if (res->Hashing)
      Hache = tables_HashPrime (n, smultin_env.HashLoad);
   else
      Hache = k;
   HacheLR = Hache;
   UnSurHache = 1.0 / HacheLR;
   res->CountSize = Hache;
   res->Count = util_Calloc ((size_t) Hache + 2, sizeof (long));
   res->Cell = util_Calloc ((size_t) Hache + 2, sizeof (smultin_CellType));
   res->NbSize = 8000;
   res->Nb = util_Calloc ((size_t) res->NbSize + 2, sizeof (smultin_CellType));

   /* Generate the points or balls */
   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 0; i <= Hache; i++)
         res->Count[i] = 0;
      res->Nb[0] = k;
      for (i = 1; i <= res->NbSize; i++)
         res->Nb[i] = 0;

      if (BitFlag) {
         /* Here, d stands for L, and t for s */
         if (res->Hashing)
            GenerAllPointsHashBits (gen, res, n, r, d, t, &CoMax, Hache,
               UnSurHache);
         else
            GenerAllPoints2Bits (gen, res, n, r, d, t);
      } else {
         if (res->Hashing)
            GenerAllPointsHash (gen, par, res, n, r, d, t, &CoMax, Hache,
               UnSurHache);
         else
            GenerAllPoints2 (gen, par, res, n, r, d, t);
      }

      if (swrite_Counters) {
         if (res->Hashing)
#ifdef USE_LONGLONG
            tables_WriteTabULL (res->Nb, 0, CoMax, 5, 12,
               "Observed numbers in res->Nb");
#else
            tables_WriteTabD (res->Nb, 0, CoMax, 5, 12, 0, 0,
               "Observed numbers in res->Nb");
#endif
         else if (!Sparse)
            tables_WriteTabL (res->Count, 0, res->CountSize - 1, 5, 10,
                              "Observed numbers in res->Count");
      }

      /* The points have been generated; now compute the statistics */
      /* if (par->bmax >= 0) */
      CalcNbCells (par, res, 0, Hache - 1, CoMax);

      for (j = 0; j < par->NbDelta; j++) {
         if (res->Hashing) {
            CalcPoDiEqHache (par, res, j, NbExp, res->Nb, CoMax, TRUE, &X);
         } else if (res->flagTab) {
            CalcPowDivEqual (par, res, j, NbExp, res->Count,
                             0, (long) k - 1, TRUE, &X);
         } else {
            CalcPowDivEqual (par, res, j, NbExp, res->Count,
                             0, (long) k - 1, FALSE, &X);
         }
         X0 = (X - res->Mu[j]) / res->Sigma[j];
         if (fabs (par->ValDelta[j] + 1.0) < EPS_LAM) {
            res->Nb[0] = k + X - n;
            res->NbCollisions += X;
            statcoll_AddObs (res->Collector[j], X);
         } else {
            statcoll_AddObs (res->Collector[j], X0);
            if (!Sparse)
               X0 = (X0 - k + 1.0) / sqrt (2.0 * k - 2.0);
         }
         /* Now, X0 is standardized, with mean 0 and variance 1.  */
         /* The following is to compute the mean and correlation. */
         SumX[j] += X0;
         SumX2[j] += X0 * X0Pre[j];
         X0Pre[j] = X0;
      }
   }

   for (j = 0; j < par->NbDelta; j++) {
      if ((Sparse == FALSE) && (fabs (par->ValDelta[j] + 1.0) < EPS_LAM))
         continue;
      CalcResultsPowDiv (par, res, j, n, Sparse, k - 1, res->Mu[j],
                         SumX, SumX2);
   }

   if (swrite_Basic) {
      for (j = 0; j < par->NbDelta; j++) {
         if ((Sparse == FALSE) && (fabs (par->ValDelta[j] + 1.0) < EPS_LAM)) {
            util_Warning (TRUE,
               "The collision test is meaningless when Sparse = FALSE");
            continue;
         }
         WriteResultsPowDiv (par, res, j, N, EColl, k - 1, Sparse, res->Mu[j]);
      }
      swrite_Final (gen, Timer);
   }
   if (localRes)
      smultin_DeleteRes (res);
}


/*=======================================================================*/

void smultin_Multinomial (unif01_Gen * gen, smultin_Param * par,
   smultin_Res * res, long N, long n, int r, long d, int t, lebool Sparse)
/* 
 * Sparse:   normal approximation for Delta != -1.
 * Non sparse:  chi-square approximation.
 * Collisions test meaningfull only in Sparse case.
 */
{
   smultin_CellType k;            /* Number of cells */
   int i;
   chrono_Chrono *Timer;
   char *TestName = "smultin_Multinomial test";

   Timer = chrono_Create ();
   if (NULL == par)
      par = &smultin_ParamDefault;

   if (par->GenerCell == smultin_GenerCellSerial ||
      par->GenerCell == smultin_GenerCellSerial2) {
      util_Assert (d > 1, "smultin_Multinomial:   d <= 1");
      util_Assert (t > 0, "smultin_Multinomial:   t < 1");
      k = d;
      for (i = 2; i <= t; i++)
         k *= d;

   } else if (par->GenerCell == smultin_GenerCellPermut) {
      util_Assert (t > 1, "Permutation... smultin_Multinomial:   t < 2");
#ifdef USE_LONGLONG
      /* longlong has more bits of precision than double */
      util_Assert (t <= 20, "smultin_GenerCellPermut:  t > 20");
      if (t == 20) {
         k = num2_Factorial (18) * 19 * 20;
      } else if (t == 19) {
         k = num2_Factorial (18) * 19;
      } else
         k = num2_Factorial (t);
#else
      util_Assert (t <= 18, "smultin_GenerCellPermut:  t > 18");
      k = num2_Factorial (t);
#endif

   } else if (par->GenerCell == smultin_GenerCellMax) {
      util_Assert (t > 1, "GenerCellMax... smultin_Multinomial:   t < 2");
      k = t;

   } else
      util_Error ("smultin_Multinomial:   par->GenerCell not initialized");

   util_Assert (k <= smultin_env.Maxk,
      "smultin_Multinomial:   k > smultin_env.Maxk");
   util_Assert (n > 4, "smultin_Multinomial:   n <= 4");
#ifndef USE_LONGLONG
   util_Assert ((double) n / k > 1.0 / num_TwoExp[31],
      "smultin_Multinomial:   NbExp <= 1/2^31");
#endif
   Multinom (gen, par, res, N, n, r, d, t, Sparse, k, TestName, Timer, FALSE);
   chrono_Delete (Timer);
}


/*=======================================================================*/

static void InitCollOver (
   smultin_Res *res, 
   long n,                    /* Number of balls */
   smultin_CellType k,        /* Number of cells = d^t */
   long d,                    /* One-dim. segment */
   int t,                     /* dimension */
   double *Esperance,         /* Expectation value */
   double *StandDev           /* Standard deviation */
   )
/*
 * Initialize the collisionOver test
 */
{
   const double Epsilon = 1.0E-20;
   const int MaxIter = 32;
   long j;
   double terme;
   double v;
   double COverDelta;

   res->NbCollisions = 0.0;
   res->NbCells[0] = 0;
   res->CollApprox = smultin_CollNotInit;
   *Esperance = -1.0;
   *StandDev = -1.0;
   COverDelta = (double) (n - t + 1)/ k;

   if (COverDelta > smultin_env.SeuilCOverNorSup &&
      COverDelta < smultin_env.SeuilCOverDense) {
      res->CollApprox = smultin_CollPoissonDense;
      *Esperance = k * exp (-COverDelta);

   } else if (COverDelta >= smultin_env.SeuilCOverNorInf &&
      COverDelta <= smultin_env.SeuilCOverNorSup) {
      res->CollApprox = smultin_CollNormal;
      *Esperance = k * (COverDelta - 1.0 + exp (-COverDelta));
      terme = k * exp (-COverDelta) * (1.0 - (1.0 + COverDelta) *
                        exp (-COverDelta));
      /* The general formula given by Marsaglia is not very good; the above
         formula used Rukhin's correction. The following values were obtained
         by Marsaglia by simulation for the special cases: */
      if (n == 2097152 && k == 1048576) {
         if ((d == 32) && (t == 4)) {
            terme = 295.0 * 295.0; /* OQSO test */
         }
         if ((d == 4) && (t == 10)) {
            terme = 339.0 * 339.0; /* DNA test */
         }
      }
      if (terme < 0.0) {
         util_Warning (TRUE, "***** InitCollOver ******* VARIANCE < 0 !!");
         *Esperance = -1.0;
         *StandDev = -1.0;
      } else
         *StandDev = sqrt (terme);

   } else if (COverDelta < smultin_env.SeuilCOverSparse) {
      res->CollApprox = smultin_CollPoissonSparse;
      if (COverDelta < 0.1) {
         /* Avoid loss of precision when COverDelta --> 0 */
         j = 3;
         v = 2.0;
         terme = COverDelta * COverDelta / 2.0;
         *Esperance = terme;
         while (fabs (terme / *Esperance) > Epsilon && j < MaxIter) {
            v += 1.0;
            terme = -terme * COverDelta / v;
            *Esperance += terme;
            ++j;
         }
         *Esperance *= k;
      } else
         *Esperance = k * (COverDelta - 1.0 + exp (-COverDelta));
   }
}


/*=======================================================================*/

static void WriteDataCollOver (
   smultin_Res *res, 
   long n,                    /* Number of balls */
   smultin_CellType k,        /* Number of cells */
   double Esperance,          /* Expectation value */
   double StandDev            /* Standard deviation */
   )
{
   double COverDelta = (double) (n) / k;
   printf ("       CollisionOver:   density = n / k = ");
   if (COverDelta >= 1.0)
      num_WriteD (COverDelta, 10, 2, 2);
   else {
      printf (" 1 / ");
      num_WriteD (1.0 / COverDelta, 10, 2, 1);
   }
   printf ("\n");
   if (res->CollApprox == smultin_CollPoissonDense) {
      printf ("       Expected number of empty cells = Mu = ");
      num_WriteD (Esperance, 10, 2, 2);
      printf ("\n");
   } else if (res->CollApprox == smultin_CollNormal) {
      printf ("       Expected number of collisions = ");
      num_WriteD (Esperance, 10, 2, 2);
      printf ("\n");
      printf ("       Expected standard deviation = ");
      num_WriteD (StandDev, 10, 2, 2);
   } else if (res->CollApprox == smultin_CollPoissonSparse) {
      printf ("       Expected number of collisions = Mu = ");
      num_WriteD (Esperance, 10, 2, 2);
   } else {
      printf ("       NO TEST FOR THIS DENSITY  n/k");
   }
   printf ("\n\n");
}


/*=======================================================================*/

static void CalcResCollOver (
   smultin_Res *res, 
   DeltaIndex s,
   long N,                    /* Number of replications */
   double Esperance,          /* Expectation value */
   double SumX,
   double SumX2
   )
/*
 * Compute results for CollisionOver test
 */
{
   double pCollLeft;              /* Left p-value of Collision test */
   double pCollRight;             /* Right p-value of Collision test */
   double racN;
   fmass_INFO W;
   statcoll_Collector *Q = res->Collector[s];

   res->Mu[s] = Esperance;
   if (Esperance < 0.0) {
      res->pVal2[s][gofw_KSP] = -1.0;
      res->pVal2[s][gofw_Mean] = -1.0;
      res->pColl = -1.0;
      return;
   }

   switch (res->CollApprox) {

   case smultin_CollNormal:
      gofw_ActiveTests1 (Q->V, Q->NObs, wdist_Normal, (double *) NULL,
         res->sVal2[s], res->pVal2[s]);
      /* This line is necessary for the array pd from module tmultin */
      res->pColl = res->pVal2[s][gofw_Mean];
      if (N > 1) {
         racN = sqrt ((double) N);
         /* Calculate the mean, the correlation and their p-values */
         res->sVal2[s][gofw_Mean] = SumX / racN;
         res->pVal2[s][gofw_Mean] = fbar_Normal1 (res->sVal2[s][gofw_Mean]);
         res->sVal2[s][gofw_Cor] = racN * SumX2 / (N - 1);
         res->pVal2[s][gofw_Cor] = fbar_Normal1 (res->sVal2[s][gofw_Cor]);
         /* This line is necessary for the array pd from module tmultin */
         res->pColl = res->pVal2[s][gofw_KSP];
      }
      break;

   case smultin_CollPoissonSparse:
      /* The sum of N Poisson obeys also a Poisson law */
      W = fmass_CreatePoisson (N * Esperance);
      pCollLeft = fdist_Poisson2 (W, (long) res->NbCollisions);
      pCollRight = fbar_Poisson2 (W, (long) res->NbCollisions);
      res->pColl = gofw_pDisc (pCollLeft, pCollRight);
      fmass_DeletePoisson (W);
      break;

   case smultin_CollPoissonDense:
      /* The sum of N Poisson obeys also a Poisson law */
#if 0
      W = fmass_CreatePoisson (N * Esperance);
      pCollLeft = fdist_Poisson2 (W, res->NbCells[0]);
      pCollRight = fbar_Poisson2 (W, res->NbCells[0]);
      res->pColl = res->pEmpty = gofw_pDisc (pCollLeft, pCollRight);
      fmass_DeletePoisson (W);
#endif
      /* Since we can have very large values of res->NbCells[0]) (2^52), we
         compute the Poisson pCollLeft, pCollRight by calling the Gamma
         distribution, since our Poisson takes a long (31 bits) argument. */

      if (res->NbCells[0] <= EMPTYLIM && N * Esperance <= EMPTYLIM) {
         pCollLeft = fbar_Gamma (res->NbCells[0] + 1.0, 15, N * Esperance);
         if ((res->NbCells[0] <= 0) || (res->NbCells[0] > res->NbCellsTotal))
            pCollRight = 1.0;
         else
            pCollRight = fdist_Gamma ((double) (res->NbCells[0]), 15,
               N * Esperance);
         res->pEmpty = res->pColl = gofw_pDisc (pCollLeft, pCollRight);
      }
      break;

   default:
      util_Error ("res->CollApprox:   Impossible case");
   }
}


/*=======================================================================*/

static void WriteResCollOver (
   smultin_Param *par,
   smultin_Res *res, 
   DeltaIndex s,
   long N,                    /* Number of replications */
   double EColl,
   double Esperance
   )
/*
 * Write results for CollisionOver test
 */
{
   int j;
   printf ("\n-----------------------------------------------\n"
           "Results of CollisionOver test:\n\n");
   if (Esperance < 0.0) {
      util_Warning (TRUE, "TEST NON IMPLEMENTED FOR THESE PARAMETERS");
      return;
   }

   switch (res->CollApprox) {

   case smultin_CollNormal:
      printf ("NORMAL approximation:\n");
      if (N == 1) {
         printf ("Value of the standardized statistic   :");
         gofw_Writep2 (res->sVal2[s][gofw_Mean], res->pVal2[s][gofw_Mean]);
      } else {
         gofw_WriteActiveTests0 (N, res->sVal2[s], res->pVal2[s]);
         printf ("Standardized empirical mean           :");
         gofw_Writep2 (res->sVal2[s][gofw_Mean], res->pVal2[s][gofw_Mean]);
         printf ("Standardized empirical correlation    :");
         gofw_Writep2 (res->sVal2[s][gofw_Cor], res->pVal2[s][gofw_Cor]);
      }
      break;

   case smultin_CollPoissonSparse:
      /* The sum of N Poisson random variables is also a Poisson r. v. */
      printf ("POISSON approximation                 :\n");
      /* "Value of N * EColl : "); num_WriteD (N * EColl, 11, 2, 2); */
      printf ("Expected number of collisions = N*Mu  : ");
      num_WriteD (N * Esperance, 11, 2, 2);
      printf ("\nObserved number of collisions         : %8ld\n",
         (long) res->NbCollisions);
      gofw_Writep1 (res->pColl);
      break;

   case smultin_CollPoissonDense:
      /* The sum of N Poisson random variables is also a Poisson r. v. */
      printf ("POISSON approximation                 :\n"
              "Expected number of empty cells = N*Mu : ");
      num_WriteD (N * Esperance, 18, 2, 2);
#ifdef USE_LONGLONG
      printf ("\nObserved number of empty cells        : %15" PRIuLEAST64 "\n",
              res->NbCells[0]);
#else
      printf ("\nObserved number of empty cells        : %15.0f\n",
              res->NbCells[0]);
#endif
      gofw_Writep1 (res->pColl);
      break;

   default:;
      util_Error ("smultin_WriteResCollOver:  IMPOSSIBLE CASE");
      break;
   }


   if (swrite_Collectors)
      statcoll_Write (res->Collector[s], 5, 14, 2, 1);

   printf ("-----------------------------\n"
      "Total number of cells containing j balls\n\n");
   for (j = 0; j <= smultin_MAXB / 2; j++) {
      printf ("  j = %2d", j);
#ifdef USE_LONGLONG
      printf ("                              : %16" PRIuLEAST64 "\n",
              res->NbCells[j]);
#else
      printf ("                              : %16.0f\n", res->NbCells[j]);
#endif
   }
   printf ("\n");
}


/*=======================================================================*/

static void OverDenseGenere (
   unif01_Gen *gen, 
   smultin_Res *res, 
   long n,                    /* Number of balls */
   int r,
   long d,                    /* Division of 1-dim interval */
   int t,                     /* Dimension */
   long k,                    /* Number of urns in t dimensions */
   long k1                    /* Number of urns in t - 1 dimensions */
   )
/*
 * Generate the n balls for smultin_MultinomialOver in the dense
 * case, and fill the counters Count and Count1
 */
{
   long element;
   long Indice;
   long j;
   long Premier[MAX_DIM];
   long *Count = res->Count;      /* Counters in t dimensions */
   long *Count1 = res->Count1;    /* Counters in t - 1 dimensions */
   smultin_CellType *Nb = res->Nb;

   util_Assert (t < MAX_DIM, "OverDenseGenere:   t > 64");
   for (j = 1; j <= res->NbSize; j++)
      Nb[j] = 0;
   Nb[0] = k;
   for (j = 0; j <= k; j++)
      Count[j] = 0;
   for (j = 0; j <= k1; j++)
      Count1[j] = 0;

   /* Generation of the first (t - 1) random numbers for the first tuple. */
   /* We shall keep them in the array Premier[] since the sequence of     */
   /* generated numbers must be circular. They will be used to build the  */
   /* last t - 1 tuples. Here, tuples are balls or points.                */
   Indice = 0;
   for (j = 1; j < t; j++) {
      element = unif01_StripL (gen, r, d);
      Premier[j] = element;
      /* Shift tuple by s bits and insert new element */
      Indice = Indice * d + element;
   }

   /* Generation of the first n - (t - 1) tuples */
   for (j = 1; j <= n - (t - 1); j++) {
      Indice %= k1;
      ++Count1[Indice];
      Indice = Indice * d + unif01_StripL (gen, r, d);
      ++Count[Indice];
   }

   /* Generation of the last (t - 1) tuples. Use the elements of Premier */
   for (j = 1; j < t; j++) {
      Indice %= k1;
      ++Count1[Indice];
      Indice = Indice * d + Premier[j];
      ++Count[Indice];
   }
}


/*=======================================================================*/

static void OverHashGenere (
   unif01_Gen *gen, 
   smultin_Res *res, 
   long n,                    /* Number of balls */
   int r,
   smultin_CellType dLR,      /* Parameter d */
   int t,                     /* Dimension */
   long Hache1,               /* Size of hashing table in t dimensions */
   long Hache11,              /* Size of hashing table in t - 1 dimensions */
   smultin_CellType k,        /* Number of urns in t dimensions */
   smultin_CellType k1,       /* Number of urns in t - 1 dimensions */
   long *CoMax,               /* Max number of balls in any cell in t dim. */ 
   long *CoMax1               /* Max number of balls in any cell in t-1 dim. */
   )
/*
 * Generate the n balls for smultin_MultinomialOver in the sparse
 * case, and fill the counters. We use hashing.
 */
{
   long j, tem;
   long d = dLR;
   smultin_CellType Indice;
   smultin_CellType element;
   double UnSurHache1;
   double UnSurHache11;
   double UnSurk1;
   smultin_CellType Premier[MAX_DIM];
   long *Count = res->Count;      /* Counters in t dimensions */
   long *Count1 = res->Count1;    /* Counters in t - 1 dimensions */
   smultin_CellType *Cell = res->Cell; /* Cell numbers in t dimensions */
   smultin_CellType *Cell1 = res->Cell1; /* Cell numbers in t - 1 dimensions */
   smultin_CellType *Nb = res->Nb;
   smultin_CellType *Nb1 = res->Nb1;

   util_Assert (t < MAX_DIM, "OverHashGenere:   t > 64");
   UnSurk1 = 1.0 / k1;
   UnSurHache1 = 1.0 / Hache1;
   UnSurHache11 = 1.0 / Hache11;

   for (j = 0; j <= Hache1; j++) {
      Count[j] = 0;
#ifdef USE_LONGLONG
      Cell[j] = MASK64;           /* Empty cells */
#else
      Cell[j] = -1.0;             /* Empty cells */
#endif
   }
   for (j = 0; j <= Hache11; j++) {
      Count1[j] = 0;
#ifdef USE_LONGLONG
      Cell1[j] = MASK64;          /* Empty cells */
#else
      Cell1[j] = -1.0;            /* Empty cells */
#endif
   }
   for (j = 1; j <= res->NbSize; j++)
      Nb[j] = 0;
   for (j = 1; j <= res->Nb1Size; j++)
      Nb1[j] = 0;
   Nb[0] = k;
   Nb1[0] = k1;
   *CoMax = 0;
   *CoMax1 = 0;

   /* Generation of the first (t - 1) elements of the first tuple. We shall
      keep them in array Premier[] since the sequence of generated numbers
      must be circular. They will be used to obtain the last t - 1 tuples.
      When we generate a random number, we keep s bits and they become the
      least significant element of the tuple. We then shift the elements so
      that the most significant element is dropped. The tuples are balls. */

   /* Generate the first (t - 1) components of the first tuple (ball) */
   Indice = 0;
   for (j = 1; j < t; j++) {
      element = unif01_StripL (gen, r, d);
      Premier[j] = element;
      /* Shift tuple by s bits and insert new element */
      Indice = Indice * dLR + element;
   }

   /* Generation of the first n - (t - 1) tuples */
   for (j = 1; j <= n - (t - 1); j++) {
      /* Operation % k1 */
#ifdef USE_LONGLONG
      Indice %= k1;
#else
      tem = Indice * UnSurk1;
      Indice -= k1 * tem;
#endif
      UpdateCountHash (res, Indice, Hache11, UnSurHache11, CoMax1, TRUE);
      Indice = Indice * dLR + unif01_StripL (gen, r, d);
      UpdateCountHash (res, Indice, Hache1, UnSurHache1, CoMax, FALSE);
   }

   /* Generate the last (t - 1) tuples. We use the elements of Premier[] */
   for (j = 1; j < t; j++) {
#ifdef USE_LONGLONG
      Indice %= k1;
#else
      tem = Indice * UnSurk1;
      Indice -= k1 * tem;
#endif
      UpdateCountHash (res, Indice, Hache11, UnSurHache11, CoMax1, TRUE);
      Indice = Indice * dLR + Premier[j];
      UpdateCountHash (res, Indice, Hache1, UnSurHache1, CoMax, FALSE);
   }
}


/*=======================================================================*/

static void OverDenseGenereBits (
   unif01_Gen *gen, 
   smultin_Res *res, 
   long n,                    /* Number of balls */
   int r,                     /* Drop first r bits of each random number */
   int L,                     /* Cells numbered with L bits */
   int s,                     /* Take s bits of each random number */
   long k,                    /* Number of urns in L dimensions */
   long k1                    /* Number of urns in L - 1 dimensions */
   )
/*
 * Generate the n balls for smultin_MultinomialBitsOver in the dense
 * case, and fill the counters Count and Count1
 */
{
   int j;
   long i;
   unsigned long Premier[MAX_DIM];
   long *Count = res->Count;      /* Counters in L dimensions */
   long *Count1 = res->Count1;    /* Counters in L - 1 dimensions */
   smultin_CellType *Nb = res->Nb;

   util_Assert (L < MAX_DIM, "OverDenseGenereBits:   L > 64");
   for (i = 1; i <= res->NbSize; i++)
      Nb[i] = 0;
   Nb[0] = k;
   for (i = 0; i <= k; i++)
      Count[i] = 0;
   for (i = 0; i <= k1; i++)
      Count1[i] = 0;

   if (L + s <= 32) {
      const unsigned long MASK = num_TwoExp[L] - 1.0;
      const unsigned long MASK1 = num_TwoExp[L - 1] - 1.0;
      const int t = (L - 1) / s + 1;
      unsigned long Z, Z0;
      int b;

      /* Generation of the first t*s random bits for the first tuple.  */
      /* We shall keep them in Premier since the sequence of generated */
      /* bits will be circular. */
      Z0 = 0;
      for (j = 0; j < t; j++) {
         Z0 <<= s;
         Premier[j] = unif01_StripB (gen, r, s);
         Z0 |= Premier[j];
      }

      /* Generation of other bits: main loop */
      for (i = 0; i < (n - t * s - 1) / s; i++) {
         Z = Z0 = (Z0 << s) | unif01_StripB (gen, r, s);
         for (j = 0; j < s; j++) {
            ++Count1[Z & MASK1];
            ++Count[Z & MASK];
            Z >>= 1;
         }
      }

      /* Generation of the last b random bits */
      Z0 = (Z0 << s) | unif01_StripB (gen, r, s);
      b = n % s;
      if (b)
         Z0 >>= (s - b);
      else
         b = s;

      Z = Z0;
      for (j = 0; j < b; j++) {
         ++Count1[Z & MASK1];
         ++Count[Z & MASK];
         Z >>= 1;
      }

      /* Must do last few bits using circular overlap with Premier */
      for (i = 0; i < t; i++) {
         Z = Z0 = (Z0 << s) | Premier[i];
         for (j = 0; j < s; j++) {
            ++Count1[Z & MASK1];
            ++Count[Z & MASK];
            Z >>= 1;
         }
      }
      return;
   }

#ifndef USE_LONGLONG
   util_Error ("OverDenseGenereBits:   L + s > 32");
#else

   /* ---------------------------------------------------------- */
   if (L + s <= 64) {
      const ulonglong MASK = num_TwoExp[L] - 1.0;
      const ulonglong MASK1 = num_TwoExp[L - 1] - 1.0;
      const int t = (L - 1) / s + 1;
      ulonglong Z, Z0;
      int b;

      /* Generation of the first t*s random bits */
      Z0 = 0;
      for (j = 0; j < t; j++) {
         Z0 <<= s;
         Premier[j] = unif01_StripB (gen, r, s);
         Z0 |= Premier[j];
      }

      /* Generation of the other random bits: main loop */
      for (i = 0; i < (n - t * s - 1) / s; i++) {
         Z = Z0 = (Z0 << s) | unif01_StripB (gen, r, s);
         for (j = 0; j < s; j++) {
            ++Count1[Z & MASK1];
            ++Count[Z & MASK];
            Z >>= 1;
         }
      }

      /* Generation of the last b random bits */
      Z0 = (Z0 << s) | unif01_StripB (gen, r, s);
      b = n % s;
      if (b)
         Z0 >>= s - b;
      else
         b = s;

      Z = Z0;
      for (j = 0; j < b; j++) {
         ++Count1[Z & MASK1];
         ++Count[Z & MASK];
         Z >>= 1;
      }

      /* Must do last few bits using circular overlap with Premier */
      for (i = 0; i < t; i++) {
         Z = Z0 = (Z0 << s) | Premier[i];
         for (j = 0; j < s; j++) {
            ++Count1[Z & MASK1];
            ++Count[Z & MASK];
            Z >>= 1;
         }
      }
      return;

      /* ---------------------------------------------------------- */
   } else {                       /* L + s > 64 */

      const ulonglong MASK = num_TwoExp[L] - 1.0;
      const ulonglong MASK1 = num_TwoExp[L - 1] - 1.0;
      const int t = (L - 1) / s + 1;
      const int q1 = 64 - L;
      const int q2 = s % q1;
      const int t2 = s / q1;
      ulonglong Z, Z0;
      unsigned long Bloc;
      int k, b;

      /* Generation of the first t*s random bits */
      Z0 = 0;
      for (j = 0; j < t; j++) {
         Z0 <<= s;
         Premier[j] = unif01_StripB (gen, r, s);
         Z0 |= Premier[j];
      }

      /* Generation of bits: main loop */
      for (i = 0; i < (n - t * s - 1) / s; i++) {
         Bloc = unif01_StripB (gen, r, s);

         /* Since L + s overflows a ulonglong, process a s-bit block in */
         /* t2 subblocks of q1 bits and one last subblock of q2 bits.   */
         for (k = 1; k <= t2; k++) {
            Z = Z0 = (Z0 << q1) | (Bloc >> (q2 + (t2 - k) * q1));
            for (j = 0; j < q1; j++) {
               ++Count1[Z & MASK1];
               ++Count[Z & MASK];
               Z >>= 1;
            }
         }
         Z = Z0 = (Z0 << q2) | Bloc;
         for (j = 0; j < q2; j++) {
            ++Count1[Z & MASK1];
            ++Count[Z & MASK];
            Z >>= 1;
         }
      }

      /* Generation of the last b random bits */
      Bloc = unif01_StripB (gen, r, s);
      b = n % s;
      if (0 == b)
         b = s;
      Bloc >>= s - b;
      {
         const int q3 = b % q1;
         const int t3 = b / q1;
         for (k = 1; k <= t3; k++) {
            Z = Z0 = (Z0 << q1) | (Bloc >> (q3 + (t3 - k) * q1));
            for (j = 0; j < q1; j++) {
               ++Count1[Z & MASK1];
               ++Count[Z & MASK];
               Z >>= 1;
            }
         }
         Z = Z0 = (Z0 << q3) | Bloc;
         for (j = 0; j < q3; j++) {
            ++Count1[Z & MASK1];
            ++Count[Z & MASK];
            Z >>= 1;
         }
      }

      /* Must do last few bits using circular overlap with Premier */
      for (i = 0; i < t; i++) {
         Bloc = Premier[i];
         for (k = 1; k <= t2; k++) {
            Z = Z0 = (Z0 << q1) | (Bloc >> (q2 + (t2 - k) * q1));
            for (j = 0; j < q1; j++) {
               ++Count1[Z & MASK1];
               ++Count[Z & MASK];
               Z >>= 1;
            }
         }
         Z = Z0 = (Z0 << q2) | Bloc;
         for (j = 0; j < q2; j++) {
            ++Count1[Z & MASK1];
            ++Count[Z & MASK];
            Z >>= 1;
         }
      }
      return;
   }
#endif
}


/*=======================================================================*/

static void OverHashGenereBits (
   unif01_Gen *gen, 
   smultin_Res *res, 
   long n,                    /* Number of balls */
   int r,
   int L,                     /* Dimension */
   int s,
   long Hache1,               /* Size of hashing table in t dimensions */
   long Hache11,              /* Size of hashing table in t - 1 dimensions */
   smultin_CellType k,        /* Number of urns in t dimensions */
   smultin_CellType k1,       /* Number of urns in t - 1 dimensions */
   long *CoMax,               /* Max number of balls in any cell in t dim. */ 
   long *CoMax1               /* Max number of balls in any cell in t-1 dim. */
   )
/*
 * Generate the n balls for smultin_MultinomialOver in the sparse
 * case, and fill the counters. We use hashing.
 */
{
   int j;
   long i;
   unsigned long Premier[MAX_DIM];
   smultin_CellType Indice;
   double UnSurHache1;
   double UnSurHache11;

   util_Assert (L < MAX_DIM, "OverHashGenereBits:   L > 64");
   UnSurHache1 = 1.0 / Hache1;
   UnSurHache11 = 1.0 / Hache11;

   for (j = 0; j <= Hache1; j++) {
      res->Count[j] = 0;
#ifdef USE_LONGLONG
      res->Cell[j] = MASK64;      /* Empty cells */
#else
      res->Cell[j] = -1.0;        /* Empty cells */
#endif
   }
   for (j = 0; j <= Hache11; j++) {
      res->Count1[j] = 0;
#ifdef USE_LONGLONG
      res->Cell1[j] = MASK64;     /* Empty cells */
#else
      res->Cell1[j] = -1.0;       /* Empty cells */
#endif
   }
   for (j = 1; j <= res->NbSize; j++)
      res->Nb[j] = 0;
   for (j = 1; j <= res->Nb1Size; j++)
      res->Nb1[j] = 0;
   res->Nb[0] = k;
   res->Nb1[0] = k1;
   *CoMax = 0;
   *CoMax1 = 0;

   if (L + s <= 32) {
      const unsigned long MASK = num_TwoExp[L] - 1.0;
      const unsigned long MASK1 = num_TwoExp[L - 1] - 1.0;
      const int t = (L - 1) / s + 1;
      unsigned long Z, Z0, b;

      /* Generation of the first t*s random bits for the first tuple. */
      /* We shall keep them in Premier since the sequence of */
      /* generated bits will be circular. */
      Z0 = 0;
      for (j = 0; j < t; j++) {
         Z0 <<= s;
         Premier[j] = unif01_StripB (gen, r, s);
         Z0 |= Premier[j];
      }

      /* Generation of all other bits: main loop */
      for (i = 0; i < (n - t * s - 1) / s; i++) {
         Z = Z0 = (Z0 << s) | unif01_StripB (gen, r, s);
         for (j = 0; j < s; j++) {
            Indice = Z & MASK1;
            UpdateCountHash (res, Indice, Hache11, UnSurHache11, CoMax1, TRUE);
            Indice = Z & MASK;
            UpdateCountHash (res, Indice, Hache1, UnSurHache1, CoMax, FALSE);
            Z >>= 1;
         }
      }

      /* Generation of the last b random bits */
      Z0 = (Z0 << s) | unif01_StripB (gen, r, s);
      b = n % s;
      if (b)
         Z0 >>= (s - b);
      else
         b = s;

      Z = Z0;
      for (j = 0; j < (int) b; j++) {
         Indice = Z & MASK1;
         UpdateCountHash (res, Indice, Hache11, UnSurHache11, CoMax1, TRUE);
         Indice = Z & MASK;
         UpdateCountHash (res, Indice, Hache1, UnSurHache1, CoMax, FALSE);
         Z >>= 1;
      }

      /* Must do last few bits using circular overlap with Premier */
      for (i = 0; i < t; i++) {
         Z = Z0 = (Z0 << s) | Premier[i];
         for (j = 0; j < s; j++) {
            Indice = Z & MASK1;
            UpdateCountHash (res, Indice, Hache11, UnSurHache11, CoMax1, TRUE);
            Indice = Z & MASK;
            UpdateCountHash (res, Indice, Hache1, UnSurHache1, CoMax, FALSE);
            Z >>= 1;
         }
      }
      return;
   }

   /* ---------------------------------------------------------- */
#ifndef USE_LONGLONG
   util_Error ("OverHashGenereBits:   L + s > 32");
#else

   if (L + s <= 64) {
      const ulonglong MASK = num_TwoExp[L] - 1.0;
      const ulonglong MASK1 = num_TwoExp[L - 1] - 1.0;
      const int t = (L - 1) / s + 1;
      ulonglong Z, Z0, b;

      /* Generation of the first t*s random bits */
      Z0 = 0;
      for (j = 0; j < t; j++) {
         Z0 <<= s;
         Premier[j] = unif01_StripB (gen, r, s);
         Z0 |= Premier[j];
      }

      /* Generation of the other random bits: main loop */
      for (i = 0; i < (n - t * s - 1) / s; i++) {
         Z = Z0 = (Z0 << s) | unif01_StripB (gen, r, s);
         for (j = 0; j < s; j++) {
            Indice = Z & MASK1;
            UpdateCountHash (res, Indice, Hache11, UnSurHache11, CoMax1, TRUE);
            Indice = Z & MASK;
            UpdateCountHash (res, Indice, Hache1, UnSurHache1, CoMax, FALSE);
            Z >>= 1;
         }
      }

      /* Generation of the last b random bits */
      Z0 = (Z0 << s) | unif01_StripB (gen, r, s);
      b = n % s;
      if (b)
         Z0 >>= (s - b);
      else
         b = s;

      Z = Z0;
      for (j = 0; j < (int) b; j++) {
         Indice = Z & MASK1;
         UpdateCountHash (res, Indice, Hache11, UnSurHache11, CoMax1, TRUE);
         Indice = Z & MASK;
         UpdateCountHash (res, Indice, Hache1, UnSurHache1, CoMax, FALSE);
         Z >>= 1;
      }

      /* Must do last few bits using circular overlap with Premier */
      for (i = 0; i < t; i++) {
         Z = Z0 = (Z0 << s) | Premier[i];
         for (j = 0; j < s; j++) {
            Indice = Z & MASK1;
            UpdateCountHash (res, Indice, Hache11, UnSurHache11, CoMax1, TRUE);
            Indice = Z & MASK;
            UpdateCountHash (res, Indice, Hache1, UnSurHache1, CoMax, FALSE);
            Z >>= 1;
         }
      }
      return;

      /* ---------------------------------------------------------- */
   } else {                       /* L + s > 64 */

      const ulonglong MASK = num_TwoExp[L] - 1.0;
      const ulonglong MASK1 = num_TwoExp[L - 1] - 1.0;
      const int t = (L - 1) / s + 1;
      const int q1 = 64 - L;
      const int t2 = s / q1;
      const int q2 = s % q1;
      ulonglong Z, Z0;
      unsigned long Bloc;
      int k, b;

      /* Generation of the first t*s random bits */
      Z0 = 0;
      for (j = 0; j < t; j++) {
         Z0 <<= s;
         Premier[j] = unif01_StripB (gen, r, s);
         Z0 |= Premier[j];
      }

      /* Generation of the other random bits: main loop */
      for (i = 0; i < (n - t * s - 1) / s; i++) {
         Bloc = unif01_StripB (gen, r, s);

         /* Since L + s overflows a ulonglong, process a s-bit block in */
         /* t2 subblocks of q1 bits and one last subblock of q2 bits.  */
         for (k = 1; k <= t2; k++) {
            Z = Z0 = (Z0 << q1) | (Bloc >> (q2 + (t2 - k) * q1));
            for (j = 0; j < q1; j++) {
               Indice = Z & MASK1;
               UpdateCountHash (res, Indice, Hache11, UnSurHache11,
                  CoMax1, TRUE);
               Indice = Z & MASK;
               UpdateCountHash (res, Indice, Hache1, UnSurHache1, CoMax,
                  FALSE);
               Z >>= 1;
            }
         }
         Z = Z0 = (Z0 << q2) | Bloc;
         for (j = 0; j < q2; j++) {
            Indice = Z & MASK1;
            UpdateCountHash (res, Indice, Hache11, UnSurHache11, CoMax1, TRUE);
            Indice = Z & MASK;
            UpdateCountHash (res, Indice, Hache1, UnSurHache1, CoMax, FALSE);
            Z >>= 1;
         }
      }

      /* Generation of the last b random bits */
      b = n % s;
      Bloc = unif01_StripB (gen, r, s);
      if (0 == b)
         b = s;
      Bloc >>= s - b;

      {
         const int q3 = b % q1;
         const int t3 = b / q1;
         for (k = 1; k <= t3; k++) {
            Z = Z0 = (Z0 << q1) | (Bloc >> (q3 + (t3 - k) * q1));
            for (j = 0; j < q1; j++) {
               Indice = Z & MASK1;
               UpdateCountHash (res, Indice, Hache11, UnSurHache11,
                  CoMax1, TRUE);
               Indice = Z & MASK;
               UpdateCountHash (res, Indice, Hache1, UnSurHache1,
                  CoMax, FALSE);
               Z >>= 1;
            }
         }
         Z = Z0 = (Z0 << q3) | Bloc;
         for (j = 0; j < q3; j++) {
            Indice = Z & MASK1;
            UpdateCountHash (res, Indice, Hache11, UnSurHache11, CoMax1, TRUE);
            Indice = Z & MASK;
            UpdateCountHash (res, Indice, Hache1, UnSurHache1, CoMax, FALSE);
            Z >>= 1;
         }
      }

      /* Must do last few bits using circular overlap with Premier */
      for (i = 0; i < t; i++) {
         Bloc = Premier[i];
         for (k = 1; k <= t2; k++) {
            Z = Z0 = (Z0 << q1) | (Bloc >> (q2 + (t2 - k) * q1));
            for (j = 0; j < q1; j++) {
               Indice = Z & MASK1;
               UpdateCountHash (res, Indice, Hache11, UnSurHache11,
                  CoMax1, TRUE);
               Indice = Z & MASK;
               UpdateCountHash (res, Indice, Hache1, UnSurHache1,
                  CoMax, FALSE);
               Z >>= 1;
            }
         }
         Z = Z0 = (Z0 << q2) | Bloc;
         for (j = 0; j < q2; j++) {
            Indice = Z & MASK1;
            UpdateCountHash (res, Indice, Hache11, UnSurHache11, CoMax1, TRUE);
            Indice = Z & MASK;
            UpdateCountHash (res, Indice, Hache1, UnSurHache1, CoMax, FALSE);
            Z >>= 1;
         }
      }
      return;
   }
#endif
}


/*=======================================================================*/

static void MultinomOver (unif01_Gen * gen, smultin_Param * par,
   smultin_Res * res, long N, long n, int r, long d, int t, lebool Sparse,
   smultin_CellType k, smultin_CellType k1, char *TestName,
   chrono_Chrono *Timer, lebool BitFlag)
{
   long Seq;
   smultin_CellType dLR = d;
   double nLR = n;
   double NbExp;              /* Expected number per cell in t dimensions */
   double NbExp1;             /* Expected number per cell in t - 1 dim. */
   double EColl;              /* Approx. expected number of collisions */
   DeltaIndex s;
   long Hache1, Hache11;      /* Hashing modules */
   long CoMax1;               /* Max number of balls in any cell: t-1 dim. */
   long CoMax;                /* Max number of balls in any cell: t dim. */
   double X, X0, X1;          /* Statistics */
   double Esperance;          /* Expected value of number of collisions */
   double StandDev;           /* Standard deviation of number of collisions */
   double V[1];               /* Number of degrees of freedom for ChiSquare */
   double SumX2[smultin_MAX_DELTA];
   double SumX[smultin_MAX_DELTA];
   double X0Pre[smultin_MAX_DELTA];
   lebool localRes = FALSE;

   NbExp = (double) n / k;
   NbExp1 = (double) n / k1;
   EColl = nLR * nLR / (2.0 * k);
   if (par == NULL)
      par = &smultin_ParamDefault;
   if (res == NULL) {
      localRes = TRUE;
      res = smultin_CreateRes (par);
   } else
      /* Clean memory from a previous call */
      CleanPD (res);

   res->NbCellsTotal = k;
   res->Over = TRUE;
   InitRes (par, res, N);
   InitPowDiv (par, res, N, Sparse, n, k - k1);
   if (swrite_Basic) {
      if (BitFlag)
         /* Here t stand for s, d for L */
         WriteDataMNBits (gen, par, res, TestName, N, n, r, d, t, Sparse, k,
                          TRUE);
      else
         WriteDataPowDiv (gen, par, res, TestName, N, n, r, d, t, Sparse, k);
   }
   for (s = 0; s < par->NbDelta; s++) {
      if (fabs (par->ValDelta[s] + 1.0) < EPS_LAM) {
         /* CollisionOver test */
         InitCollOver (res, n, k, d, t, &Esperance, &StandDev);
         if (swrite_Basic)
            WriteDataCollOver (res, n, k, Esperance, StandDev);
      }
   }
   for (s = 0; s < par->NbDelta; s++) {
      SumX[s] = 0.0;
      SumX2[s] = 0.0;
      X0Pre[s] = 0.0;
   }
   CalcTabFj (par, res, Sparse, (double) k, NbExp);
   if (res->Hashing) {
      Hache1 = tables_HashPrime (n, smultin_env.HashLoad);
      if ((unsigned) Hache1 > k1)
         Hache11 = k1;
      else
         Hache11 = Hache1;
      res->Cell = util_Calloc ((size_t) Hache1 + 2, sizeof (smultin_CellType));
      res->Cell1 = util_Calloc ((size_t) Hache11 + 2,
         sizeof (smultin_CellType));
   } else {
      Hache1 = k;
      Hache11 = k1;
   }
   res->CountSize = Hache1;
   res->Count1Size = Hache11;
   res->Count = util_Calloc ((size_t) Hache1 + 2, sizeof (long));
   res->Count1 = util_Calloc ((size_t) Hache11 + 2, sizeof (long));
   res->NbSize = res->Nb1Size = 8000;
   res->Nb = util_Calloc ((size_t) res->NbSize + 2, sizeof (smultin_CellType));
   res->Nb1 = util_Calloc ((size_t) res->Nb1Size + 2,
      sizeof (smultin_CellType));

   /* Generate the points or balls */
   for (Seq = 1; Seq <= N; Seq++) {
      if (BitFlag) {
         /* Here, d stands for L, and t for s */
         if (res->Hashing) {
            OverHashGenereBits (gen, res, n, r, d, t, Hache1, Hache11, k, k1,
               &CoMax, &CoMax1);
         } else {
            OverDenseGenereBits (gen, res, n, r, d, t, Hache1, Hache11);
         }
      } else {
         if (res->Hashing) {
            OverHashGenere (gen, res, n, r, dLR, t, Hache1, Hache11, k, k1,
               &CoMax, &CoMax1);
         } else {
            OverDenseGenere (gen, res, n, r, d, t, Hache1, Hache11);
         }
      }

      if (swrite_Counters) {
         if (res->Hashing) {
#ifdef USE_LONGLONG
            tables_WriteTabULL (res->Nb, 0, CoMax, 5, 12,
               "Observed numbers in res->Nb");
            tables_WriteTabULL (res->Nb1, 0, CoMax1, 5, 12,
               "Observed numbers in res->Nb1");
#else
            tables_WriteTabD (res->Nb, 0, CoMax, 5, 12, 0, 0,
               "Observed numbers in res->Nb");
            tables_WriteTabD (res->Nb1, 0, CoMax1, 5, 12, 0, 0,
               "Observed numbers in res->Nb1");
#endif
         } else if (!Sparse) {
            tables_WriteTabL (res->Count, 0, res->CountSize - 1, 5,
               10, "Observed numbers in res->Count");
            tables_WriteTabL (res->Count1, 0, res->Count1Size - 1, 5,
               10, "Observed numbers in res->Count1");
         }
      }

      /* The balls have been generated; now compute the statistics */
      for (s = 0; s < par->NbDelta; s++) {
         /* Compute the stat. X */
         if (res->Hashing) {
            CalcPoDiEqHache (par, res, s, NbExp, res->Nb, CoMax, TRUE, &X);

         } else if (res->flagTab) {
            CalcPowDivEqual (par, res, s, NbExp,
               res->Count, 0, Hache1 - 1, TRUE, &X);

         } else {
            CalcPowDivEqual (par, res, s, NbExp,
               res->Count, 0, Hache1 - 1, FALSE, &X);
         }

         if (fabs (par->ValDelta[s] + 1.0) < EPS_LAM && Esperance >= 0.0) {
            /* CollisionOver test */
            switch ((unsigned) res->CollApprox) {
            case smultin_CollPoissonDense: /* Number of empty cells */
               X0 = k + X - nLR;
               break;
            case smultin_CollPoissonSparse: /* Number of collisions */
               X0 = X;
               break;
            case smultin_CollNormal: /* Standardized number of collisions */
               X0 = (X - Esperance) / StandDev;
               break;
            default:
               util_Error
                  ("smultin_MultinomialOver: Computing X0 with CollNotInit");
               break;
            }
            res->NbCollisions += X;
            res->Nb[0] = k + X - nLR;
            statcoll_AddObs (res->Collector[s], X0);
            CalcNbCells (par, res, 0, Hache1 - 1, CoMax);

         } else {
            /* In the case delta = 1, X-X1 is approx. a chi-square with
               k - k1 degrees of freedom, or a normal in the sparse case */
            /* Compute X1 */
            if (res->Hashing) {
               CalcPoDiEqHache (par, res, s, NbExp1, res->Nb1,
                  CoMax1, FALSE, &X1);
            } else {
               CalcPowDivEqual (par, res, s, NbExp1,
                  res->Count1, 0, Hache11 - 1, FALSE, &X1);
            }
            X0 = (X - X1 - res->Mu[s]) / res->Sigma[s];
            statcoll_AddObs (res->Collector[s], X0);
            if (!Sparse)
               X0 = (X0 - k + k1) / sqrt (2.0 * (k - k1));
            /* Now, X0 is standardized, with mean 0 and variance 1.  */
         }

         /* The following is to compute the mean and correlation */
         SumX[s] += X0;
         SumX2[s] += X0 * X0Pre[s];
         X0Pre[s] = X0;
      }
   }

   /* For now, we understand only the cases delta = 1 and Collision */
   for (s = 0; s < par->NbDelta; s++) {
      statcoll_Collector *Q = res->Collector[s];
      double racN = sqrt ((double) N);

      if (par->ValDelta[s] > -1.0 + EPS_LAM) {
         /* Not Collisions test */
         if (Sparse) {
            util_Warning (fabs (par->ValDelta[s] - 1.0) > EPS_LAM,
  "The theoretical distribution for the overlapping case\nis known only for Delta = 1");
            gofw_ActiveTests1 (Q->V, Q->NObs, wdist_Normal,
               (double *) NULL, res->sVal2[s], res->pVal2[s]);
         } else {
            V[0] = k - k1;
            gofw_ActiveTests1 (Q->V, Q->NObs, wdist_ChiSquare, V,
               res->sVal2[s], res->pVal2[s]);
         }
         /* Compute the mean, the correlation, and their p-values */
         if (Q->NObs > 1) {
            res->sVal2[s][gofw_Mean] = SumX[s] / racN;
            res->pVal2[s][gofw_Mean] = fbar_Normal1 (res->sVal2[s][gofw_Mean]);
            res->sVal2[s][gofw_Cor] = racN * SumX2[s] / (N - 1);
            res->pVal2[s][gofw_Cor] = fbar_Normal1 (res->sVal2[s][gofw_Cor]);
         }
         if (swrite_Basic) {
            WriteResultsPowDiv (par, res, s, N, EColl, k - k1, Sparse,
               res->Mu[s]);
         }

      } else if (fabs (par->ValDelta[s] + 1.0) < EPS_LAM) {
         /* Collisions test */
         CalcResCollOver (res, s, N, Esperance, SumX[s], SumX2[s]);
         if (swrite_Basic) {
            WriteResCollOver (par, res, s, N, EColl, Esperance);
         }
      }
   }
   if (swrite_Basic)
      swrite_Final (gen, Timer);

   if (localRes)
      smultin_DeleteRes (res);
}


/*=======================================================================*/

void smultin_MultinomialOver (unif01_Gen * gen, smultin_Param * par,
   smultin_Res * res, long N, long n, int r, long d, int t, lebool Sparse)
{
   int i;
   smultin_CellType k1;       /* Number of urns in t - 1 dimensions */
   smultin_CellType k;        /* Number of urns in t dimensions */
   double NbExp;              /* Expected number per cell in t dimensions */
   chrono_Chrono *Timer;
   char *TestName = "smultin_MultinomialOver test";

   Timer = chrono_Create ();
   if (NULL == par)
      par = &smultin_ParamDefault;
   k1 = 1;
   for (i = 1; i < t; i++)
      k1 *= d;
   k = k1 * d;
   NbExp = (double) n / k;
   util_Assert (n > 4, "smultin_MultinomialOver:   n <= 4");
   util_Assert (t > 1, "smultin_MultinomialOver:   t < 2");
   if (par->GenerCell != smultin_GenerCellPermut)
      util_Assert (d > 1, "smultin_MultinomialOver:   d <= 1");
   util_Assert (k <= smultin_env.Maxk,
      "smultin_MultinomialOver:   d^t > Maxk");
#ifndef USE_LONGLONG
   util_Assert (NbExp > 1.0 / num_TwoExp[31],
      "smultin_MultinomialOver:   NbExp <= 1/2^31");
#endif
   MultinomOver (gen, par, res, N, n, r, d, t, Sparse, k, k1,
                 TestName, Timer, FALSE);
   chrono_Delete (Timer);
}


/*=======================================================================*/

void smultin_MultinomialBits (unif01_Gen *gen, smultin_Param *par,
   smultin_Res *res, long N, long n, int r, int s, int L, lebool Sparse)
{
/* 
 * Sparse:   normal approximation for Delta != -1.
 * Non sparse:  chi-square approximation.
 * Collisions test meaningfull only in Sparse case.
 */
   smultin_CellType k;            /* Number of cells */
   chrono_Chrono *Timer;
   char *TestName = "smultin_MultinomialBits test";

   Timer = chrono_Create ();
   k = num_TwoExp[L];
   if (NULL == par)
      par = &smultin_ParamDefault;
   if (L >= s) {
      long d = num_TwoExp[s];
      int t = L / s;
      if (swrite_Basic) {
         printf
            ("***********************************************************\n"
            "Test smultin_MultinomialBits calling smultin_Multinomial\n\n");
         printf ("   N = %2ld,  n = %2ld,  r = %1d", N, n, r);
         printf (",   s = %2d,   L = %2d,   Sparse = ", s, L);
         util_WriteBool (Sparse, 5);
         printf ("\n\n   Number of bits = n*L = %.0f\n\n\n", (double) n * L);
      }
      if ((t == 1) && (s > 30)) {
         util_Warning (TRUE, "smultin_MultinomialBits:   L = s  and  s > 30");
         return;
      }
      util_Assert (L % s == 0, "smultin_MultinomialBits:   L Mod s > 0");
      par->GenerCell = smultin_GenerCellSerial;
      smultin_Multinomial (gen, par, res, N, n, r, d, t, Sparse);
      return;
   }

   util_Assert (s % L == 0, "smultin_MultinomialBits:   s Mod L > 0");
   util_Assert (k <= smultin_env.Maxk,
      "smultin_MultinomialBits:   k > Maxk");
   util_Assert (n > 4, "smultin_MultinomialBits:   n <= 4");
#ifndef USE_LONGLONG
   util_Assert ((double) n / k > 1.0 / num_TwoExp[31],
      "smultin_MultinomialBits:   NbExp <= 1/2^31");
#endif
   Multinom (gen, par, res, N, n, r, L, s, Sparse, k, TestName, Timer, TRUE);
   chrono_Delete (Timer);
}


/*=======================================================================*/

void smultin_MultinomialBitsOver (unif01_Gen * gen, smultin_Param * par,
   smultin_Res * res, long N, long n, int r, int s, int L, lebool Sparse)
{
   smultin_CellType k1;       /* Number of urns in L - 1 dimensions */
   smultin_CellType k;        /* Number of urns in L dimensions */
   double NbExp;              /* Expected number per cell in L dimensions */
   chrono_Chrono *Timer;
   char *TestName = "smultin_MultinomialBitsOver test";

   Timer = chrono_Create ();
   if (NULL == par)
      par = &smultin_ParamDefault;
   util_Assert (L <= 64, "smultin_MultinomialBitsOver:   L > 64");
   k1 = num_TwoExp[L - 1];
   k = num_TwoExp[L];
   NbExp = (double) n / k;
   util_Assert (n > 4, "smultin_MultinomialBitsOver:   n <= 4");
   util_Assert (L > 1, "smultin_MultinomialBitsOver:   L < 2");
   util_Assert (s > 0, "smultin_MultinomialBitsOver:   s < 1");
   util_Assert (k <= smultin_env.Maxk,
      "smultin_MultinomialBitsOver:   L too large");
#ifndef USE_LONGLONG
   util_Assert (NbExp > 1.0 / num_TwoExp[31],
      "smultin_MultinomialBitsOver:   NbExp <= 1/2^31");
#endif
   MultinomOver (gen, par, res, N, n, r, L, s, Sparse, k, k1,
                 TestName, Timer, TRUE);

   chrono_Delete (Timer);
}
