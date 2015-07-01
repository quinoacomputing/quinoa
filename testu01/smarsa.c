/*************************************************************************\
 *
 * Package:        TestU01
 * File:           smarsa.c
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

#include "smarsa.h"
#include "smultin.h"
#include "wdist.h"
#include "swrite.h"
#include "unif01.h"

#include "vectorsF2.h"

#include "gofs.h"
#include "gofw.h"
#include "fdist.h"
#include "fbar.h"
#include "fmass.h"
#include "statcoll.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define LENGTH 200

/* MAXK = 2^64 */
#define STR_MAXK "18446744073709551616"



/*---------------------------- Extern variables ---------------------------*/

#ifdef USE_LONGLONG
double smarsa_Maxk = 18446744073709551616.0;   /* 2^64 */ 
#else
double smarsa_Maxk = num_MaxIntDouble;        /* 2^53 */
#endif




/*------------------------------- Functions -------------------------------*/

static void WriteResultsPoisson (sres_Poisson *res, long N)
{
   printf ("\n----------------------------------------------------"
           "\nTotal expected number = N*Lambda      : ");
   num_WriteD (N * res->Lambda, 10, 2, 2);
   printf ("\nTotal observed number                 : %7ld\n",
      (long) res->sVal2);
   gofw_Writep1 (res->pVal2);
   printf ("\n");
}


/*=========================================================================*/

static void InitRes (
   smarsa_Res *res,            /* Results holder */
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

smarsa_Res * smarsa_CreateRes (void)
{
   smarsa_Res *res;
   res = util_Malloc (sizeof (smarsa_Res));
   res->Bas = sres_CreateBasic ();
   res->Pois = sres_CreatePoisson ();
   res->Pois->pLeft = -1.0;
   res->Pois->pRight = -1.0;
   return res;
}


/*-------------------------------------------------------------------------*/

void smarsa_DeleteRes (smarsa_Res *res)
{
   if (res == NULL)
      return;
   sres_DeleteBasic (res->Bas);
   sres_DeletePoisson (res->Pois);
   util_Free (res);
}


/*=========================================================================*/

static void InitRes2 (
   smarsa_Res2 *res,          /* Results holder */
   long N,                    /* Number of replications */
   int jmax,                  /* Max class index for GCD */
   int tmax                   /* Max class index for NumIter */
)
/* 
 * Initializes the smarsa_Res2 structure
 */
{
   sres_InitChi2 (res->GCD, N, jmax, "smarsa_GCD:   GCD");
   sres_InitChi2 (res->NumIter, N, tmax, "smarsa_GCD:   NumIter");
}


/*-------------------------------------------------------------------------*/

smarsa_Res2 *smarsa_CreateRes2 (void)
{
   smarsa_Res2 *res;
   res = util_Malloc (sizeof (smarsa_Res2));
   res->GCD = sres_CreateChi2 ();
   res->NumIter = sres_CreateChi2 ();
   return res;
}


/*-------------------------------------------------------------------------*/

void smarsa_DeleteRes2 (smarsa_Res2 *res)
{
   if (res == NULL)
      return;
   sres_DeleteChi2 (res->GCD);
   sres_DeleteChi2 (res->NumIter);
   util_Free (res);
}


/*=========================================================================*/

void smarsa_SerialOver (unif01_Gen *gen, sres_Basic *res,
   long N, long n, int r, long d, int t)
{
   double ValDelta[] = { 1.0 };
   smultin_Param *par;

   if (swrite_Basic)
      printf ("***********************************************************\n"
         "Test smarsa_SerialOver calling smultin_MultinomialOver\n\n");

   par = smultin_CreateParam (1, ValDelta, smultin_GenerCellSerial, 0);
   if (NULL == res) {
      smultin_MultinomialOver (gen, par, NULL, N, n, r, d, t, FALSE);
   } else {
      smultin_Res *resm;
      resm = smultin_CreateRes (par);
      smultin_MultinomialOver (gen, par, resm, N, n, r, d, t, FALSE);
      sres_InitBasic (res, N, "smarsa_SerialOver");
      statcoll_SetDesc (res->sVal1, "SerialOver sVal1");
      res->sVal1->NObs = resm->Collector[0]->NObs;
      tables_CopyTabD (resm->Collector[0]->V, res->sVal1->V, 1, N);
      tables_CopyTabD (resm->sVal2[0], res->sVal2, 0, gofw_NTestTypes - 1);
      tables_CopyTabD (resm->pVal2[0], res->pVal2, 0, gofw_NTestTypes - 1);
      smultin_DeleteRes (resm);
   }
   smultin_DeleteParam (par);
}


/*=========================================================================*/

void smarsa_CollisionOver (unif01_Gen *gen, smarsa_Res *res,
   long N, long n, int r, long d, int t)
{
   double ValDelta[] = { -1.0 };
   smultin_Param *par;

   if (swrite_Basic)
      printf ("***********************************************************\n"
         "Test smarsa_CollisionOver calling smultin_MultinomialOver\n\n");

   par = smultin_CreateParam (1, ValDelta, smultin_GenerCellSerial, 3);
   if (NULL == res) {
      smultin_MultinomialOver (gen, par, NULL, N, n, r, d, t, TRUE);
   } else {
      smultin_Res *resm;
      resm = smultin_CreateRes (par);
      smultin_MultinomialOver (gen, par, resm, N, n, r, d, t, TRUE);
      InitRes (res, N, resm->Mu[0], "smarsa_CollisionOver");
      statcoll_SetDesc (res->Bas->sVal1, "CollisionOver sVal1");
      statcoll_SetDesc (res->Pois->sVal1, "CollisionOver sVal1");
      res->Pois->sVal1->NObs = resm->Collector[0]->NObs;
      res->Bas->sVal1->NObs = resm->Collector[0]->NObs;
      tables_CopyTabD (resm->Collector[0]->V, res->Bas->sVal1->V, 1, N);
      tables_CopyTabD (resm->Collector[0]->V, res->Pois->sVal1->V, 1, N);
      res->Pois->pVal2 = resm->pColl;
      if (resm->CollApprox == smultin_CollPoissonSparse)
         res->Pois->sVal2 = resm->NbCollisions;
      else
         res->Pois->sVal2 = resm->NbCells[0];
      tables_CopyTabD (resm->sVal2[0], res->Bas->sVal2, 0,
         gofw_NTestTypes - 1);
      tables_CopyTabD (resm->pVal2[0], res->Bas->pVal2, 0,
         gofw_NTestTypes - 1);
      smultin_DeleteRes (resm);
   }
   smultin_DeleteParam (par);
}


/*=========================================================================*/

void smarsa_Opso (unif01_Gen * gen, smarsa_Res * res, long N, int r, int p)
{
   int d;
   long NBalls;

   switch (p) {
   case 1:
      NBalls = 2097152;
      d = 1024;
      break;
   case 2:
      NBalls = 4194304;
      d = 2048;
      break;
   case 3:
      NBalls = 8388608;
      d = 2048;
      break;
   default:
      util_Error ("smarsa_Opso:  p must be in {1, 2, 3}");
   }

   if (swrite_Basic)
      printf ("***********************************************************\n"
         "Test smarsa_Opso calling smarsa_CollisionOver\n\n");
   smarsa_CollisionOver (gen, res, N, NBalls, r, d, 2);
}


/*=========================================================================*/
/*
 * The CPU time needed for BirthdaySpacings is 6 times longer when I used
 * the standard function qsort of stdlib.h. Thus we use our own QuickSort.
 */

#undef QSORT
#ifdef QSORT
static int compareD (const void *p0, const void *q0)
{
   double x = *((const double *) p0);
   double y = *((const double *) q0);
   return (x < y) ? -1 : (x > y) ? 1 : 0;
}
/* qsort ((void *)(DatDiff + 1), (size_t) n, sizeof (double), compareD); */
#endif


/*=========================================================================*/

static void WriteDataBirth (unif01_Gen * gen, char *TestName, long N, long n,
   int r, long d, int t, int p, double k, smultin_CellType kc,
   double Lambda)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",    d = %1ld,    t = %1d,    p = %1d\n\n", d, t, p);
#ifdef USE_LONGLONG
   if (kc == 0 && d > 1)    /* kc = 2^64 */
      printf ("\n      Number of cells = d^t = " STR_MAXK "\n");
   else
      printf ("\n      Number of cells = d^t = %18" PRIuLEAST64 "\n", kc);
#else
   printf ("\n      Number of cells = d^t = %16.0f\n", k);
#endif
   printf ("      Lambda = Poisson mean = ");
   num_WriteD (Lambda, 12, 4, 2);
   printf ("\n\n");
}


/*-------------------------------------------------------------------------*/

void smarsa_BirthdaySpacings (unif01_Gen *gen, sres_Poisson *res,
   long N, long n, int r, long d, int t, int Order)
{
   long Seq;                      /* Replication number */
   long j;
   long Sum;
   double Y;                      /* Number of collisions */
   double k;
   smultin_CellType kc;
   double Lambda;                 /* Poisson mean */
   smultin_CellType *Dates, *DatDiff;
   fmass_INFO Mass;
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "smarsa_BirthdaySpacings test";

   Timer = chrono_Create ();
   kc = k = d;
   for (j = 2; j <= t; j++) {
      k *= d;
      kc *= d;
   }
   Lambda = (double) n * n / k * (n / 4.0);

   if (swrite_Basic)
      WriteDataBirth (gen, TestName, N, n, r, d, t, Order, k, kc, Lambda);

   if (d <= 1) {
      util_Warning (TRUE,
                    "smarsa_BirthdaySpacings:   d <= 1.  The test is not done.");
      return;
   }
   if (k > smarsa_Maxk) {
      util_Warning (TRUE,
        "smarsa_BirthdaySpacings:   d^t > smarsa_Maxk.  The test is not done.");
      return;
   }
   if (8.0 * N * Lambda > sqrt (sqrt (k))) {
      util_Warning (TRUE,
        "smarsa_BirthdaySpacings:   8N Lambda > k^(1/4).  The test is not done.");
      return;
   }
   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreatePoisson ();
   }
   sres_InitPoisson (res, N, Lambda, "smarsa_BirthdaySpacings");

   Dates = util_Calloc (1 + (size_t) n, sizeof (smultin_CellType));
   DatDiff = util_Calloc (1 + (size_t) n, sizeof (smultin_CellType));

   sprintf (str, "The N statistic values (a Poisson with mean %g):", Lambda);
   statcoll_SetDesc (res->sVal1, str);

   Sum = 0;
   for (Seq = 1; Seq <= N; Seq++) {
      /* Generate and sort the "birth dates" */
      if (Order == 2) {
         for (j = 1; j <= n; j++) {
            Dates[j] = smultin_GenerCellSerial2 (gen, r, t, d);
         }
      } else {
         for (j = 1; j <= n; j++) {
            Dates[j] = smultin_GenerCellSerial (gen, r, t, d);
         }
      }
#ifdef USE_LONGLONG 
      tables_QuickSortULL (Dates, 1, n);
      /* Compute the differences between adjacent dates */
      gofs_DiffULL (Dates, DatDiff, 1, n, 0ULL, 1ULL);
      /* The last cell is a special case */
      DatDiff[n] = kc - Dates[n] + Dates[1];
      tables_QuickSortULL (DatDiff, 1, n);
#else
      tables_QuickSortD (Dates, 1, n);
      /* Compute the differences between adjacent dates */
      gofs_DiffD (Dates, DatDiff, 1, n, 0.0, 1.0);
      /* The last cell is a special case */
      DatDiff[n] = kc - Dates[n] + Dates[1];
      tables_QuickSortD (DatDiff, 1, n);
#endif

      /* Count the number of collisions in DatDiff */
      Y = 0.0;
      for (j = 2; j <= n; j++) {
         if (DatDiff[j] == DatDiff[j - 1])
            Y += 1.0;
      }
      Sum += Y;
      statcoll_AddObs (res->sVal1, Y);
      if (swrite_Counters) {
#ifdef USE_LONGLONG
         tables_WriteTabULL (Dates, 1, n, 3, 21, "Birthdates:");
         tables_WriteTabULL (DatDiff, 1, n, 3, 21, "Birthdate differences:");
#else
         tables_WriteTabD (Dates, 1, n, 4, 17, 0, 0, "Birthdates:");
         tables_WriteTabD (DatDiff, 1, n, 4, 17, 0, 0,
            "Birthdate differences:");
#endif
      }
   }

   res->sVal2 = Sum;
   Mass = fmass_CreatePoisson (N * Lambda);
   res->pLeft = fdist_Poisson2 (Mass, Sum);
   res->pRight = fbar_Poisson2 (Mass, Sum);
   fmass_DeletePoisson (Mass);
   res->pVal2 = gofw_pDisc (res->pLeft, res->pRight);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 1, 1);
   if (swrite_Basic) {
      WriteResultsPoisson (res, N);
      swrite_Final (gen, Timer);
   }
   util_Free (Dates);
   util_Free (DatDiff);
   if (localRes)
      sres_DeletePoisson (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataCAT (unif01_Gen *gen, char *TestName,
   long N, long n, int r, long d, int t, long S[], double Lambda)
{
   int i;
   swrite_Head (gen, TestName, N, n, r);
   printf (",    d = %1ld,    t = %1d\n\n", d, t);
   for (i = 0; i < t; i++) {
      printf ("      S[%1d] =  %1ld\n", i, S[i]);
   }
   printf ("\n      Lambda = Poisson mean = ");
   num_WriteD (Lambda, 12, 4, 2);
   printf ("\n\n");
}


/*-------------------------------------------------------------------------*/

static void TestCATData (long d, int t, long S1[])
/*
 * Test that the key to search for has no overlap, that is cannot be
 * written as ABA, where A and B are parts of the key.
 */
{
   int i, j, s;
   long k1, k2;
   i = 0;
   j = t - 1;
   k1 = k2 = 0;
   while (i < j) {
      k1 = k1 * d + S1[i];
      k2 = 0;
      for (s = j; s < t; s++)
         k2 = k2 * d + S1[s];
      util_Assert (k1 != k2,
         "CATData:   target cell number of the form ABA");
      i++;
      j--;
   }
}


/*-------------------------------------------------------------------------*/
#if 0
static void CATGenere1 (
   unif01_Gen *gen, 
   long n,               /* Number of points */
   int r,                /* Drop the first r bits of each U01 */
   long d,               /* Number of segments on the 1-dim. line */
   int t,                /* Dimension */
   long Key,             /* Key to search for */
   long k1,              /* = d^(t-1) */
   long *Count           /* Number of times Key appears */
   )
/*
 * Generate the n points in the dense case and count the number of times
 * cell Key appears. This is the circular version with n points. It also
 * correspond to the case of aperiodic Key.
 */
{
   int j, i;
   long Indice = 0;
   long Y = 0;                    /* Counter */
   long Premier[32];

   util_Assert (t <= 32, "smarsa_CAT.Genere:   t > 32");

   /* Generation of the first (t - 1) elements of the first tuple */
   for (j = 1; j < t; j++) {
      Premier[j] = unif01_StripL (gen, r, d);
      Indice = Indice * d + Premier[j];
   }

   /* Generation of the n - (tt - 1) tuples */
   for (j = 1; j <= n - (t - 1); j++) {
      /* Remove the leftmost component ... */
      Indice %= k1;
      /* ... shift and get another for the rightmost one */
      Indice = Indice * d + unif01_StripL (gen, r, d);
      if (Indice == Key) {
         ++Y;
         /* Key found: jump over the whole Indice and restart */
         Indice = 0;
         for (i = 1; i < t; i++) {
            Indice = Indice * d + unif01_StripL (gen, r, d);
            j++;
         }
      }
   }

   /* Generation of the last (t - 1) tuples. We use numbers in array */
   /* Premier[] so that the sequence is in fact circular */
   for (j = 1; j < t; j++) {
      Indice %= k1;
      Indice = Indice * d + Premier[j];
      if (Indice == Key)
         ++Y;
   }

   *Count = Y;
}
#endif

/*-------------------------------------------------------------------------*/

static void CATGenere (
   unif01_Gen *gen, 
   long n,               /* Number of points */
   int r,                /* Drop the first r bits of each U01 */
   long d,               /* Number of segments on the 1-dim. line */
   int t,                /* Dimension */
   long Key,             /* Key to search for */
   long k1,              /* = d^(t-1) */
   long *Count           /* Number of times Key appears */
   )
/*
 * Generate the n points in the dense case and count the number of times
 * cell Key appears. This is the non-circular version with n - t + 1 points.
 * It also correspond to the case of aperiodic Key.
 */
{
   int j, i;
   long Indice;
   long Y = 0;                    /* Counter */

   /* Generation of the first (t - 1) elements of the first tuple */
   Indice = 0;
   for (j = 1; j < t; j++)
      Indice = Indice * d + unif01_StripL (gen, r, d);

   /* Generation of the n - (tt - 1) tuples */
   for (j = 1; j <= n - (t - 1); j++) {
      /* Remove the leftmost component ... */
      Indice %= k1;
      /* ... shift and get another for the rightmost one */
      Indice = Indice * d + unif01_StripL (gen, r, d);
      if (Indice == Key) {
         ++Y;
         /* Key found: jump over the whole Indice and restart */
         Indice = 0;
         for (i = 1; i < t; i++) {
            Indice = Indice * d + unif01_StripL (gen, r, d);
            j++;
         }
      }
   }

   *Count = Y;
}


/*-------------------------------------------------------------------------*/

void smarsa_CAT (unif01_Gen *gen, sres_Poisson *res,
   long N, long n, int r, long d, int t, long S[])
{
   long Seq;
   long i;
   double k;
   long k1;                       /* d^(t-1) */
   long Key;                      /* Cell number to search for */
   double Lambda;                 /* Poisson mean */
   long Sum;
   long Co;
   fmass_INFO Mass;
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "smarsa_CAT test";

   Timer = chrono_Create ();
   k1 = d;
   for (i = 2; i < t; i++)
      k1 *= d;
   k = k1 * d;
   Lambda = (n - t + 1) / k;
   if (swrite_Basic)
      WriteDataCAT (gen, TestName, N, n, r, d, t, S, Lambda);
   util_Assert (d > 1, "smarsa_CAT:   d <= 1");

   Key = 0;
   for (i = 0; i < t; i++) {
      if (S[i] < 0 || S[i] >= d) {
         util_Error ("smarsa_CAT:   S[i] must be in [0, d - 1]");
      }
      Key = Key * d + S[i];
   }
   TestCATData (d, t, S);
   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreatePoisson ();
   }
   sres_InitPoisson (res, N, Lambda, "smarsa_CAT");
   sprintf (str, "The N statistic values (a Poisson with mean %g):", Lambda);
   statcoll_SetDesc (res->sVal1, str);

   Sum = 0;
   for (Seq = 1; Seq <= N; Seq++) {
      CATGenere (gen, n, r, d, t, Key, k1, &Co);
      statcoll_AddObs (res->sVal1, (double) Co);
      Sum += Co;
   }

   res->sVal2 = Sum;
   Mass = fmass_CreatePoisson (res->Mu);
   res->pLeft = fdist_Poisson2 (Mass, Sum);
   res->pRight = fbar_Poisson2 (Mass, Sum);
   fmass_DeletePoisson (Mass);
   res->pVal2 = gofw_pDisc (res->pLeft, res->pRight);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 1, 1);
   if (swrite_Basic) {
      WriteResultsPoisson (res, N);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeletePoisson (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataCATBits (unif01_Gen *gen, char *TestName,
   long N, long n, int r, int s, int L, unsigned long Key, double Lambda)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   s = %1d,   L = %1d,   Key = %lu\n\n", s, L, Key);
   printf ("      Lambda = Poisson mean = ");
   num_WriteD (Lambda, 12, 4, 2);
   printf ("\n\n");
}


/*-------------------------------------------------------------------------*/

static void TestCATBitsData (int L, unsigned long Key)
/*
 * Test that the key to search for has no overlap, that is cannot be
 * written as ABA, where A and B are parts of the key.
 */
{
   int i;
   unsigned long mask = 1, shift = L - 1;
   i = 0;
   while (i < L / 2) {
      if ((Key & mask) == (Key >> shift)) {
         bitset_WriteSet ("Key =  ", Key, L);
         util_Error ("CATBitsData:   Key of the form ABA");
      }
      i++;
      shift--;
      mask = num_TwoExp[i + 1] - 1.0;
   }
}


/*-------------------------------------------------------------------------*/

static void CATGenerBits (unif01_Gen *gen, long n, int r, int s, 
   int L, unsigned long KEY0, long *Count)
{
/*
 * Generate the bits in the CATBits test. Points are generated with
 * overlapping. We have a window of size L bits, and we slide it 1 bit
 * forward at each step to generate a point. We then check whether it
 * equals the L bits Key.
 */
   const unsigned long MASK0 = num_TwoExp[L] - 1.0;
   unsigned long Mask, Key, Z0, Z;
   int j0, j, k;
   long i;
   long co;

   util_Assert (L <= 32, "CATBits:   GenerBits:   L > 32");
   co = 0;

   if ((s >= L) && (L <= 16)) {
      const int q = s - L;

      /* Make sure to skip the first half of the loop for the first number
         since there is no previous Z */
      j0 = L;
      Z = 0;

      for (i = 0; i < n / s; i++) {
         Z0 = unif01_StripB (gen, r, s);

         /* The last L - j0 bits of the previous number */
         Mask = MASK0 << (L - j0);
         Key = KEY0 << (L - j0);
         Z |= (Z0 >> q);
         j = j0;
         while (j < L) {
            if (Key == (Z & Mask)) {
               co++;
               j += L;
               Mask >>= L;
               Key >>= L;
            } else {
               j++;
               Mask >>= 1;
               Key >>= 1;
            }
         }
         j0 = j % L;

         /* The first s - L bits of the current number */
         Z = Z0;
         Mask = MASK0 << (q - j0);
         Key = KEY0 << (q - j0);
         j = j0;
         while (j < q) {
            if (Key == (Z & Mask)) {
               co++;
               j += L;
               Mask >>= L;
               Key >>= L;
            } else {
               j++;
               Mask >>= 1;
               Key >>= 1;
            }
         }
         j0 = j - q;
         Z = Z0 << L;
      }

   } else if (s >= L) {
#ifdef USE_LONGLONG
      const ulonglong MASK0 = num_TwoExp[L] - 1.0;
      ulonglong Z, Z0;
      ulonglong Mask, Key;
      const int q = s - L;

      /* Make sure to skip the first half of the loop for the first number
         since there is no previous Z */
      j0 = L;
      Z = 0;

      for (i = 0; i < n / s; i++) {
         Z0 = unif01_StripB (gen, r, s);

         /* The last L - j0 bits of the previous number */
         Mask = MASK0 << (L - j0);
         Key = KEY0 << (L - j0);
         Z |= (Z0 >> q);
         j = j0;
         while (j < L) {
            if (Key == (Z & Mask)) {
               co++;
               j += L;
               Mask >>= L;
               Key >>= L;
            } else {
               j++;
               Mask >>= 1;
               Key >>= 1;
            }
         }
         j0 = j % L;

         /* The first s - L bits of the current number */
         Z = Z0;
         Mask = MASK0 << (q - j0);
         Key = KEY0 << (q - j0);
         j = j0;
         while (j < q) {
            if (Key == (Z & Mask)) {
               co++;
               j += L;
               Mask >>= L;
               Key >>= L;
            } else {
               j++;
               Mask >>= 1;
               Key >>= 1;
            }
         }
         j0 = j - q;
         Z = Z0 << L;
      }
#else
      if (L <= s)
         util_Error ("CATGenerBits:   L <= s and L > 16");
#endif

   } else if ((s < L) && (L + s <= 32)) {
      const int t = L / s;
      util_Assert (L % s == 0, "CATBits:   L > s but L % s not 0");

      /* Generation of the first L random bits */
      Z = 0;
      for (j = 0; j < t; j++) {
         Z <<= s;
         Z |= unif01_StripB (gen, r, s);
      }
      j0 = 0;

      /* Generation of the rest of the random bits */
      for (i = 0; i < (n - L) / s; i++) {
         Z = (Z << s) | unif01_StripB (gen, r, s);
         Mask = MASK0 << (s - j0);
         Key = KEY0 << (s - j0);
         j = j0;
         while (j < s) {
            if (Key == (Z & Mask)) {
               co++;
               j += L;
               i += t - 1;
               for (k = 1; k < t; k++) {
                  Z <<= s;
                  Z |= unif01_StripB (gen, r, s);
               }
            } else {
               j++;
               Mask >>= 1;
               Key >>= 1;
            }
         }
         j0 = j % s;
      }

   } else {
#ifdef USE_LONGLONG
      const ulonglong MASK0 = num_TwoExp[L] - 1.0;
      const int t = L / s;
      ulonglong Z;
      ulonglong Mask, Key, Key0 = KEY0;

      if (L > s) {
         util_Assert (L % s == 0, "CATBits:   L > s but L % s not 0");
      }

      /* Generation of the first L random bits */
      Z = 0;
      for (j = 0; j < t; j++) {
         Z <<= s;
         Z |= unif01_StripB (gen, r, s);
      }
      j0 = 0;

      /* Generation of the rest of the random bits */
      for (i = 0; i < (n - L) / s; i++) {
         Z = (Z << s) | unif01_StripB (gen, r, s);
         Mask = MASK0 << (s - j0);
         Key = Key0 << (s - j0);
         j = j0;
         while (j < s) {
            if (Key == (Z & Mask)) {
               co++;
               j += L;
               i += t - 1;
               for (k = 1; k < t; k++) {
                  Z <<= s;
                  Z |= unif01_StripB (gen, r, s);
               }
            } else {
               j++;
               Mask >>= 1;
               Key >>= 1;
            }
         }
         j0 = j % s;
      }
#else
      if (L == s)
         util_Error ("CATGenereBits:   L = s and s > 16");
      else
         util_Error ("CATGenereBits:   L > s and L + s > 32");
#endif
   }

   *Count = co;
}


/*-------------------------------------------------------------------------*/

void smarsa_CATBits (unif01_Gen *gen, sres_Poisson *res,
   long N, long n, int r, int s, int L, unsigned long Key)
{
   long Seq;
   double Lambda;                 /* Poisson mean */
   long Sum;
   long Co;
   fmass_INFO Mass;
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "smarsa_CATBits test";

   Timer = chrono_Create ();
   Lambda = (n - L + 1) / num_TwoExp[L];
   if (swrite_Basic)
      WriteDataCATBits (gen, TestName, N, n, r, s, L, Key, Lambda);
   util_Assert (L > 1, "smarsa_CATBits:   L <= 1");

   TestCATBitsData (L, Key);
   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreatePoisson ();
   }
   sres_InitPoisson (res, N, Lambda, "smarsa_CATBits");
   sprintf (str, "The N statistic values (a Poisson with mean %g):", Lambda);
   statcoll_SetDesc (res->sVal1, str);

   Sum = 0;
   for (Seq = 1; Seq <= N; Seq++) {
      CATGenerBits (gen, n, r, s, L, Key, &Co);
      statcoll_AddObs (res->sVal1, (double) Co);
      Sum += Co;
   }

   res->sVal2 = Sum;
   Mass = fmass_CreatePoisson (res->Mu);
   res->pLeft = fdist_Poisson2 (Mass, Sum);
   res->pRight = fbar_Poisson2 (Mass, Sum);
   fmass_DeletePoisson (Mass);
   res->pVal2 = gofw_pDisc (res->pLeft, res->pRight);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 1, 1);
   if (swrite_Basic) {
      WriteResultsPoisson (res, N);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeletePoisson (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataMatRank (unif01_Gen * gen, char *TestName,
   long N, long n, int r, int s, int L, int k)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",    s = %1d,    L = %1d,    k = %1d\n\n", s, L, k);
}


/*-------------------------------------------------------------------------*/
#if 0

static int RankOfBitMatrix (bitset_BitSet M[], int maxrow)
/*
 * Calculation of the rank of the bit-matrix M
 */
{
   const int MaxBit = 31;         /* number of bits in a word - 1 */
   bitset_BitSet Swap;
   int rank = 0;
   int i;
   int CL = 1;

   while (CL <= MaxBit) {
      /* All components of M shift their bits 1 position to the left */
      for (i = 0; i < maxrow; i++)
         M[i] <<= 1;

      /* Search of the first M[i] with 1 as the major bit */
      i = rank;
      for (;;) {
         if ((bitset_TestBit (M[i], MaxBit)) || (i == maxrow - 1))
            break;
         ++i;
      }
      /* Diagonalization of matrix M */
      if (i < maxrow - 1) {
         Swap = M[rank];
         M[rank] = M[i];
         M[i] = Swap;
         for (i = rank + 1; i < maxrow; i++) {
            if (bitset_TestBit (M[i], MaxBit))
               M[i] ^= M[rank];
         }
         ++rank;
         if (rank == MaxBit)
            return rank;
      }
      ++CL;
   }
   return rank;
}


/*-------------------------------------------------------------------------*/

#define lmax 64

void smarsa_MatrixRank (unif01_Gen *gen, sres_Chi2 *res,
   long N, long n, int r, int s, int l, int k)
{
   long Seq;
   long Rep;
   int j;
   int i;
   long L;                        /* One line of bits */
   int c;                         /* Number-1 of U01 used to build a line */
   long d;                        /* Get s bits of a generated U01 */
   long a;                        /* Get b bits of a generated U01 */
   int b;                         /* Number of bits of last U01 of a line */
   int Minkl;                     /* Min (k, l) */
   long NbGroups;                 /* Number of classes for ChiSquare */
   long jhigh;                    /* Index of the highest class */
   long jlow;                     /* Index of the lowest class */
   int Rank;                      /* Rank of matrix */
   double X2;
   double Prod;
   long *Loca;                    /* Redirections in merging Chi2 classes */
   long *Count;                   /* Observed numbers */
   double *NbExp;                 /* Expected numbers */
   bitset_BitSet M[lmax];         /* Matrix */
   double V[1];                   /* Number of degrees of freedom for Chi2 */
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "smarsa_MatrixRank test";

   Timer = chrono_Create ();
   /* We shall need c + 1 random numbers to build a line of the matrix */
   c = k / s;
   b = k % s;
   a = num_TwoExp[b];
   d = num_TwoExp[s];
   if (swrite_Basic)
      WriteDataMatRank (gen, TestName, N, n, r, s, l, k);
   if (k <= l)
      Minkl = k;
   else
      Minkl = l;
   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, Minkl, "smarsa_MatrixRank");
   NbExp = res->NbExp;
   Count = res->Count;
   Loca = res->Loc;

   Prod = n * pow (2.0, -(double) (l * k));
   NbExp[0] = Prod;
   for (j = 1; j <= Minkl; j++) {
      Prod = Prod * pow (2.0,  (double) (l + k - 2*j + 1)) *
                (1.0 - 1.0 / num_TwoExp[l - j + 1]) *
                (1.0 - 1.0 / num_TwoExp[k - j + 1]) /
                (1.0 - 1.0 / num_TwoExp[j]);
      NbExp[j] = Prod;
   }

   jlow = 0;
   jhigh = Minkl;
   if (swrite_Classes)
      gofs_WriteClasses (NbExp, Loca, jlow, jhigh, 0);
   gofs_MergeClasses (NbExp, Loca, &jlow, &jhigh, &NbGroups);
   if (swrite_Classes)
      gofs_WriteClasses (NbExp, Loca, jlow, jhigh, NbGroups);
   res->jmin = jlow;
   res->jmax = jhigh;
   res->degFree = NbGroups - 1;

   util_Assert (n > 2.0 * gofs_MinExpected,
      "smarsa_MatrixRank:    n <= 2*gofs_MinExpected");
   util_Assert (k <= 31, "smarsa_MatrixRank:   k > 31");
   util_Assert (l <= lmax, "smarsa_MatrixRank:   L > 64");
   util_Assert (l * k <= 1020, "smarsa_MatrixRank:   L*k > 1020");
   util_Assert (NbGroups > 1,
      "smarsa_MatrixRank:   number of classes = 1."
      "   Increase  n  or decrease  |L - k|");

   sprintf (str, "The N statistic values (a ChiSquare with %1ld degrees"
                 " of freedom):", NbGroups - 1);
   statcoll_SetDesc (res->sVal1, str);

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = jlow; i <= jhigh; i++)
         Count[i] = 0;
      for (Rep = 1; Rep <= n; Rep++) {
         /* Generate the l x k matrix and compute its rank */
         for (i = 0; i < l; i++) {
            /* Build one line of bits L */
            L = 0;
            for (j = 1; j <= c; j++)
               /* Generate s bits */
               L = d * L + unif01_StripB (gen, r, s);
            /* The last b bits of a line of the matrix */
            if (a > 1)
               L = a * L + unif01_StripB (gen, r, b);
            M[i] = L;
         }
         /* Set all remaining lines to 0 */
         for (i = l; i < lmax; i++)
            M[i] = 0;

         Rank = RankOfBitMatrix (M, lmax);
         ++Count[Loca[Rank]];
      }

      X2 = gofs_Chi2 (NbExp, Count, jlow, jhigh);
      statcoll_AddObs (res->sVal1, X2);
      if (swrite_Counters)
         tables_WriteTabL (Count, jlow, jhigh, 5, 12, "Observed Numbers");
   }

   V[0] = NbGroups - 1;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
                      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetChi2SumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   /* !!!! Attention, this Write must use the right pVal */
   if (swrite_Basic) {
      swrite_AddStrChi (str, LENGTH + 1, res->degFree);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTest (N, res);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
}

#endif

/*=========================================================================*/
#if 0
static void ZeroMat (Matrix * M)
{
   int i;

   for (i = 0; i < M->nblignes; i++)
      PutBVToZero (&(M->lignes[i][0]));
}
#endif

/*-------------------------------------------------------------------------*/

void smarsa_MatrixRank (unif01_Gen *gen, sres_Chi2 *res,
   long N, long n, int r, int s, int l, int k)
{
   long Seq;
   long Rep;
   int j;
   int i;
   int c;                         /* Number-1 of U01 used to build a line */
   int b;                         /* Number of bits of last U01 of a line */
   unsigned long bmask;           /* b bits mask */
   unsigned long smask;           /* s bits mask */
   int Minkl;                     /* Min (k, l) */
   long NbGroups;                 /* Number of classes for ChiSquare */
   long jhigh;                    /* Index of the highest class */
   long jlow;                     /* Index of the lowest class */
   int Rank;                      /* Rank of matrix */
   double X2;
   double temp;
   long *Loca;                    /* Redirections in merging Chi2 classes */
   long *Count;                   /* Observed numbers */
   double *NbExp;                 /* Expected numbers */
   double Par[1];                 /* Number of degrees of freedom for Chi2 */
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "smarsa_MatrixRank test";
   Matrix *M;
   BitVect *V;

   Timer = chrono_Create ();
   /* We shall need ceiling(c) random numbers to build a line of the matrix */
   c = k / s;
   b = k % s;
   bmask = num_TwoExp[b] - 1.0;
   /* The b most significant bits are set */
   bmask <<= vectorsF2_WL - b;
   smask = num_TwoExp[s] - 1.0;
   /* The s most significant bits are set */
   smask <<= vectorsF2_WL - s;

   if (swrite_Basic)
      WriteDataMatRank (gen, TestName, N, n, r, s, l, k);
   Minkl = util_Min (k, l);
   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, Minkl, "smarsa_MatrixRank");
   NbExp = res->NbExp;
   Count = res->Count;
   Loca = res->Loc;

   temp = num_Log2((double) n) - l * k;
   NbExp[0] = pow (2.0, temp);
   for (j = 1; j <= Minkl; j++) {
      temp += l + k - 2*j + 1 +
	      num_Log2(1.0 - pow (2.0, -(double) (l - j + 1))) +
	      num_Log2(1.0 - pow (2.0, -(double) (k - j + 1))) -
	      num_Log2(1.0 - pow (2.0, -(double) j));
      NbExp[j] = pow (2.0, temp);
   }

   jlow = 0;
   jhigh = Minkl;
   if (swrite_Classes)
      gofs_WriteClasses (NbExp, Loca, jlow, jhigh, 0);
   gofs_MergeClasses (NbExp, Loca, &jlow, &jhigh, &NbGroups);
   if (swrite_Classes)
      gofs_WriteClasses (NbExp, Loca, jlow, jhigh, NbGroups);
   res->jmin = jlow;
   res->jmax = jhigh;
   res->degFree = NbGroups - 1;

   util_Warning (NbGroups <= 1,
      "smarsa_MatrixRank:   number of Chi2 classes = 1.\n"
      "   Increase  n  or decrease  |L - k|.");
   if (res->degFree < 1) {
      if (localRes)
         sres_DeleteChi2 (res);
      return;
   }
   util_Assert (n >= 2.0 * gofs_MinExpected,
      "smarsa_MatrixRank:    n <= 2*gofs_MinExpected");

   sprintf (str, "The N statistic values (a ChiSquare with %1ld degrees"
                 " of freedom):", NbGroups - 1);
   statcoll_SetDesc (res->sVal1, str);

   M = util_Malloc (sizeof (Matrix));
   AllocMat (M, l, k, 1);

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = jlow; i <= jhigh; i++)
         Count[i] = 0;

      for (Rep = 1; Rep <= n; Rep++) {
         /* Generate the l x k matrix and compute its rank */
         for (i = 0; i < l; i++) {
            V = &(M->lignes[i][0]);
            /* Build one line of bits */
            for (j = 0; j < c; j++) {
               /* Shift by s and generate s new bits */
               BVRShiftSelf (V, s);
               V->vect[0] |= (smask &
                  (gen->GetBits (gen->param, gen->state) << r));
            }
            /* The last b bits of a line of the matrix */
            if (b > 0) {
               BVRShiftSelf (V, b);
               V->vect[0] |= (bmask &
                  (gen->GetBits (gen->param, gen->state) << r));
            }
         }
         Rank = GaussianElimination (M, l, k, 1);
         ++Count[Loca[Rank]];
      }

      X2 = gofs_Chi2 (NbExp, Count, jlow, jhigh);
      statcoll_AddObs (res->sVal1, X2);
      if (swrite_Counters)
         tables_WriteTabL (Count, jlow, jhigh, 5, 12, "Observed Numbers");
   }

   FreeMat (M);
   util_Free (M);

   Par[0] = NbGroups - 1;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, Par,
                      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetChi2SumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   /* !!!! Attention, this Write must use the right pVal */
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

static void WriteDataSavir2 (unif01_Gen * gen, char *TestName,
   long N, long n, int r, long m, int t)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",    m = %1ld,    t = %1d\n\n", m, t);
}


/*-------------------------------------------------------------------------*/

void smarsa_Savir2 (unif01_Gen *gen, sres_Chi2 *res,
   long N, long n, int r, long m, int t)
{
   const double eps = 1.0E-15;
   long I;
   long msup;                     /* Dimension - 1 of arrays */
   long Seq;
   long Rep;
   long j;
   int i;
   long NbGroups;                 /* Number of classes for ChiSquare */
   long jhigh;                    /* Index of the highest class */
   long jlow;                     /* Index of the lowest class */
   double X2;                     /* ChiSquare Statistic */
   double UnSurm = 1.0 / m;
   double *Prob;                  /* Probabilities */
   long *Loca;
   double V[1];                   /* Number degrees of freedom for Chi2 */
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "smarsa_Savir2 test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataSavir2 (gen, TestName, N, n, r, m, t);

   Prob = util_Calloc ((size_t) m + 2, sizeof (double));
   Prob[m + 1] = 0.0;
   for (j = 1; j <= m; j++)
      Prob[j] = UnSurm;
   for (i = 2; i <= t; i++) {
      for (j = m; j >= 1; j--)
         Prob[j] = Prob[j + 1] + Prob[j] / j;
   }
   j = 1;
   while (Prob[j] > eps)
      ++j;
   msup = j - 1;

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, msup, "smarsa_Savir2");

   for (j = 1; j <= msup; j++)
      res->NbExp[j] = Prob[j] * n;
   util_Free (Prob);
   Loca = res->Loc;

   jlow = 1;
   jhigh = msup;
   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, Loca, jlow, jhigh, 0);
   gofs_MergeClasses (res->NbExp, Loca, &jlow, &jhigh, &NbGroups);
   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, Loca, jlow, jhigh, NbGroups);
   res->jmin = jlow;
   res->jmax = jhigh;
   res->degFree = NbGroups - 1;

   util_Warning (NbGroups < 2,
      "smarsa_Savir2:   Number of classes = 1.\n   Decrease t or increase n.");
   if (res->degFree < 1) {
      if (localRes)
         sres_DeleteChi2 (res);
      return;
   }
   util_Assert (n >= 2.0 * gofs_MinExpected,
      "smarsa_Savir2:    n <= 2*gofs_MinExpected");

   sprintf (str, "The N statistic values (a ChiSquare with %1ld degrees"
                 " of freedom):", NbGroups - 1);
   res->sVal1 = statcoll_Create (N, str);

   for (Seq = 1; Seq <= N; Seq++) {
      for (j = jlow; j <= jhigh; j++)
         res->Count[j] = 0;
      for (Rep = 1; Rep <= n; Rep++) {
         I = m;
         for (i = 1; i <= t; i++)
            I = 1 + unif01_StripD (gen, r) * I;
         if (I > msup)
            ++res->Count[Loca[msup]];
         else
            ++res->Count[Loca[I]];
      }

      if (swrite_Counters)
         tables_WriteTabL (res->Count, jlow, jhigh, 5, 12,
            "Observed Numbers");

      X2 = gofs_Chi2 (res->NbExp, res->Count, jlow, jhigh);
      statcoll_AddObs (res->sVal1, X2);
   }

   V[0] = NbGroups - 1;
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

static void WriteDataGCD (unif01_Gen * gen, char *TestName,
   long N, long n, int r, int s)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   s = %1d\n\n", s);
}


/*-------------------------------------------------------------------------*/

void smarsa_GCD (unif01_Gen *gen, smarsa_Res2 *res,
                 long N, long n, int r, int s)
{
  /*
   The theoretical distribution for the number of iterations is unknown.
   The binomial is a very rough approximation: thus the printing of the
   results is commented out.
   */
   const double C1 = 6 / (num_Pi * num_Pi);
   const double P1 = 0.376;
   const int KMAX = 50;
   unsigned long U, V, temp;
   double X;
   double Param[1];
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "smarsa_GCD test";
   sres_Chi2 *GCD;
   sres_Chi2 *NumIter;
   int jmax, j, k;
   long Seq, i;
   double *NbExp;
   long *Loc;
   long NbClasses;
   fmass_INFO Q;

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataGCD (gen, TestName, N, n, r, s);
   if (n < 30) {
      util_Warning (TRUE, "n < 30");
      return;
   }
   if (n > pow (2.0, 1.5*s)) {
      util_Warning (TRUE, "n > 2^(1.5s)");
      return;
   } 
   if (res == NULL) {
      localRes = TRUE;
      res = smarsa_CreateRes2 ();
   }
   jmax = 1 + sqrt (C1 * n / gofs_MinExpected);
   util_Assert (jmax > 1, "smarsa_GCD:   jmax < 2");
   InitRes2 (res, N, jmax, KMAX);

   GCD = res->GCD;
   GCD->jmin = 1;
   GCD->jmax = jmax;
   GCD->degFree = jmax - 1;
   sprintf (str, "GCD; the N statistic values (a ChiSquare with %1d degrees"
                 " of freedom):", jmax - 1);
   statcoll_SetDesc (GCD->sVal1, str);

   /* Compute the probabilities for the GCD values */
   NbExp = GCD->NbExp;
   Loc = GCD->Loc;
   X = 0.0;
   for (j = 1; j < jmax; j++) {
      NbExp[j] = n * C1 / ((double) j * j);
      X += NbExp[j];
      Loc[j] = j;
   }
   NbExp[jmax] = n - X;

   if (swrite_Classes) {
      printf ("Classes for the GCD values:\n");
      gofs_WriteClasses (GCD->NbExp, GCD->Count, 1, jmax, 0);
   }

   NumIter = res->NumIter;
   /* Compute expected numbers for number of iterations */
   Q = fmass_CreateBinomial (KMAX, P1, 1.0 - P1);
   for (i = 0; i <= KMAX; i++)
      NumIter->NbExp[i] = n * fmass_BinomialTerm2 (Q, i);
   fmass_DeleteBinomial (Q);

   NumIter->jmin = 0;
   NumIter->jmax = KMAX;
   if (swrite_Classes) {
      printf ("\nClasses for the number of iterations:\n");
      gofs_WriteClasses (NumIter->NbExp, NumIter->Loc, NumIter->jmin,
                         NumIter->jmax, 0);
   }
   gofs_MergeClasses (NumIter->NbExp, NumIter->Loc, &NumIter->jmin,
                      &NumIter->jmax, &NbClasses);

   if (swrite_Classes)
      gofs_WriteClasses (NumIter->NbExp, NumIter->Loc, NumIter->jmin,
                         NumIter->jmax, NbClasses);

   sprintf (str, "NumIter; the N statistic values (a ChiSquare with %1ld"
                 " degrees of freedom):", NbClasses - 1);
   statcoll_SetDesc (NumIter->sVal1, str);
   NumIter->degFree = NbClasses - 1;
   util_Assert (NumIter->degFree >= 1, "NumIter->degFree < 1");

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 0; i <= KMAX; i++)
         NumIter->Count[i] = 0;
      for (i = 0; i <= GCD->jmax; i++)
         GCD->Count[i] = 0;
      for (i = 1; i <= n; i++) {
         k = 0;
         do {
            U = unif01_StripB (gen, r, s);
            V = unif01_StripB (gen, r, s);
         } while (0 == U || 0 == V);
         do {
            temp = U % V;
            U = V;
            V = temp;
            k++;
         } while (V > 0);
         if ((long) U > GCD->jmax)
            U = GCD->jmax;
         (GCD->Count[U])++;
         if (k > KMAX)
            k = KMAX;
         (NumIter->Count[NumIter->Loc[k]])++;
      }
      if (swrite_Counters) {
         tables_WriteTabL (GCD->Count, GCD->jmin, GCD->jmax, 5, 10,
                           "Observed numbers for GCD values:");
 /*         tables_WriteTabL (NumIter->Count, NumIter->jmin, NumIter->jmax, 5,
                           10, "Observed numbers for number of iterations:");
 */
      }

      X = gofs_Chi2 (GCD->NbExp, GCD->Count, GCD->jmin, GCD->jmax);
      statcoll_AddObs (GCD->sVal1, X);
      X = gofs_Chi2 (NumIter->NbExp, NumIter->Count, NumIter->jmin,
                     NumIter->jmax);
      statcoll_AddObs (NumIter->sVal1, X);
   }

   Param[0] = GCD->degFree;
   gofw_ActiveTests2 (GCD->sVal1->V, GCD->pVal1->V, N, wdist_ChiSquare,
                      Param, GCD->sVal2, GCD->pVal2);
   GCD->pVal1->NObs = N;
   sres_GetChi2SumStat (GCD);
/*
   Param[0] = NumIter->degFree;
   gofw_ActiveTests2 (NumIter->sVal1->V, NumIter->pVal1->V, N,
      wdist_ChiSquare, Param, NumIter->sVal2, NumIter->pVal2);
   NumIter->pVal1->NObs = N;
*/

   if (swrite_Basic) {
      if (swrite_Collectors)
         statcoll_Write (GCD->sVal1, 5, 14, 4, 3);
      printf ("\n-----------------------------------------------\n");
      if (N == 1) {
         printf ("Number of degrees of freedom          : %4ld\n",
                  GCD->degFree);
         printf ("Chi2 statistic for GCD values         :");
         gofw_Writep2 (GCD->sVal2[gofw_Mean], GCD->pVal2[gofw_Mean]);
      } else {
         printf ("Test results for GCD values:\n");
         gofw_WriteActiveTests0 (N, GCD->sVal2, GCD->pVal2);
         swrite_Chi2SumTest (N, GCD);
      }
      /*
      if (swrite_Collectors)
         statcoll_Write (NumIter->sVal1, 5, 14, 4, 3);
      printf ("\n-----------------------------------------------\n");
      if (N == 1) {
         printf ("Number of degrees of freedom          : %4ld\n",
                  NumIter->degFree);
	 printf ("Chi2 statistic for NumIter            :");
	 gofw_Writep2 (NumIter->sVal2[gofw_Mean], NumIter->pVal2[gofw_Mean]);
      } else {
	 printf ("Test results for NumIter:\n");
	 gofw_WriteActiveTests0 (N, NumIter->sVal2, NumIter->pVal2);
         swrite_SumTest (N, NumIter->sVal2[gofw_Sum], NumIter->pVal2[gofw_Sum], 
                         N*NumIter->degFree);
      }
      */
      printf ("\n\n");
      swrite_Final (gen, Timer);
   }

   if (localRes)
      smarsa_DeleteRes2 (res);
   chrono_Delete (Timer);

}


/*=========================================================================*/
