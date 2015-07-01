/*************************************************************************\
 *
 * Package:        TestU01
 * File:           swalk.c
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
#include "bitset.h"

#include "swalk.h"
#include "wdist.h"
#include "swrite.h"
#include "unif01.h"

#include "fmass.h"
#include "fbar.h"
#include "gofw.h"
#include "gofs.h"

#include <math.h>
#include <stddef.h>
#include <string.h>




/*------------------------------- Constants -------------------------------*/

#define PREC 52                   /* Max number of bits in a number */
#define LENGTH 200                 /* Max length of strings */


typedef enum {
   swalk_rwH,                         /* H statistic */
   swalk_rwM,                         /* M statistic */
   swalk_rwJ,                         /* J statistic */
   swalk_rwR,                         /* R statistic */
   swalk_rwC,                         /* C statistic */
   swalk_rw_N                         /* Total number of statistics here */
} swalk_rwType;

/* The name of each type of statistic in swalk_rwType (predefined). */
static const char *swalk_rwName[swalk_rw_N] = {
   "Statistic H",
   "Statistic M",
   "Statistic J",
   "Statistic R",
   "Statistic C"
};


/*--------------------------------- Types ---------------------------------*/

typedef struct {
   long X;
   long S;
   long S_2;
   long M;
   long R;
   long J;
   long C;
} WorkType;


/* Type of algorithm used in swalk_VarGeo */
typedef enum {
   swalk_AlgoP,
   swalk_AlgoN
   } swalk_AlgoType;




/*-------------------------------- functions ------------------------------*/

static void CalcNbExp (
   long n,                    /* Sample size */
   long L0,
   long k,
   swalk_Res *res
   )
/* 
 ************   IMPORTANT: we assume that L is even  ************
 * Compute the expected numbers for the different statistics in the
 * swalk_RandomWalk1 and swalk_RandomWalk1a tests. We start from the
 * maximum term and compute on each side all terms larger than epsilon.
 * We set all others to 0.
 */
{
   const double epsilon = 1.0E-16;
   double *NbExp;
   long L1, L2;
   const long L = L0 + k;
   long i;
   double nLR = n;
   double epsn;

   util_Assert (!(L & 1), "CalcNbExp:   L is odd");
   L2 = L / 2;
   epsn = epsilon * nLR;

   /*----------- statistic H -----------*/
   NbExp = res->H[k]->NbExp;
   for (i = 0; i <= L; i++)
      NbExp[i] = 0.0;

   NbExp[L2] = nLR * fmass_BinomialTerm1 (L, 0.5, 0.5, L2);
   i = L2;
   while (i > 0 && NbExp[i] > epsn) {
      NbExp[i - 1] = NbExp[i] * i / (L - i + 1);
      --i;
   }
   i = L2;
   while (i < L && NbExp[i] > epsn) {
      NbExp[i + 1] = NbExp[i] * (L - i) / (i + 1);
      ++i;
   }

   /*----------- statistic M -----------*/
   NbExp = res->M[k]->NbExp;
   for (i = 0; i <= L; i++)
      NbExp[i] = 0.0;

   NbExp[0] = res->H[k]->NbExp[L2];
   i = 0;
   while (i < L && NbExp[i] > epsn) {
      NbExp[i + 1] = NbExp[i] * ((L - i) / 2) / ((L + i) / 2 + 1);
      NbExp[i + 2] = NbExp[i + 1];
      i += 2;
   }

   /*----------- statistic J -----------*/
   NbExp = res->J[k]->NbExp;
   for (i = 0; i <= L; i++)
      NbExp[i] = 0.0;

   NbExp[0] = res->M[k]->NbExp[0];
   NbExp[L] = NbExp[0];
   i = 0;
   while (i < L2 && NbExp[i] > epsn) {
      NbExp[i + 2] = NbExp[i] * ((L - i) / 2) *
         (1 + i) / ((double)(i / 2 + 1) * (L - i - 1));
      NbExp[L - i - 2] = NbExp[i + 2];
      i += 2;
   }

   /*----------- statistic R -----------*/
   NbExp = res->R[k]->NbExp;
   for (i = 0; i <= L; i++)
      NbExp[i] = 0.0;

   NbExp[0] = res->J[k]->NbExp[0];
   i = 0;
   while (i < L2 && NbExp[i] > epsn) {
      NbExp[i + 1] = NbExp[i] * (L - 2 * i) / (L - i);
      ++i;
   }

   /*----------- statistic C -----------*/
   NbExp = res->C[k]->NbExp;
   for (i = 0; i <= L; i++)
      NbExp[i] = 0.0;

   NbExp[0] = 2.0 * nLR * fmass_BinomialTerm1 (L - 1, 0.5, 0.5, L2);
   i = 0;
   L1 = L2 - 1;
   while (i < L1 && NbExp[i] > epsn) {
      NbExp[i + 1] = NbExp[i] * (L2 - i - 1) / (L2 + i + 1);
      ++i;
   }
}


/*-------------------------------------------------------------------------*/

static void WriteTabWalk (
   swalk_Res *res,
   long N                 /* Number of replications */
   )
/*
 * Write the values of the statistics and their p-values in a table
 * for all values of the walk length L from L0 to L1. These values have
 * all been written before, but it is nice to have them all together.
 * When L0 = L1, it is not useful. When the number of replications N > 1, it 
 * is not called either because there would be so many statistics to write.
 */
{
   swalk_rwType m;
   long k;
   double p;
   long L0 = res->L0;             /* Shortest walk length considered */
   long L1 = res->L1;             /* Longest walk length considered */

   if (L1 == L0)
      return;
   if (N > 1)
      return;
   printf ("\n\n***********************************************"
           "\nTABLES FOR THE RESULTS ABOVE");
   for (m = 0; m < swalk_rw_N; m++) {
      printf ("\n\n===============================================\n");
      printf ("Test on the values of the ");
      printf ("%s", swalk_rwName[m]);
      printf ("\n\n  Walk length      Chi-square        p-value\n\n");

      for (k = 0; k <= L1 - L0; k += 2) {
         printf ("%8ld", L0 + k);
         switch (m) {
         case swalk_rwH:
            num_WriteD (res->H[k]->sVal2[gofw_Mean], 18, 3, 2);
            p = res->H[k]->pVal2[gofw_Mean];
            break;
         case swalk_rwM:
            num_WriteD (res->M[k]->sVal2[gofw_Mean], 18, 3, 2);
            p = res->M[k]->pVal2[gofw_Mean];
            break;
         case swalk_rwJ:
            num_WriteD (res->J[k]->sVal2[gofw_Mean], 18, 3, 2);
            p = res->J[k]->pVal2[gofw_Mean];
            break;
         case swalk_rwR:
            num_WriteD (res->R[k]->sVal2[gofw_Mean], 18, 3, 2);
            p = res->R[k]->pVal2[gofw_Mean];
            break;
         case swalk_rwC:
            num_WriteD (res->C[k]->sVal2[gofw_Mean], 18, 3, 2);
            p = res->C[k]->pVal2[gofw_Mean];
            break;
         default:
            util_Error ("swalk:  WriteTabWalk: no such case");
         }
         num_WriteD (p, 18, 3, 2);
         if (p < gofw_Suspectp || p > 1.0 - gofw_Suspectp) {
            printf ("     *****");
         }
         printf ("\n");
      }
   }
   printf ("\n\n");
}


/*-------------------------------------------------------------------------*/

static void WriteResultWalk (
   swalk_Res *res,
   long N                 /* Number of replications */
   )
/*
 * Write the basic results of the swalk_RandomWalk1 and swalk_RandomWalk1a
 * tests for all walks of length L between L0 and L1. May write the
 * statistical collectors too. 
 */
{
   swalk_rwType m;
   long L0 = res->L0;             /* Shortest walk length considered */
   long L1 = res->L1;             /* Longest walk length considered */
   long k;
   sres_Chi2 *Q;

   printf ("\n");
   for (k = 0; k <= L1 - L0; k += 2) {
      if (L1 > L0) {
	 printf ("\n\n==============================================="
	         "\nWALK OF %3ld STEPS\n", L0 + k);
      }
      for (m = 0; m < swalk_rw_N; m++) {
	 printf ("-----------------------------------------------\n"
                 "Test on the values of the ");
         printf ("%s", swalk_rwName[m]);
         printf ("\n\n");
         switch (m) {
         case swalk_rwH:
            Q = res->H[k];
            break;
         case swalk_rwM:
            Q = res->M[k];
            break;
         case swalk_rwJ:
            Q = res->J[k];
            break;
         case swalk_rwR:
            Q = res->R[k];
            break;
         case swalk_rwC:
            Q = res->C[k];
            break;
         default:
            util_Error ("swalk:  WriteResultWalk: no such case");
         }
         if (N == 1) {
            printf ("Number of degrees of freedom          : %4ld\n",
                 Q->degFree);
            printf ("ChiSquare statistic                   :");
            gofw_Writep2 (Q->sVal2[gofw_Mean], Q->pVal2[gofw_Mean]);
         } else {
            gofw_WriteActiveTests0 (N, Q->sVal2, Q->pVal2);
            swrite_Chi2SumTest (N, Q);
         }
         printf ("\n");
         if (swrite_Collectors)
            statcoll_Write (Q->sVal1, 5, 14, 4, 3);
      }
   }
   WriteTabWalk (res, N);
}


/*-------------------------------------------------------------------------*/

static void WriteDetailsWalk (
   swalk_Res *res,
   long k,               /* Walk length L = k + L0 */
   long n                /* sample size */
   )
/*
 * Write detailed results for the different statistics of a walk of length
 * L = k + L0: the expected numbers, the observed numbers (the counters),
 * and the normalized values. If the normalized value is outside the interval
 * [-3, 3], we indicate it explicitly.

 */
{
   swalk_rwType m;
   long i;
   double iObs;                   /* Weighted sum of the observed numbers */
   double Obs;                    /* Observed numbers */
   double iEsp;                   /* Weighted sum of the expected numbers */
   double Esp;                    /* Expected numbers */
   double Z;                      /* Normalized value */
   double Var;                    /* Variance */
   long L0 = res->L0;             /* Shortest length of walks considered */
   double nLR = n;
   sres_Chi2 *Q;

   printf ("================================================\n");
   printf ("Walk of %3ld steps\n", L0 + k);

   for (m = 0; m < swalk_rw_N; m++) {
      printf ("------------------------------------------------\n"
              "Counters of the ");
      printf ("%s", swalk_rwName[m]);
      printf
         ("\n\n  i     Expected num. Observed num.  (Exp. - Obs.)/sigma\n\n");
      iEsp = 0.0;
      iObs = 0.0;
      switch (m) {
      case swalk_rwH:
         Q = res->H[k];
         break;
      case swalk_rwM:
         Q = res->M[k];
         break;
      case swalk_rwJ:
         Q = res->J[k];
         break;
      case swalk_rwR:
         Q = res->R[k];
         break;
      case swalk_rwC:
         Q = res->C[k];
         break;
      default:
         util_Error ("swalk:  WriteDetailsWalk: no such case");
      }
      i = Q->jmin - 1;
      do {
         i = Q->Loc[i + 1];
         Esp = Q->NbExp[i];
         Obs = Q->Count[i];
         iEsp = iEsp + Esp * i;   /* Expected value of the statistic */
         iObs = iObs + Obs * i;   /* Observed mean of the statistic */

         /* If Esp = 0, this is a class that has been redirected to another
            for the ChiSquare test; we shall not print it since the counters 
            have been redirected also: they are necessarily 0. */
         if (Esp > 0.0) {
            printf ("%4ld", i);
            num_WriteD (Esp, 14, 2, 0);
            num_WriteD (Obs, 12, 0, 0);
            Var = Esp * (1.0 - Esp / nLR);
            if (Var <= 0.0)
               Z = (Obs - Esp) * 1.E100;
            else
               Z = (Obs - Esp) / sqrt (Var);
            num_WriteD (Z, 18, 4, 3);
            if (Z > 3.0 || Z < -3.0)
               printf ("    *****");
            printf ("\n");
         }
      } while (i != Q->jmax);

      printf ("\nExpected mean  = ");
      num_WriteD (iEsp / nLR, 10, 2, 0);
      printf ("\nEmpirical mean = ");
      num_WriteD (iObs / nLR, 10, 2, 0);
      printf ("\n\n");
   }
   printf ("\n");
}


/*-------------------------------------------------------------------------*/

static void WriteDataWalk1 (unif01_Gen *gen, char *TestName,
   long N, long n, int r, int s, long L0, long L1)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   s = %1d,   L0 = %4ld,   L1 = %4ld\n\n\n", s, L0, L1);
}


/*=========================================================================*/

swalk_Res * swalk_CreateRes (void)
{
   swalk_Res *res;

   res = util_Malloc (sizeof (swalk_Res));
   memset (res, 0, sizeof (swalk_Res));
   res->H = util_Calloc (1, sizeof (sres_Chi2 *));
   res->M = util_Calloc (1, sizeof (sres_Chi2 *));
   res->J = util_Calloc (1, sizeof (sres_Chi2 *));
   res->R = util_Calloc (1, sizeof (sres_Chi2 *));
   res->C = util_Calloc (1, sizeof (sres_Chi2 *));
   res->H[0] = sres_CreateChi2 ();
   res->M[0] = sres_CreateChi2 ();
   res->J[0] = sres_CreateChi2 ();
   res->R[0] = sres_CreateChi2 ();
   res->C[0] = sres_CreateChi2 ();
   res->imax = 0;
   res->name = util_Calloc (1, sizeof (char));
   return res;
}


/*-------------------------------------------------------------------------*/

void swalk_DeleteRes (swalk_Res *res)
{
   long i;

   if (res == NULL)
      return;
   util_Free (res->name);
   for (i = 0; i <= res->imax; i += 2) {
      sres_DeleteChi2 (res->H[i]);
      sres_DeleteChi2 (res->M[i]);
      sres_DeleteChi2 (res->R[i]);
      sres_DeleteChi2 (res->J[i]);
      sres_DeleteChi2 (res->C[i]);
   }
   util_Free (res->H);
   util_Free (res->R);
   util_Free (res->M);
   util_Free (res->J);
   util_Free (res->C);
   util_Free (res);
}


/*-------------------------------------------------------------------------*/

static void InitRes (
   swalk_Res *res,
   WorkType *work,
   long N,                /* Number of replications */
   long L0,               /* Shortest walk length considered */
   long L1,               /* Longest walk length considered */
   char *nam
)
/*
 * Allocates memory for arrays to be used in the swalk_RandomWalk1a and
 * swalk_RandomWalk1 tests. Arrays for walk length L will be indexed
 * by i = L - L0.
 */
{
   long i, imax, L;

   util_Assert (!(L0 & 1), "InitRes:   L0 is odd");
   if (L1 & 1)
      L1--;
   util_Assert (L1 >= L0, "InitRes:   L1 < L0");
   imax = L1 - L0;

   for (i = imax + 2; i <= res->imax; i += 2) {
      sres_DeleteChi2 (res->H[i]);
      sres_DeleteChi2 (res->M[i]);
      sres_DeleteChi2 (res->R[i]);
      sres_DeleteChi2 (res->J[i]);
      sres_DeleteChi2 (res->C[i]);
   }
   res->H = util_Realloc (res->H, ((size_t) imax + 1) * sizeof(sres_Chi2 *));
   res->R = util_Realloc (res->R, ((size_t) imax + 1) * sizeof(sres_Chi2 *));
   res->M = util_Realloc (res->M, ((size_t) imax + 1) * sizeof(sres_Chi2 *));
   res->J = util_Realloc (res->J, ((size_t) imax + 1) * sizeof(sres_Chi2 *));
   res->C = util_Realloc (res->C, ((size_t) imax + 1) * sizeof(sres_Chi2 *));

   for (i = res->imax + 2; i <= imax; i += 2) {
      res->H[i] = sres_CreateChi2 ();
      res->M[i] = sres_CreateChi2 ();
      res->J[i] = sres_CreateChi2 ();
      res->R[i] = sres_CreateChi2 ();
      res->C[i] = sres_CreateChi2 ();
   }

   for (i = 0; i <= imax; i += 2) {
      L = i + L0;
      sres_InitChi2 (res->H[i], N, L, "");
      sres_InitChi2 (res->M[i], N, L, "");
      sres_InitChi2 (res->R[i], N, L, "");
      sres_InitChi2 (res->J[i], N, L, "");
      sres_InitChi2 (res->C[i], N, L, "");
      res->R[i]->jmax = L / 2;
      res->C[i]->jmax = L / 2;
   }
   res->L1 = L1;
   res->L0 = L0;
   res->imax = imax;
   res->work = work;
   res->name = util_Realloc (res->name, 1 + strlen (nam) * sizeof (char));
   strcpy (res->name, nam);
}


/*-------------------------------------------------------------------------*/

static void Steps (
   unif01_Gen *gen, 
   swalk_Res *res,
   long j,               /* will generate j-th random number */
   int r,                /* drop the first r bits of each random number */
   int s                 /* keep only s bits of each random number */
   )
/*
 * Generates a random number Z; each of the s bits of Z kept is a step of
 * the random walk. Updates all statistics for these s steps.
 */
{
   int i;
   long k;
   unsigned long Z, iBit;
   const unsigned long SBIT = 1UL << (s - 1);
   WorkType *work = res->work;

   Z = unif01_StripB (gen, r, s);
   iBit = SBIT;

   for (i = s - 1; i >= 0; i--) {
      ++res->L;
      if (Z & iBit)               /* If i bit of Z is 1 */
         work->X = 1;
      else
         work->X = -1;
      work->S += work->X;
      if (work->S > work->M)
         work->M = work->S;
      if (work->S == 0)
         ++work->R;
      if ((s * j - i) & 1) {
         if (work->S > 0)
            ++work->J;
         if (work->S * work->S_2 < 0)
            ++work->C;
         work->S_2 = work->S;
      }

      if ((res->L >= res->L0) && !(res->L & 1)) {
         k = res->L - res->L0;
         ++res->H[k]->Count[res->H[k]->Loc[(res->L + work->S) / 2]];
         ++res->M[k]->Count[res->M[k]->Loc[work->M]];
         ++res->J[k]->Count[res->J[k]->Loc[2 * work->J]];
         ++res->R[k]->Count[res->R[k]->Loc[work->R]];
         ++res->C[k]->Count[res->C[k]->Loc[work->C]];
      }
      iBit >>= 1;
   }
}


/*-------------------------------------------------------------------------*/

void swalk_RandomWalk1 (unif01_Gen *gen, swalk_Res *res,
   long N, long n, int r, int s, long L0, long L1)
{
   swalk_rwType m;
   long DeltaL;
   int LMS;
   long LDS;
   long i, j, k, Rep, Seq;
   double khi;
   double V[1];               /* Number degrees of freedom for ChiSquare */
   long NbClasses;
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "swalk_RandomWalk1 test";
   WorkType work;
   sres_Chi2 *Q;

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataWalk1 (gen, TestName, N, n, r, s, L0, L1);

   util_Assert (L0 > 0,    "swalk_RandomWalk1:   L0 <= 0");
   util_Assert (!(L0 & 1), "swalk_RandomWalk1:   L0 must be even");
   util_Assert (!(L1 & 1), "swalk_RandomWalk1:   L1 must be even");
   util_Assert (L1 >= L0,  "swalk_RandomWalk1:   L0 > L1");
   util_Assert (r + s <= PREC, "swalk_RandomWalk1:   r + s > 32");
   if (n < 3.0 * gofs_MinExpected) {
      util_Warning (TRUE, "swalk_RandomWalk1:   n < 3*gofs_MinExpected");
      return;
   }

   DeltaL = L1 - L0;
   LDS = L1 / s;
   LMS = L1 % s;
   if (res == NULL) {
      localRes = TRUE;
      res = swalk_CreateRes ();
   }
   InitRes (res, &work, N, L0, L1, "swalk_RandomWalk1");

   /* Compute the expected numbers and merge classes for the ChiSquare */
   for (k = 0; k <= DeltaL; k += 2) {
      CalcNbExp (n, L0, k, res);

      for (m = 0; m < swalk_rw_N; m++) {
         switch (m) {
         case swalk_rwH:
            Q = res->H[k];
            break;
         case swalk_rwM:
            Q = res->M[k];
            break;
         case swalk_rwJ:
            Q = res->J[k];
            break;
         case swalk_rwR:
            Q = res->R[k];
            break;
         case swalk_rwC:
            Q = res->C[k];
            break;
         default:
            util_Error ("swalk_RandomWalk1:   no such case");
         }
         if (swrite_Classes) {
            if (L1 > L0) {
               printf ("===============================================\n");
               printf ("Walk of %3ld steps\n", L0 + k);
            }
            printf ("===============================================\nThe ");
            printf ("%s", swalk_rwName[m]);
            printf ("\n");
            gofs_WriteClasses (Q->NbExp, Q->Loc, Q->jmin, Q->jmax, 0);
         }
         gofs_MergeClasses (Q->NbExp, Q->Loc, &(Q->jmin), &(Q->jmax),
            &NbClasses);

         if (swrite_Classes) {
            gofs_WriteClasses (Q->NbExp, Q->Loc, Q->jmin, Q->jmax, NbClasses);
         }

         /* Set description for second level statistical collectors */
         sprintf (str, "The N statistic values (a ChiSquare with %ld degrees"
                       " of freedom) ", NbClasses - 1);
         statcoll_SetDesc (Q->sVal1, str);
         Q->degFree = NbClasses - 1;
      }
   }

   /* Beginning of test */
   for (Seq = 1; Seq <= N; Seq++) {

      /* Reset counters to zero for each sequence */
      for (k = 0; k <= DeltaL; k += 2) {
         for (i = res->H[k]->jmin; i <= res->H[k]->jmax; i++)
            res->H[k]->Count[i] = 0;
         for (i = res->R[k]->jmin; i <= res->R[k]->jmax; i++)
            res->R[k]->Count[i] = 0;
         for (i = res->J[k]->jmin; i <= res->J[k]->jmax; i++)
            res->J[k]->Count[i] = 0;
         for (i = res->M[k]->jmin; i <= res->M[k]->jmax; i++)
            res->M[k]->Count[i] = 0;
         for (i = res->C[k]->jmin; i <= res->C[k]->jmax; i++)
            res->C[k]->Count[i] = 0;
      }

      /* A ChiSquare sample of size n */
      for (Rep = 1; Rep <= n; Rep++) {
         work.S = 0;
         work.S_2 = 0;
         work.M = 0;
         work.R = 0;
         work.J = 0;
         work.C = 0;
         res->L = 0;

         /* One random walk of length L */
         for (j = 1; j <= LDS; j++)
            Steps (gen, res, j, r, s);
         /* the last LMS steps of L */
         if (LMS > 0)
            Steps (gen, res, LDS + 1, r, LMS);
      }

      for (k = 0; k <= DeltaL; k += 2) {
         khi = gofs_Chi2 (res->H[k]->NbExp, res->H[k]->Count, res->H[k]->jmin,
            res->H[k]->jmax);
         statcoll_AddObs (res->H[k]->sVal1, khi);
         khi = gofs_Chi2 (res->M[k]->NbExp, res->M[k]->Count, res->M[k]->jmin,
            res->M[k]->jmax);
         statcoll_AddObs (res->M[k]->sVal1, khi);
         khi = gofs_Chi2 (res->R[k]->NbExp, res->R[k]->Count, res->R[k]->jmin,
            res->R[k]->jmax);
         statcoll_AddObs (res->R[k]->sVal1, khi);
         khi = gofs_Chi2 (res->J[k]->NbExp, res->J[k]->Count, res->J[k]->jmin,
            res->J[k]->jmax);
         statcoll_AddObs (res->J[k]->sVal1, khi);
         khi = gofs_Chi2 (res->C[k]->NbExp, res->C[k]->Count, res->C[k]->jmin,
            res->C[k]->jmax);
         statcoll_AddObs (res->C[k]->sVal1, khi);
         if (swrite_Counters)
            WriteDetailsWalk (res, k, n);
      }
   }

   for (k = 0; k <= DeltaL; k += 2) {
      for (m = 0; m < swalk_rw_N; m++) {
         switch (m) {
         case swalk_rwH:
            Q = res->H[k];
            break;
         case swalk_rwM:
            Q = res->M[k];
            break;
         case swalk_rwJ:
            Q = res->J[k];
            break;
         case swalk_rwR:
            Q = res->R[k];
            break;
         case swalk_rwC:
            Q = res->C[k];
            break;
         default:
            util_Error ("swalk_RandomWalk1:   no such case2");
         }
         V[0] = Q->degFree;
         Q->pVal1->NObs = Q->sVal1->NObs;
         gofw_ActiveTests2 (Q->sVal1->V, Q->pVal1->V, N, wdist_ChiSquare,
            V, Q->sVal2, Q->pVal2);
         sres_GetChi2SumStat (Q);
     }
   }

   if (swrite_Basic) {
      WriteResultWalk (res, N);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      swalk_DeleteRes (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataWalk1a (unif01_Gen *gen, char *TestName,
   long N, long n, int r, int s, int t, long L, bitset_BitSet maskc)
{
   int i;

   swrite_Head (gen, TestName, N, n, r);
   printf (",   s = %1d,  t =  %1d,   L = %1ld\n\n", s, t, L);
   printf ("   C = { ");

   util_Assert (t <= 31, "swalk_RandomWalk1a:   t > 31");
   for (i = 0; i < t; i++) {
      if (bitset_TestBit (maskc, i)) {
         printf ("%1d", i);
         if (i < t - 1)
            printf (", ");
      }
   }
   printf (" }\n\n\n");
}


/*-------------------------------------------------------------------------*/

void swalk_RandomWalk1a (unif01_Gen *gen, swalk_Res *res,
   long N, long n, int r, int s, int t, long L, bitset_BitSet maskc)
{
   swalk_rwType m;
   long z2, z1, y;
   long C, J, R, M, S_2, S, X;    /* Statistics */
   long i, j, pas, Rep, Seq;      /* Indices */
   bitset_BitSet ens;
   double khi;                    /* ChiSquare value */
   double V[1];                   /* Number deg. of freedom for ChiSquare */
   long NbClasses;
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "swalk_RandomWalk1a test";
   sres_Chi2 *Q;

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataWalk1a (gen, TestName, N, n, r, s, t, L, maskc);

   util_Assert (!(L & 1), "swalk_RandomWalk1a:   L is odd");
   util_Assert (L > 0, "swalk_RandomWalk1a:   L <= 0");
   util_Assert (r + s <= PREC, "swalk_RandomWalk1a:   r + s > 32");
   util_Assert (s <= PREC, "swalk_RandomWalk1a:   s > 32");
   if (n < 3.0 * gofs_MinExpected) {
      util_Warning (TRUE, "swalk_RandomWalk1a:   n < 3*gofs_MinExpected");
      return;
   }
   if (res == NULL) {
      localRes = TRUE;
      res = swalk_CreateRes ();
   }
   InitRes (res, NULL, N, L, L, "swalk_RandomWalk1a");

   /* Compute the expected numbers */
   CalcNbExp (n, L, 0, res);

   /* Merge classes for the ChiSquare */
   for (m = 0; m < swalk_rw_N; m++) {
      switch (m) {
      case swalk_rwH:
         Q = res->H[0];
         break;
      case swalk_rwM:
         Q = res->M[0];
         break;
      case swalk_rwJ:
         Q = res->J[0];
         break;
      case swalk_rwR:
         Q = res->R[0];
         break;
      case swalk_rwC:
         Q = res->C[0];
         break;
      default:
         util_Error ("swalk_RandomWalk1a:   no such case");
      }
      if (swrite_Classes) {
         printf ("===============================================\nThe ");
         printf ("%s", swalk_rwName[m]);
         printf ("\n");
         gofs_WriteClasses (Q->NbExp, Q->Loc, Q->jmin, Q->jmax, 0);
      }
      gofs_MergeClasses (Q->NbExp, Q->Loc, &(Q->jmin), &(Q->jmax), &NbClasses);

      if (swrite_Classes)
         gofs_WriteClasses (Q->NbExp, Q->Loc, Q->jmin, Q->jmax, NbClasses);

      /* Set description for second level statistical collectors */
      sprintf (str, "The N statistic values (a ChiSquare with %ld degrees of"
                    " freedom) ", NbClasses - 1);
      statcoll_SetDesc (Q->sVal1, str);
      Q->degFree = NbClasses - 1;
   }

   /* Generate the first t bits */
   z1 = 0;
   for (i = 0; i <= (t - 1) / s; i++) {
      z2 = unif01_StripB (gen, r, s);
      for (j = 1; j <= s; j++) {
         z1 = 2 * z1 + (z2 & 1);
         z2 /= 2;
      }
   }
   j = 0;
   z2 = unif01_StripB (gen, r, s);

   /* Beginning of test */
   for (Seq = 1; Seq <= N; Seq++) {

      /* Reset counters to zero for each sequence */
      for (i = res->H[0]->jmin; i <= res->H[0]->jmax; i++)
         res->H[0]->Count[i] = 0;
      for (i = res->R[0]->jmin; i <= res->R[0]->jmax; i++)
         res->R[0]->Count[i] = 0;
      for (i = res->J[0]->jmin; i <= res->J[0]->jmax; i++)
         res->J[0]->Count[i] = 0;
      for (i = res->M[0]->jmin; i <= res->M[0]->jmax; i++)
         res->M[0]->Count[i] = 0;
      for (i = res->C[0]->jmin; i <= res->C[0]->jmax; i++)
         res->C[0]->Count[i] = 0;

      for (Rep = 1; Rep <= n; Rep++) {
         C = J = R = M = S_2 = S = 0;
         pas = 0;

         /* Generate a random walk of L steps */
         do {
            do {
               ++pas;
               ++j;
               z1 = 2 * z1 + (z2 & 1);
               z2 /= 2;
               ens = maskc & z1;
               y = 0;
               for (i = 0; i < t; i++) {
                  if (bitset_TestBit (ens, i))
                     ++y;
               }
               if (y & 1)
                  X = 1;
               else
                  X = -1;
               S += X;
               if (S > M)
                  M = S;
               if (S == 0)
                  ++R;
               if (pas & 1) {
                  if (S > 0)
                     ++J;
                  if (S * S_2 < 0)
                     ++C;
                  S_2 = S;
               }
            } while (!(j == s || pas == L));
            if (j == s) {
               j = 0;
               z2 = unif01_StripB (gen, r, s);
            }
         } while (pas != L);

         /* Update counters */
         ++res->H[0]->Count[res->H[0]->Loc[(L + S) / 2]];
         ++res->M[0]->Count[res->M[0]->Loc[M]];
         ++res->J[0]->Count[res->J[0]->Loc[2 * J]];
         ++res->R[0]->Count[res->R[0]->Loc[R]];
         ++res->C[0]->Count[res->C[0]->Loc[C]];
      }

      khi = gofs_Chi2 (res->H[0]->NbExp, res->H[0]->Count, res->H[0]->jmin,
         res->H[0]->jmax);
      statcoll_AddObs (res->H[0]->sVal1, khi);
      khi = gofs_Chi2 (res->M[0]->NbExp, res->M[0]->Count, res->M[0]->jmin,
         res->M[0]->jmax);
      statcoll_AddObs (res->M[0]->sVal1, khi);
      khi = gofs_Chi2 (res->R[0]->NbExp, res->R[0]->Count, res->R[0]->jmin,
         res->R[0]->jmax);
      statcoll_AddObs (res->R[0]->sVal1, khi);
      khi = gofs_Chi2 (res->J[0]->NbExp, res->J[0]->Count, res->J[0]->jmin,
         res->J[0]->jmax);
      statcoll_AddObs (res->J[0]->sVal1, khi);
      khi = gofs_Chi2 (res->C[0]->NbExp, res->C[0]->Count, res->C[0]->jmin,
         res->C[0]->jmax);
      statcoll_AddObs (res->C[0]->sVal1, khi);

      if (swrite_Counters)
         WriteDetailsWalk (res, 0, n);
   }

   for (m = 0; m < swalk_rw_N; m++) {
      switch (m) {
      case swalk_rwH:
         Q = res->H[0];
         break;
      case swalk_rwM:
         Q = res->M[0];
         break;
      case swalk_rwJ:
         Q = res->J[0];
         break;
      case swalk_rwR:
         Q = res->R[0];
         break;
      case swalk_rwC:
         Q = res->C[0];
         break;
      default:
         util_Error ("swalk_RandomWalk1a:   no such case2");
      }
      V[0] = Q->degFree;
      Q->pVal1->NObs = Q->sVal1->NObs;
      gofw_ActiveTests2 (Q->sVal1->V, Q->pVal1->V, N, wdist_ChiSquare, V,
                         Q->sVal2, Q->pVal2);
      sres_GetChi2SumStat (Q);
   }

   if (swrite_Basic) {
      WriteResultWalk (res, N);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      swalk_DeleteRes (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataGeo (unif01_Gen *gen, char *TestName, 
   long N, long n, int r, double Mu, swalk_AlgoType Algo)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   Mu = %10.8f,   Algo = ", Mu);
   if (Algo == swalk_AlgoP)
      printf ("AlgoP\n\n");
   else
      printf ("AlgoN\n\n");
   printf ("   Expected length of a walk = %14.2f\n\n\n", 1.0 / (1.0 - Mu));
}


/*-------------------------------------------------------------------------*/

static void WriteNbExpCount (sres_Chi2 *res, double Prob[])
/* 
 * Writes the expected numbers, the observed numbers, and the normalized
 * values in swalk_VarGeo.
 */
{
   long L;
   double Ecart;
   double y;

   printf ("--------------------------------------------------\n"
      "Length  NumExpected  NumObserved  Normalized value\n\n");
   for (L = res->jmin; L < res->jmax; L = res->Loc[L + 1]) {
      printf ("%4ld %14.2f %10ld ", L, res->NbExp[L], res->Count[L]);
      Ecart = sqrt (res->NbExp[L] * (1.0 - Prob[L]));
      y = (res->Count[L] - res->NbExp[L]) / Ecart;
      printf ("%14.2f\n", y);
   }
   L = res->jmax;
   printf ("%4ld %14.2f %10ld ", L, res->NbExp[L], res->Count[L]);
   Ecart = sqrt (res->NbExp[L] * (1.0 - Prob[L]));
   y = (res->Count[L] - res->NbExp[L]) / Ecart;
   printf ("%14.2f\n\n\n", y);
}


/*-------------------------------------------------------------------------*/

static void AlgorithmP (unif01_Gen *gen, sres_Chi2 *res, double Prob[],
   long N, long n, int r, double Mu)
{
   long j;
   long L;
   long Seq;
   double X;
   double U;

   for (Seq = 1; Seq <= N; Seq++) {
      for (L = res->jmin; L <= res->jmax; L++)
         res->Count[L] = 0;

      for (j = 1; j <= n; j++) {
         L = 1;
         U = unif01_StripD (gen, r);
         while (U < Mu) {
            ++L;
            U = unif01_StripD (gen, r);
         }
         if (L >= res->jmax)
            ++res->Count[res->Loc[res->jmax]];
         else
            ++res->Count[res->Loc[L]];
      }
      if (swrite_Counters)
         WriteNbExpCount (res, Prob);

      X = gofs_Chi2 (res->NbExp, res->Count, res->jmin, res->jmax);
      statcoll_AddObs (res->sVal1, X);
   }
}


/*-------------------------------------------------------------------------*/

static void AlgorithmN (unif01_Gen *gen, sres_Chi2 *res, double Prob[],
   long N, long n, int r, double Mu)
{
   long j;
   long L;
   long Seq;
   double X;
   double U;

   Mu = 1.0 - Mu;
   for (Seq = 1; Seq <= N; Seq++) {
      for (L = res->jmin; L <= res->jmax; L++)
         res->Count[L] = 0;

      for (j = 1; j <= n; j++) {
         L = 1;
         U = unif01_StripD (gen, r);
         while (U >= Mu) {
            ++L;
            U = unif01_StripD (gen, r);
         }
         if (L >= res->jmax)
            ++res->Count[res->Loc[res->jmax]];
         else
            ++res->Count[res->Loc[L]];
      }
      if (swrite_Counters)
         WriteNbExpCount (res, Prob);

      X = gofs_Chi2 (res->NbExp, res->Count, res->jmin, res->jmax);
      statcoll_AddObs (res->sVal1, X);
   }
}


/*-------------------------------------------------------------------------*/

static void swalk_VarGeo (unif01_Gen *gen, sres_Chi2 *res,
   long N, long n, int r, double Mu, swalk_AlgoType Algo)
{
   const double epsilon = 1.0E-10;
   long L;
   double nLR = n;
   double V[1];                /* Number degrees of freedom for ChiSquare */
   char str[LENGTH + 1];
   long tt;
   long NbClasses;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "swalk_VarGeo test";
   double *Prob;

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataGeo (gen, TestName, N, n, r, Mu, Algo);

   util_Assert (r < PREC, "swalk_VarGeo:   r > 52");
   util_Assert (Mu > 0.0 && Mu < 1.0, "swalk_VarGeo:   Mu not in (0,1)");

   /* We consider only the terms of the geometric law with */
   /* probability > epsilon */
   tt = 1 + (log (epsilon) - num2_log1p (-Mu)) / log (Mu);
   Prob = util_Calloc (1 + (size_t) tt, sizeof (double));

   /* The probabilities and the expected numbers: NbExp = n*Prob */
   Prob[1] = 1.0 - Mu;
   for (L = 1; L <= tt - 2; L++)
      Prob[L + 1] = Mu * Prob[L];
   Prob[tt] = fbar_Geometric (1.0 - Mu, tt);

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, tt, "swalk_VarGeo");

   for (L = 1; L <= tt; L++)
      res->NbExp[L] = nLR * Prob[L];

   res->jmin = 1;
   res->jmax = tt;
   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, res->Loc, res->jmin, res->jmax, 0);
   gofs_MergeClasses (res->NbExp, res->Loc, &res->jmin, &res->jmax,
                      &NbClasses);
   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, res->Loc, res->jmin, res->jmax,
                         NbClasses);

   sprintf (str, "The N statistic values (a ChiSquare with %1ld degrees"
                 " of freedom):", NbClasses - 1);
   statcoll_SetDesc (res->sVal1, str);
   res->degFree = NbClasses - 1;
   if (res->degFree < 1) {
      if (localRes)
         sres_DeleteChi2 (res);
      return;
   }

   if (Algo == swalk_AlgoP)
      AlgorithmP (gen, res, Prob, N, n, r, Mu);
   else
      AlgorithmN (gen, res, Prob, N, n, r, Mu);

   V[0] = res->degFree;
   res->pVal1->NObs = N;
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
   util_Free (Prob);
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
}


/*-------------------------------------------------------------------------*/

void swalk_VarGeoP (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r, double Mu)
{
   swalk_VarGeo (gen, res, N, n, r, Mu, swalk_AlgoP);
}

/*-------------------------------------------------------------------------*/

void swalk_VarGeoN (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r, double Mu)
{
   swalk_VarGeo (gen, res, N, n, r, Mu, swalk_AlgoN);
}
