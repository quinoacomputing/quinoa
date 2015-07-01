/*************************************************************************\
 *
 * Package:        TestU01
 * File:           scomp.c
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

#include "scomp.h"
#include "sres.h"
#include "swrite.h"
#include "unif01.h"

#include "fbar.h"
#include "wdist.h"
#include "gofw.h"
#include "gofs.h"
#include "statcoll.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>


#define LENGTH 100

/* Empirical Mean for Lempel-Ziv test for 2^3 <= n <= 2^28. Obtained by
   simulation with N = 1000  */
static const double LZMu[] = {
   0.0,   0.0,   0.0,   4.44,   7.64,
   12.5,   20.8,   34.8,   58.9,   101.1,
   176.0,   310.0,   551.9,   992.3,   1799.,
   3286.2,   6041.5,   11171.5,   20761.8,   38760.4,
   72654.,   136677.,   257949.,   488257.,   926658.,
   1762965.,   3361490.,   6422497.,   12293930.
};

/* Empirical Standard Deviation for Lempel-Ziv test for 2^3 <= n <= 2^28 */
static const double LZSigma[] = {
   0.0,   0.0,   0.0,   0.49,   0.51,
   0.62,   0.75,   0.78,   0.86,   0.94,
   1.03,   1.19,   1.43,   1.68,   2.09,
   2.46,   3.36,   4.2,   5.4,   6.8,
   9.1,   10.9,   14.7,   19.1,   25.2,
   33.5,   44.546,   58.194,   75.513
};



/*--------------------------------- Types ---------------------------------*/

/* Bit trie used in Lempel-Ziv test. If left != NULL, this means a 0 bit.
   If right != NULL, this means a 1 bit. The word is the sequence obtained
   by following the tree until a NULL pointer is met. */
   
struct BitTrie_t {
   struct BitTrie_t *left;
   struct BitTrie_t *right;
};
typedef struct BitTrie_t BitTrie_t;




/*-------------------------------- Functions ------------------------------*/

static void DeleteBitTrie (BitTrie_t *tree)
{
   if (tree == NULL)
      return;
   DeleteBitTrie (tree->left);
   DeleteBitTrie (tree->right);
   util_Free (tree);
}


/*=========================================================================*/


static void InitRes (
   scomp_Res *res,            /* Results holder */
   long N,                    /* Number of replications */
   int jmax,                  /* Max class index for size of jumps */
   int tmax                   /* Max class index for linear complexity */
)
/* 
 * Initializes the scomp_Res structure
 */
{
   sres_InitBasic (res->JumpNum, N,
      "scomp_LinearComp:   Number of Jumps");
   sres_InitChi2 (res->JumpSize, N, jmax,
      "scomp_LinearComp:   Size of Jumps");
   sres_InitChi2 (res->LinComp, N, tmax,
      "scomp_LinearComp:   Linear Complexity");
}


/*-------------------------------------------------------------------------*/

scomp_Res * scomp_CreateRes (void)
{
   scomp_Res *res;
   res = util_Malloc (sizeof (scomp_Res));
   res->JumpNum = sres_CreateBasic ();
   res->JumpSize = sres_CreateChi2 ();
   res->LinComp = sres_CreateChi2 ();
   return res;
}


/*-------------------------------------------------------------------------*/

void scomp_DeleteRes (scomp_Res *res)
{
   if (res == NULL)
      return;
   sres_DeleteBasic (res->JumpNum);
   sres_DeleteChi2 (res->JumpSize);
   sres_DeleteChi2 (res->LinComp);
   util_Free (res);
}


/*=========================================================================*/

static void WriteDataJumps (unif01_Gen *gen, char *TestName, long N, long n,
   int r, int s, double muComp, double mu, double sigma)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",    s = %1d\n", s);
   if (swrite_Parameters) {
      printf ("\n      muComp = ");
      num_WriteD (muComp, 12, 4, 2);
      printf ("\n      Mu     = ");
      num_WriteD (mu, 12, 4, 2);
      printf ("\n      Sigma  = ");
      num_WriteD (sigma, 12, 4, 2);
   }
   printf ("\n\n");
}


/*-------------------------------------------------------------------------*/

static void BerlekampMassey (
   scomp_Res *res,
   long n,                  /* Number of bits */
   double *pComp,           /* Linear complexity */
   double *pNumJ,           /* Number of jumps */
   int *Bits,
   int *Polyb,
   int *Polyc,
   int *PolycOld
   )
/*
 * Berlekamp-Massey algorithm to calculate the linear complexity.
 */
{
   int b;
   long Loc;
   long i;
   long m;
   long k;
   long L;                  /* Linear complexity */
   long NumJ;               /* Number of jumps */
   sres_Chi2 *resl = res->JumpSize;

   for (k = 0; k <= resl->jmax; k++)
      resl->Count[k] = 0;

   Polyc[0] = 1;
   Polyb[0] = 1;
   L = 0;
   NumJ = 0;
   k = 0;
   m = -1;

   while (k < n) {
      /* Return the value of the current polynomial to see if it can
         generate the next bit */
      b = 0;
      for (i = 1; i <= L; i++)
         /* b ^= Polyc[i] * Bits[k + 1 - i]; */
         b = (b + Polyc[i] * Bits[k + 1 - i]) & 1;

      if (Bits[k + 1] != b) {
         /* Update c(x)_old and c(x) */
         for (i = 0; i <= L; i++)
            PolycOld[i] = Polyc[i];

         for (i = 0; i <= L; i++) {
            if (Polyb[i] == 1)
               Polyc[k - m + i] = (++Polyc[k - m + i]) & 1;
         }
         if (2 * L <= k) {
            L = k + 1 - L;
            NumJ++;
            Loc = labs (k + 1 - 2 * L);
            if (Loc <= resl->jmax)
               ++resl->Count[Loc];
            else
               ++resl->Count[resl->jmax];
            /* Update B */
            for (i = 0; i <= L; i++)
               Polyb[i] = PolycOld[i];
            m = k;
         }
      }
      ++k;
   }
   *pNumJ = NumJ;
   *pComp = L;
}


/*=========================================================================*/

void scomp_LinearComp (unif01_Gen *gen, scomp_Res *res,
                       long N, long n, int r, int s)
{
   const double epsilon = 1.0E-10;
   const int tt = 1 - 0.5 * num_Log2 (3 * epsilon); /* Dimension */
   const long K0 = n/s;
   long i, Seq;
   int j, k;
   int M0;
   double NumJ;                   /* Number of Jumps */
   double sigma, mu;              /* Parameters of number of jumps */
   double comp;                   /* Linear complexity */
   double muComp;                 /* Mean of linear complexity */
   unsigned long Nombre;          /* Random number */
   int Parite;
   double *Prob;
   long *Loca;
   long tmin, tmax, NbClasses;
   double X2;
   double temp;
   int *Bits;                     /* 4 Arrays of bits */
   int *Polyb;
   int *Polyc;
   int *PolycOld;
   double Param[1];
   char str[LENGTH + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "scomp_LinearComp test";
   sres_Basic *resJN;
   sres_Chi2 *resJL;
   sres_Chi2 *resLC;
   lebool JL_OK = TRUE;          /* If TRUE do the JL test, otherwise not */
   lebool LC_OK = FALSE;         /* If TRUE do the LC test, otherwise not */

   Timer = chrono_Create ();
   n = K0 * s;
   Parite = n & 1;
   if (n >= DBL_MAX_EXP)
      temp = 0.0;
   else
      temp = pow (2.0, -(double) n);

   mu = n / 4.0 + (4 + Parite) / 12.0 - temp / 3.0;
   sigma = n / 8.0 - (2 - Parite)/(9.0 - Parite) + n * temp / 6.0
           + (6 + Parite) * temp / 18.0 - temp * temp / 9.0;
   sigma = sqrt (sigma);
   muComp = n / 2.0 + (4 + Parite) / 18.0;
   M0 = num_Log2 (mu / gofs_MinExpected);
   if (M0 < 2) {
      /* 0 degree of freedom for the chi2, do not do the test JL. */
      JL_OK = FALSE;
   }

   if (swrite_Basic)
      WriteDataJumps (gen, TestName, N, n, r, s, muComp, mu, sigma);

   /* util_Assert (M0 > 1, "scomp_LinearComp:   n*s is too small"); */
   util_Assert (M0 <= num_MaxTwoExp, "scomp_LinearComp:   M0 > num_MaxTwoExp");

   Prob = util_Calloc (1 + (size_t) tt, sizeof (double));
   Bits = util_Calloc ((size_t) n + 1, sizeof (int));
   Polyb = util_Calloc ((size_t) n + 1, sizeof (int));
   Polyc = util_Calloc ((size_t) n + 1, sizeof (int));
   PolycOld = util_Calloc ((size_t) n + 1, sizeof (int));

   if (res == NULL) {
      localRes = TRUE;
      res = scomp_CreateRes ();
   }
   M0 = util_Max (M0, 1);
   InitRes (res, N, M0, tt);
   resJN = res->JumpNum;
   resJL = res->JumpSize;
   resLC = res->LinComp;
   Loca = resLC->Loc;

   if (N > 2.0 * gofs_MinExpected) {
      /* Compute the expected probabilities for the linear complexity. */
      /* We put in Prob[k] the probabilities for i = k and i = -k of   */
      /* the statistic defined in the NIST document 800-22, p. 86. */
      temp = Prob[0] = 0.5;
      for (k = 1; k < tt; k++) {
	 Prob[k] = 1.5 * pow (4.0, -(double) k);
	 temp += Prob[k];
      }
      Prob[tt] = 1.0 - temp;
      for (k = 0; k <= tt; k++) {
	 resLC->Count[k] = 0;
	 resLC->NbExp[k] = N * Prob[k];
      }
      tmin = 0;
      tmax = tt;
      if (swrite_Classes) {
	 printf ("Classes for the linear complexity:\n");
	 gofs_WriteClasses (resLC->NbExp, Loca, tmin, tmax, 0);
      }
      gofs_MergeClasses (resLC->NbExp, Loca, &tmin, &tmax, &NbClasses);
      resLC->jmax = tmax;
      resLC->jmin = tmin;
      resLC->degFree = NbClasses - 1;
      if (NbClasses < 2) {
	 /* 0 degree of freedom for the chi2, do not do the test LC. */
	 LC_OK = FALSE;
      } else
	 LC_OK = TRUE;
   }

   statcoll_SetDesc (resJN->sVal1,
      "The number of jumps: the N statistic values (a standard normal):");
   sprintf (str, "The jumps size: the N statistic values (a ChiSquare"
                 " with %1d degrees of freedom):", M0 - 1);
   statcoll_SetDesc (resJL->sVal1, str);

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 0; i < K0; i++) {
         Nombre = unif01_StripB (gen, r, s);
         for (j = s; j >= 1; j--) {
            Bits[s * i + j] = Nombre & 1;
            Nombre >>= 1;
         }
      }

      BerlekampMassey (res, n, &comp, &NumJ, Bits, Polyb, Polyc, PolycOld);

      /* Value of the statistic for the linear complexity */
      if (LC_OK) {
	 comp = comp - muComp;
	 if (Parite)
	    comp = -comp;
	 comp += 2.0 / 9.0;
         /* comp is now an integer: truncate correctly and avoid off-by-1
            error because of small floating-point inaccuracies. */
	 if (comp >= 0.0)
	    k = comp + 0.5;
	 else
	    k = comp - 0.5;
	 if (k < 0)
	    k = -k;
	 if (k >= tt)
	    ++resLC->Count[Loca[tt]];
	 else
	    ++resLC->Count[Loca[k]];
      }

      /* Value of the normal statistic for the number of jumps */
      statcoll_AddObs (resJN->sVal1, (NumJ - mu) / sigma);

      /* Value of the statistic for the size of the jumps */
      if (JL_OK) {
	 for (k = 1; k < M0; k++) {
	    resJL->NbExp[k] = NumJ / num_TwoExp[k];
	    resJL->Loc[k] = k;
	 }
	 resJL->NbExp[M0] = NumJ / num_TwoExp[M0 - 1];
	 resJL->Loc[M0] = M0;
	 resJL->jmax = M0;
	 resJL->jmin = 1;
	 resJL->degFree = M0 - 1;

	 X2 = gofs_Chi2 (resJL->NbExp, resJL->Count, 1, M0);
	 statcoll_AddObs (resJL->sVal1, X2);
	 if (swrite_Classes) {
	    printf ("\n\nClasses for the size of the jumps:\n");
	    gofs_WriteClasses (resJL->NbExp, (long *) NULL, 1, M0, 0);
	 }
         if (swrite_Counters)
            tables_WriteTabL (resJL->Count, 1, M0, 5, 10,
                "Size of the jumps:   observed numbers");
      }
   }

   gofw_ActiveTests2 (resJN->sVal1->V, resJN->pVal1->V, N, wdist_Normal,
      (double *) NULL, resJN->sVal2, resJN->pVal2);
   resJN->pVal1->NObs = N;
   sres_GetNormalSumStat (resJN);

   if (JL_OK) {
      Param[0] = M0 - 1;
      gofw_ActiveTests2 (resJL->sVal1->V, resJL->pVal1->V, N, 
         wdist_ChiSquare, Param, resJL->sVal2, resJL->pVal2);
      resJL->pVal1->NObs = N;
      sres_GetChi2SumStat (resJL);
   }

   if (LC_OK) {
      X2 = gofs_Chi2 (resLC->NbExp, resLC->Count, tmin, tmax);
      resLC->sVal2[gofw_Mean] = X2;
      resLC->pVal2[gofw_Mean] = fbar_ChiSquare2 (NbClasses - 1, 8, X2);
   }

   if (swrite_Basic) {
      if (JL_OK) {
	 printf ("\n-----------------------------------------------\n");
	 if (N == 1) {
            printf ("Number of degrees of freedom          : %4ld\n",
                    resJL->degFree);
	    printf ("Chi2 statistic for size of jumps      :");
	    gofw_Writep2 (resJL->sVal2[gofw_Mean], resJL->pVal2[gofw_Mean]);
	 } else {
	    printf ("Test results for the size of jumps:\n");
	    gofw_WriteActiveTests0 (N, resJL->sVal2, resJL->pVal2);
            swrite_Chi2SumTest (N, resJL);
	 }
	 if (swrite_Collectors)
	    statcoll_Write (resJL->sVal1, 5, 14, 4, 3);
      }

      printf ("\n-----------------------------------------------\n");
      if (N == 1) {
         printf ("Normal statistic for number of jumps  :");
         gofw_Writep2 (resJN->sVal2[gofw_Mean], resJN->pVal2[gofw_Mean]);
      } else {
         printf ("Test results for the number of jumps:\n");
         gofw_WriteActiveTests0 (N, resJN->sVal2, resJN->pVal2);
         swrite_NormalSumTest (N, resJN);
      }
      if (swrite_Collectors)
         statcoll_Write (resJN->sVal1, 5, 14, 4, 3);

      if (LC_OK) {
         printf ("\n-----------------------------------------------\n");
         printf ("Test results for the linear complexity:\n\n");
         printf ("Number of degrees of freedom          : %4ld\n",
                 resLC->degFree);
         printf ("Chi2 statistic on the N replications  :");
         gofw_Writep2 (resLC->sVal2[gofw_Mean], resLC->pVal2[gofw_Mean]);
         if (swrite_Classes)
	    gofs_WriteClasses (resLC->NbExp, Loca, tmin, tmax, NbClasses);
         if (swrite_Counters)
            tables_WriteTabL (resLC->Count, tmin, tmax, 5, 10,
               "Linear Complexity:   observed numbers");
      }

      printf ("\n\n");
      swrite_Final (gen, Timer);
   }

   util_Free (Prob);
   util_Free (Bits);
   util_Free (Polyb);
   util_Free (Polyc);
   util_Free (PolycOld);
   if (localRes)
      scomp_DeleteRes (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataLZ (
   unif01_Gen *gen,      /* generator */
   char *Test,           /* Test name */
   long N,               /* Number of replications */
   int k,                /* Sample size n = 2^k */
   int r,                /* r first bits of each random number dropped */
   int s                 /* s bits of each random number used */
)
{
   long n;
   n = num_TwoExp[k];
   swrite_Head (gen, Test, N, n, r);
   printf (",   s = %4d,   k = %4d\n\n", s, k);
}


/*-------------------------------------------------------------------------*/

static long LZ78 (unif01_Gen * gen, long n, int r, int s)
/*
 * The parameters are the same as in scomp_LempelZiv. The trie contains
 * a left (right) branch if a word with a 0 (1) bit after the prefix has
 * been seen before. We descend one level in the trie with each bit until
 * a leaf is met. Add a branch for a new word, and restart at root.
 */
{
   const unsigned long kMAX = 1UL << (s - 1);
   unsigned long Y, k;
   long i;                        /* Count the number of bits overall */
   long W;                        /* Count the number of words */
   lebool done = FALSE;          /* Start a new word */
   BitTrie_t *trie, *root;

   W = i = 0;
   trie = root = util_Malloc (sizeof (BitTrie_t));
   trie->left = trie->right = NULL;
   Y = unif01_StripB (gen, r, s);
   k = kMAX;

   while (i < n) {
      /* Start a new word: match it as far as possible in the trie */
      done = FALSE;
      trie = root;
      while (!done) {
         if ((Y & k) == 0) {      /* Bit 0 */
            if (trie->left) {
	       /* We have seen it before: descend in branch */
               trie = trie->left;
            } else {
	       /* A leaf: this is a new word */
               W++;
               done = TRUE;
               trie->left = util_Malloc (sizeof (BitTrie_t));
               trie = trie->left;
               trie->left = trie->right = NULL;
            }

         } else {                 /* Bit 1 */
            if (trie->right) {
               trie = trie->right;
            } else {
               W++;
               done = TRUE;
               trie->right = util_Malloc (sizeof (BitTrie_t));
               trie = trie->right;
               trie->left = trie->right = NULL;
            }
         }
         i++;
         if (i >= n) {
            done = TRUE;
            if ((trie->left != NULL) || (trie->right != NULL))
               W++;
            break;
	 }
         k >>= 1;
         if (k == 0) {
	    /* Have used the s bits in the number; generate a new number */
            Y = unif01_StripB (gen, r, s);
            k = kMAX;
         }
      }
   }
   DeleteBitTrie (root);
   return W;
}


/*-------------------------------------------------------------------------*/

void scomp_LempelZiv (unif01_Gen *gen, sres_Basic *res,
   long N, int t, int r, int s)
{
   long Seq, n;
   double X;
   /*   const double lg_n = num_Log2 ((double) n); */
   long W;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "scomp_LempelZiv test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataLZ (gen, TestName, N, t, r, s);
   util_Assert (r + s <= 32, "scomp_LempelZiv:   r + s > 32");
   util_Assert (t <= 28, "scomp_LempelZiv:   k > 28");
   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateBasic ();
   }
   n = num_TwoExp[t];
   sres_InitBasic (res, N, "scomp_LempelZiv");
   statcoll_SetDesc (res->sVal1, "sVal1:   a standard normal");

   for (Seq = 1; Seq <= N; Seq++) {
      W = LZ78 (gen, n, r, s);
      /*  X = (W - n / lg_n) / sqrt (0.266 * n / (lg_n * lg_n * lg_n)); */
      X = (W - LZMu[t]) / LZSigma[t];
      statcoll_AddObs (res->sVal1, X);
      if (swrite_Counters) {
         printf ("%12ld ", W);
         if (Seq % 5 == 0)
            printf ("\n");
         if (Seq >= N)
            printf ("\n\n");
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
   if (localRes)
      sres_DeleteBasic (res);
   chrono_Delete (Timer);
}
