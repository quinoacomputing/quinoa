/*************************************************************************\
 *
 * Package:        TestU01
 * File:           sstring.c
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
#include "bitset.h"

#include "sstring.h"
#include "unif01.h"
#include "wdist.h"
#include "swrite.h"
#include "sres.h"

#include "gofs.h"
#include "gofw.h"
#include "fbar.h"
#include "statcoll.h"

#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>





/*------------------------------ Constants --------------------------------*/

/* Minimal length (number of bits) of a sequence for LongestHeadRun */
#define LMIN 1000

/* Max string length for the correlations in PeriodsInStrings */
#define MAX_CORR 31

/* Max dimension of arrays */
#define DIM 1000

/* Max string lengths */
#define LEN1 200
#define LEN2 200





/*-------------------------------- Types ----------------------------------*/

typedef struct InfoListC *ListC;     /* A correlation list */

struct InfoListC {
   long Nb;                          /* Number of bits of a correlation in
                                        the initial computations; then
                                        number of occurences */
   bitset_BitSet C;                  /* A correlation */
   long Pop;                         /* Population related to C (and c) */
   ListC Ext;                        /* The smallest extension of C */
   ListC Ext0;                       /* The smallest extension of D longer
                                        than C, if C is an extension of D.
                                        Initially NULL */
   ListC Next;                       /* Next correlation of same length */
};

/* Corr contains the lists of correlations of each length for s in [0..smax]
 */
typedef struct {
   ListC Corr[MAX_CORR + 1];
   int smax;
} sstring_Corr;


/*----------------------------- Variables --------------------------------*/

lebool sstring_CorrFlag = FALSE;
lebool sstring_Counters = FALSE;




/*----------------------------- Functions --------------------------------*/

static void InitRes3 (
   sstring_Res3 *res,         /* Results holder */
   long N,                    /* Number of replications */
   int jmax                   /* Max class index for chi2 */
)
/* 
 * Initializes the sstring_Res3 structure
 */
{
   sres_InitBasic (res->NBits, N, "sstring_Run:   Number of Bits");
   sres_InitChi2 (res->NRuns, N, jmax, "sstring_Run:   Number of Runs");
   res->Count0 = util_Realloc (res->Count0, (jmax + 1) * sizeof (long));
   res->Count1 = util_Realloc (res->Count1, (jmax + 1) * sizeof (long));
   res->NRuns->jmin = 1;
   res->NRuns->degFree = jmax - 1;
}


/*-------------------------------------------------------------------------*/

sstring_Res3 * sstring_CreateRes3 (void)
{
   sstring_Res3 *res;
   res = util_Malloc (sizeof (sstring_Res3));
   res->NBits = sres_CreateBasic ();
   res->NRuns = sres_CreateChi2 ();
   res->Count0 = util_Calloc (1, sizeof (long));
   res->Count1 = util_Calloc (1, sizeof (long));
   return res;
}


/*-------------------------------------------------------------------------*/

void sstring_DeleteRes3 (sstring_Res3 *res)
{
   if (res == NULL)
      return;
   res->Count0 = util_Free (res->Count0);
   res->Count1 = util_Free (res->Count1);
   sres_DeleteBasic (res->NBits);
   sres_DeleteChi2 (res->NRuns);
   util_Free (res);
}


/*=========================================================================*/

static void InitRes2 (
   sstring_Res2 *res,         /* Results holder */
   long N,                    /* Number of replications */
   int jhigh                  /* Max class index for chi2 */
)
/* 
 * Initializes the sstring_Res2 structure
 */
{
   sres_InitDisc (res->Disc, N,
      "sstring_LongestHeadRun:   Global longest run of 1's");
   sres_InitChi2 (res->Chi, N, jhigh,
      "sstring_LongestHeadRun:   Block longest runs of 1's");
}


/*-------------------------------------------------------------------------*/

sstring_Res2 * sstring_CreateRes2 (void)
{
   sstring_Res2 *res;
   res = util_Malloc (sizeof (sstring_Res2));
   res->Chi = sres_CreateChi2 ();
   res->Disc = sres_CreateDisc ();
   return res;
}


/*-------------------------------------------------------------------------*/

void sstring_DeleteRes2 (sstring_Res2 *res)
{
   if (res == NULL)
      return;
   sres_DeleteChi2 (res->Chi);
   sres_DeleteDisc (res->Disc);
   util_Free (res);
}


/*=========================================================================*/

static void InitRes (
   sstring_Res *res,          /* Results holder */
   long N,                    /* Number of replications */
   int L,                     /* Size of blocks (number of bits) */
   int d,                     /* Parameter for sub-matrices */
   char *nam
)
/* 
 * Initializes res
 */
{
   int i;
   sres_InitBasic (res->Bas, N, nam);

   if (res->L > 0) {
      tables_DeleteMatrixL (&res->Counters);
      tables_DeleteMatrixD (&res->ZCounters);
   }
   res->Counters = tables_CreateMatrixL (L + 2, L + 1);
   res->ZCounters = tables_CreateMatrixD (L + 2, L + 1);

   if (d < 0)
      d = 0;
   for (i = d + 1; i <= res->d; i++)
      sres_DeleteBasic (res->Block[i]);

   for (i = res->d + 1; i <= d; i++)
      res->Block[i] = sres_CreateBasic ();

   for (i = 1; i <= d; i++)
      sres_InitBasic (res->Block[i], N, nam);

   res->L = L;
   res->d = d;
}


/*-------------------------------------------------------------------------*/

sstring_Res * sstring_CreateRes (void)
{
   sstring_Res *res;
   res = util_Malloc (sizeof (sstring_Res));
   memset (res, 0, sizeof (sstring_Res));
   res->Bas = sres_CreateBasic ();
   res->Style = tables_Plain;
   res->L = -1;
   res->d = 0;
   return res;
}


/*-------------------------------------------------------------------------*/

void sstring_DeleteRes (sstring_Res *res)
{
   if (res == NULL)
      return;

   if (res->d > 0) {
      int i;
      for (i = 1; i <= res->d; i++) {
         sres_DeleteBasic (res->Block[i]);
      }
   }
   if (res->L > 0) {
      tables_DeleteMatrixD (&res->ZCounters);
      tables_DeleteMatrixL (&res->Counters);
   }
   sres_DeleteBasic (res->Bas);
   util_Free (res);
}


/*=========================================================================*/

static long Psi (
   bitset_BitSet C,      /* Correlation making up tail of correlation k */
   long j,               /* Length of correlation C */
   long k                /* Correlation made up of 1 followed by 0's 
                            until C, i.e. k = 100...000C */
)
{
   /* j <=> c & k <=> k in the article */
   if (k > j)
      return 0;
   if (k <= 0)
      return (long) num_TwoExp[-k];
   if (bitset_TestBit (C, j - k))
      return 1;
   else
      return 0;
}


/*-------------------------------------------------------------------------*/

static void DeleteCorr (sstring_Corr *corr)
/*
 * Delete all correlations.
 */
{
   ListC Ci, OldCi;
   int i;

   if (corr == NULL)
      return;
   for (i = 0; i <= corr->smax; i++) {
      Ci = corr->Corr[i];
      while (Ci) {
         OldCi = Ci;
         Ci = Ci->Next;
         util_Free (OldCi);
      }
   }
   util_Free (corr);
}


/*-------------------------------------------------------------------------*/

static sstring_Corr * CreateCorr (int s)
/*
 * Compute all possible correlations for strings of length s
 */
{
   ListC CjE, Cj, Ci, OldCi, XS;
   sstring_Corr *corr;
   int j, i, k, Tmax;
   long p;

   corr = util_Malloc (sizeof (sstring_Corr));
   memset (corr, 0, sizeof (sstring_Corr));

   corr->smax = s;

   XS = corr->Corr[0] = util_Malloc (sizeof (struct InfoListC));
   XS->Nb = 0;
   XS->Pop = 1;
   XS->Ext = NULL;
   XS->Ext0 = NULL;
   XS->Next = NULL;

   XS = corr->Corr[1] = util_Malloc (sizeof (struct InfoListC));
   XS->Nb = 1;
   bitset_SetBit (XS->C, 0);
   XS->Pop = 2;
   XS->Ext = NULL;
   XS->Ext0 = NULL;
   XS->Next = NULL;

   for (i = 2; i <= s; i++) {            /* i is the string length */
      /* Count and build the list of correlations of length i. */
      Ci = corr->Corr[i] = util_Malloc (sizeof (struct InfoListC));

      for (j = 0; j <= i - 2; j++) {
         /* j is the length of correlation C in 100...00C */
         Cj = corr->Corr[j];

	 while (Cj) {
	    /* Compute the number of strings of length i and with corre- */
	    /* lation "10...0C", and if > 0, add this corr. to Corr[i] */
	    p = Cj->Pop * Psi (Cj->C, j, 2*j - i);
	    /* Check if 1C may be a correlation. Possible only if
	       C = 111....1, i.e. the last of the list */
	    if (Cj->Next == NULL)
	       p -= 2 * Psi (Cj->C, j, 2*j + 2 - i);
	    CjE = Cj->Ext;
	    /* Note: j <= j-2, i.e.  j+1 <= (i+j) / 2  */
	    Tmax = (i + j) / 2;
	    while (CjE && CjE->Nb <= Tmax) {
	       p -= CjE->Pop * Psi (Cj->C, j, 2*CjE->Nb - i);
	       CjE = CjE->Ext0;
	    }
	    /* p = number of strings looked for */
	    if (p > 0) {
	       /* Put this correlation in Ci */
	       Ci->Nb = i;
	       Ci->Pop = p;
	       Ci->Ext = NULL;
	       Ci->Ext0 = NULL;
	       /* Ci->C becomes Cj->C shifted right by i-j  */
	       /* positions, with a 1 in first position.    */
	       Ci->C = 0;
	       bitset_SetBit (Ci->C, 0);
	       if (j > 0) {
		  for (k = 0; k < j; k++) {
		     if (bitset_TestBit (Cj->C, k)) {
			bitset_SetBit (Ci->C, k + i - j);
		     }
		  }
	       }
	       if (Cj->Ext == NULL)
		  Cj->Ext = Ci;
	       else {
		  CjE = Cj->Ext;
		  while (CjE->Ext0)
		     CjE = CjE->Ext0;
		  CjE->Ext0 = Ci;
	       }
	       OldCi = Ci;
	       Ci = util_Malloc (sizeof (struct InfoListC));
	       OldCi->Next = Ci;
	    }
	    Cj = Cj->Next;
	 }
      }
      /* For j = i-1, we have the correlation "11...1" */
      Ci->C = 0;
      for (k = 0; k < i; k++) {
         bitset_SetBit (Ci->C, k);
      }
      Ci->Nb = i;
      Ci->Pop = 2;
      Ci->Ext = NULL;
      Ci->Ext0 = NULL;
      Ci->Next = NULL;
   }
   return corr;
}


/*=========================================================================*/

static void sstring_WriteCorr (sstring_Corr *corr, int s)
{
   ListC Cs;
   int k;
   char str [LEN1 + 1];

   if (corr == NULL) {
      util_Warning (TRUE,
         "sstring_WriteCorr:   corr is a NULL pointer");
      return;
   }

   if (corr->smax < s) {
      sprintf (str, "sstring_WriteCorr:   invalid s = %d", s);
      util_Error (str);
   }

   Cs = corr->Corr[s];
   if (Cs == NULL)
      return;
   printf ("\n-----------------------------------------------------\n"
      "List of correlations of length %d and their population\n\n", s);
   while (Cs) {
      for (k = 0; k < s; k++) {
         if (bitset_TestBit (Cs->C, k))
            printf ("1");
         else
            printf ("0");
      }
      printf ("%12ld\n", Cs->Pop);
      Cs = Cs->Next;
   }
   printf ("\n\n");
}


/*-------------------------------------------------------------------------*/

static bitset_BitSet GenerateC (
   unif01_Gen *gen,      /* Generator */
   int r,                /* Drop first r bits of each random number */
   int s                 /* Keep next s bits of each random number */
)
/*
 * Generate a string of s bits and return its correlation. To determine the
 * correlation of a bit string, compare strings g and d (initially equal).
 * For each iteration:
 *     1) drop leftmost bit of string g
 *     2) drop rightmost bit of string d
 *     3) if g = d, bit k of the correlation c is 1, otherwise 0.
 */
{
   int k;
   unsigned long g, d, lbit = s - 1;
   bitset_BitSet c = 0;

   /* Generate a random number; drop r most significant bits; keep s next
      bits */
   g = d = unif01_StripB (gen, r, s);

   /* Initialization of correlation, trivial case */
   bitset_SetBit (c, 0);
   for (k = 1; k < s; k++) {
      /* drop leftmost bit of string g */
      bitset_ClearBit (g, lbit);
      /* drop rightmost bit of string d	*/
      d >>= 1;
      /* if g = d, bit k of the correlation c is 1, otherwise 0 */
      if (g == d)
         bitset_SetBit (c, k);
      lbit--;
   }
   return c;
}


/*=========================================================================*/

static void WriteDataPeriod (
   unif01_Gen *gen,      /* generator */
   char *Test,           /* Test name */
   long N,               /* Number of replications */
   long n,               /* Sample size */
   int r,                /* r first bits of each random number dropped */
   int s                 /* s bits of each random number used */
)
{
   swrite_Head (gen, Test, N, n, r);
   printf (",   s = %4d\n\n", s);
}


/*-------------------------------------------------------------------------*/

void sstring_PeriodsInStrings (unif01_Gen *gen, sres_Chi2 *res,
   long N, long n, int r, int s)
{
   ListC XS, Ci;
   sstring_Corr *corr;
   long jhigh,                     /* Highest class for ChiSquare */
        jlow,                      /* Lowest class for ChiSquare */
        NbGroups;                  /* Number of classes for ChiSquare */
   long j, i;
   long Seq;                       /* One replication of the test */
   double Fraction, X2;
   bitset_BitSet D;                /* A correlation */
   double V[1];                    /* Number of ChiSquare degrees of freedom */
   char str [LEN1 + 1];
   double NbExp [DIM + 1];
   long Loca [DIM + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sstring_PeriodsInStrings test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataPeriod (gen, TestName, N, n, r, s);

   util_Assert (r >= 0, "sstring_PeriodsInStrings:   r < 0");
   util_Assert (r <= 31, "sstring_PeriodsInStrings:   r > 31");
   util_Assert (r + s <= 31, "sstring_PeriodsInStrings:   r + s > 31");
   util_Assert (s <= 31, "sstring_PeriodsInStrings:   s > 31");
   util_Assert (s >= 2, "sstring_PeriodsInStrings:   s < 2");
   /*   util_Assert (n > 2.0 * gofs_MinExpected,
	"sstring_PeriodsInStrings:    n <= 2*gofs_MinExpected"); */

   Fraction = n / num_TwoExp[s];
   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   corr = CreateCorr (s);
   if (sstring_CorrFlag)
      sstring_WriteCorr (corr, s);

   /* Get the expected numbers of the population count */
   XS = Ci = corr->Corr[s];
   j = 1;
   while (Ci) {
      NbExp[j] = Fraction * Ci->Pop;
      Ci = Ci->Next;
      ++j;
      util_Assert (j <= DIM, "sstring_PeriodsInStrings:   DIM too small");
   }
   jlow = 1;
   jhigh = j - 1;

   if (swrite_Classes)
      gofs_WriteClasses (NbExp, Loca, jlow, jhigh, 0);

   /* Merge classes for the chi-square test */
   gofs_MergeClasses (NbExp, Loca, &jlow, &jhigh, &NbGroups);

   if (swrite_Classes)
      gofs_WriteClasses (NbExp, Loca, jlow, jhigh, NbGroups);

   res->degFree = NbGroups - 1;
   if (res->degFree < 1) {
      if (localRes)
         sres_DeleteChi2 (res);
      return;
   }
   sres_InitChi2 (res, N, jhigh, "sstring_PeriodsInStrings");
   res->jmin = jlow;
   tables_CopyTabD (NbExp, res->NbExp, jlow, jhigh);
   tables_CopyTabL (Loca, res->Loc, jlow, jhigh);

   sprintf (str, "The N statistic values (a ChiSquare with %1ld degrees"
                 " of freedom):", NbGroups - 1);
   statcoll_SetDesc (res->sVal1, str);

   /* Test begins */
   for (Seq = 1; Seq <= N; Seq++) {
      /* Zero the population counters */
      Ci = XS;
      while (Ci) {
         Ci->Nb = 0;
         Ci = Ci->Next;
      }
      for (i = 1; i <= n; i++) {
         D = GenerateC (gen, r, s);
         /* Find the correlation */
         Ci = XS;
         while (Ci->C != D)
            Ci = Ci->Next;
         ++Ci->Nb;
      }

      /* Keep the observed numbers in sstring_Count */
      for (j = jlow; j <= jhigh; j++)
         res->Count[j] = 0;
      Ci = XS;
      j = 1;
      while (Ci) {
         if (j >= res->jmax)
            res->Count[res->jmax] += Ci->Nb;
         else
            res->Count[Loca[j]] += Ci->Nb;
         Ci = Ci->Next;
         ++j;
      }

      X2 = gofs_Chi2 (res->NbExp, res->Count, jlow, jhigh);
      statcoll_AddObs (res->sVal1, X2);
      if (swrite_Counters)
         tables_WriteTabL (res->Count, jlow, jhigh, 5, 10,
                           "Observed population counts");
   }

   res->degFree = V[0] = NbGroups - 1;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
                      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetChi2SumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);
   if (swrite_Basic) {
      swrite_AddStrChi (str, LEN1, NbGroups - 1);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTest (N, res);
      swrite_Final (gen, Timer);
   }
   DeleteCorr (corr);
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static double ProbabiliteLHR (long j, double Lnl)
/*
 * Returns the probability that the longest series of successive 1 has
 * length = j.
 */
{
   double x, temp;
   temp = (j + 1) * num_Ln2 - Lnl;
   x = exp (-exp (-temp));
   temp += num_Ln2;
   x = exp (-exp (-temp)) - x;
   return x;
}


/*-------------------------------------------------------------------------*/

static void WriteDataLongHead (unif01_Gen *gen, char *Test,
   long N, long n, int r, int s, long L)
{
   swrite_Head (gen, Test, N, n, r);
   printf (",   s = %1d,   L = %1ld\n\n", s, L);
}


/*-------------------------------------------------------------------------*/

void sstring_LongestHeadRun (unif01_Gen *gen, sstring_Res2 *res,
   long N, long n, int r, int s, long L)
{
   const double eps = DBL_EPSILON;
   const long K = L/s;             /* Number of iterations */
   long Rep;                       /* Current replication number */
   long Seq;                       /* Current sequence number */
   long i;
   double LnLen;                   /* log (L) or log (NnL) */
   double X2, temp;
   int j;
   long longest;                   /* Longest serie of 1 in a block */
   long longest2;                  /* Longest serie of 1 over all blocks */
   long longest3;                  /* Longest serie of 1 over all Replic. */
   long longueur;                  /* Run length in a block */
   long longueur2;                 /* Run length in a sequence */
   long longueur3;                 /* Run length in a replication */
   long jhigh;                     /* Highest class for Chi2 */
   long jhigh2;                    /* Highest index for CDF[j] */
   long NbGroups;                  /* Number of classes for Chi2 */
   bitset_BitSet ensemble;         /* Chosen bits in each generated number */
   double V[1];                    /* Number degrees of freedom for Chi2 */
   char str [LEN1 + 1];
   double NbExp [DIM + 1];         /* Expected numbers */
   double CDF [DIM + 1];           /* Cumulative probabilities */
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sstring_LongestHeadRun test";
   sres_Chi2 *Chi;
   sres_Disc *Disc;

   Timer = chrono_Create ();
   L = K * s;
   if (swrite_Basic)
      WriteDataLongHead (gen, TestName, N, n, r, s, L);
   util_Assert (L >= LMIN, "sstring_LongestHeadRun:   L < 1000");
   if (res == NULL) {
      localRes = TRUE;
      res = sstring_CreateRes2 ();
   }
   jhigh = DIM;

   /* Get the expected numbers for the chi-square for blocks of L bits */
   LnLen = log ((double) L);
   CDF[0] = ProbabiliteLHR (0, LnLen);
   NbExp[0] = n * CDF[0];
   for (j = 1; j < DIM; j++) {
      temp = ProbabiliteLHR (j, LnLen);
      NbExp[j] = n * temp;
      CDF[j] = temp + CDF[j-1];
      if  ((temp <= eps) && (CDF[j] > 0.5)) {
         jhigh = j;
         break;
      }
   }
   util_Assert (jhigh > 0, "sstring_LongestHeadRun:   jhigh = 0");
   NbExp[jhigh] = n * (1.0 - CDF[jhigh - 1]);

   /* Get the probabilities for the global run over the N*n*L bits */
   LnLen = log (N * (double) n * (double) L);   /* Avoid overflow of long */
   CDF[0] = ProbabiliteLHR (0, LnLen);
   for (j = 1; j < DIM; j++) {
      temp = ProbabiliteLHR (j, LnLen);
      CDF[j] = temp + CDF[j-1];
      if ((temp <= eps) && (CDF[j] > 0.5)) {
         jhigh2 = j;
         break;
      }
   }

   InitRes2 (res, N, jhigh);
   Disc = res->Disc;
   Chi = res->Chi;
   tables_CopyTabD (NbExp, Chi->NbExp, 0, jhigh);

   if (swrite_Classes)
      gofs_WriteClasses (Chi->NbExp, Chi->Loc, 0, jhigh, 0);
   gofs_MergeClasses (Chi->NbExp, Chi->Loc, &Chi->jmin, &Chi->jmax, &NbGroups);
   if (swrite_Classes)
      gofs_WriteClasses (Chi->NbExp, Chi->Loc, Chi->jmin, Chi->jmax, NbGroups);
   Chi->degFree = NbGroups - 1;
   if (Chi->degFree <= 0) {
      util_Assert (1, "sstring_LongestHeadRun:   Chi->degFree = 0");
      if (localRes)
         sstring_DeleteRes2 (res);
      return;
   }

   sprintf (str, "The N statistic values (a ChiSquare with %1ld degrees"
                 " of freedom):", Chi->degFree);
   statcoll_SetDesc (Chi->sVal1, str);
   statcoll_SetDesc (Disc->sVal1,
        "The longest run of 1 for each replication ");

   /* Beginning of test */
   longest3 = longueur3 = 0;
   for (Rep = 1; Rep <= N; Rep++) {
      for (i = Chi->jmin; i <= Chi->jmax; i++)
         Chi->Count[i] = 0;

      longest2 = -1;            /* -1 at the beginning of a new replication */
      longueur2 = 0;
      for (Seq = 1; Seq <= n; Seq++) {
	 longest = -1;            /* -1 at the beginning of a new sequence */
         longueur = 0;
         for (i = 1; i <= K; i++) {
            /* Now build a block of L bits */
            ensemble = unif01_StripB (gen, r, s);
            /* Examine each bit of a number */
            for (j = s - 1; j >= 0; j--) {
               if (bitset_TestBit (ensemble, j))
                  ++longueur;
               else {
		  /* Beginning of a sequence: merge last block of 1's of   */
		  /* last sequence with first block of 1's of new sequence */
		  if (longest < 0) {
		     /* Beginning of a replication: merge last sequence of */
		     /* 1's of last replication with first sequence of 1's */
		     /* of new replication */
                     if (longest2 < 0) {
                        longueur3 += longueur;
                        if (longueur3 > longest3)
                           longest3 = longueur3;
		     }
                     longueur2 += longueur;
                     if (longueur2 > longest2)
                        longest2 = longueur2;
		  }
                  if (longueur > longest)
                     longest = longueur;
                  longueur = 0;
               }
            }
         }
         if (longueur > longest)
            longest = longueur;
         if (longest >= Chi->jmax)
            ++Chi->Count[Chi->jmax];
         else if (longest <= Chi->jmin)
            ++Chi->Count[Chi->jmin];
         else
            ++Chi->Count[Chi->Loc[longest]];
         if (longest > longest2)
            longest2 = longest;
         longueur3 = longueur2 = longueur;
      }

      X2 = gofs_Chi2 (Chi->NbExp, Chi->Count, Chi->jmin, Chi->jmax);
      statcoll_AddObs (Chi->sVal1, X2);
      statcoll_AddObs (Disc->sVal1, (double) longest2);
      if (longest2 > longest3)
	 longest3 = longest2;
      if (swrite_Counters)
         tables_WriteTabL (Chi->Count, Chi->jmin, Chi->jmax, 5, 10,
                           "Observed numbers");
      longueur = 0;
      for (j = Chi->jmin; j <= Chi->jmax; j++)
         longueur += Chi->Count[j];
      util_Warning (longueur != n, "Total Count != n");
   }
   Disc->sVal2 = longest3;
   if (longest3 > jhigh2) {
      Disc->pLeft = 1.0;
      Disc->pRight = 0.0;
   } else {
      Disc->pLeft = CDF[longest3];
      if (longest3 > 0)
         Disc->pRight = 1.0 - CDF[longest3 - 1];
      else
         Disc->pRight = 1.0;
   }
   Disc->pVal2 = gofw_pDisc (Disc->pLeft, Disc->pRight);

   V[0] = Chi->degFree;
   gofw_ActiveTests2 (Chi->sVal1->V, Chi->pVal1->V, N, wdist_ChiSquare, V,
      Chi->sVal2, Chi->pVal2);
   Chi->pVal1->NObs = N;
   sres_GetChi2SumStat (Chi);

   if (swrite_Collectors) {
      statcoll_Write (Chi->sVal1, 5, 14, 4, 3);
      statcoll_Write (Disc->sVal1, 5, 14, 0, 0);
   }
   if (swrite_Basic) {
      swrite_AddStrChi (str, LEN1, Chi->degFree);
      gofw_WriteActiveTests2 (N, Chi->sVal2, Chi->pVal2, str);
      swrite_Chi2SumTest (N, Chi);
      printf ("-----------------------------------------------\n");
      printf ("Global longest run of 1               :");
      gofw_Writep2 (Disc->sVal2, Disc->pVal2);
      printf ("\n\n");
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sstring_DeleteRes2 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void HammingWeight2_L (unif01_Gen * gen, sres_Basic * res,
   long N, int r, int s, long L, long K)
/*
 * Generate all the n bits for the HammingWeight2 test in the case L > s. 
 * For the last number generated in a block of L bits, we keep its first
 * LMods bits and discard the other bits.
 */
{
   const int LDivs = L / s;         /* A block uses LDivs numbers ... */
   const int LMods = L % s;         /* + 1 if LMods > 0 */
   const double L2 = L / 2.0;
   int co, j;
   long i, Seq;
   unsigned long Z;
   double X2;

   for (Seq = 1; Seq <= N; Seq++) {
      X2 = 0.0;
      for (i = 0; i < K; i++) {
         /* Generate a block of L bits */
         co = 0;
         for (j = 0; j < LDivs; j++) {
            Z = unif01_StripB (gen, r, s);
            while (Z > 0) {               /* Count the number of 1 bits */
               Z &= Z - 1;                /* Clear lowest 1 bit */
               ++co;
            }
         }
         /* The last bits of the block */
         if (LMods > 0) {
            Z = unif01_StripB (gen, r, LMods);
            while (Z > 0) {
               Z &= Z - 1;
               ++co;
            }
         }
         X2 += (co - L2)*(co - L2);
      }
      X2 *= 4.0 / L;
      statcoll_AddObs (res->sVal1, X2);
   }
}


/*-------------------------------------------------------------------------*/

static void HammingWeight2_S (unif01_Gen * gen, sres_Basic * res,
   long N, int r, int s, long L, long K)
/*
 * Generate all the n bits for the HammingWeight2 test in the case L <= s. 
 * A number generates sDivL blocks. If s % L == 0, we use all s bits of the
 * number.
 */
{
   const int sDivL = s / L;         /* A number generates sDivL blocks */
   const long Q = K / sDivL + (K % sDivL > 0);
   const unsigned long MASK = num_TwoExp[L] - 1.0;
   const double L2 = L / 2.0;
   int co, j;
   long i, Seq;
   unsigned long Z, Y;
   double X2;

   for (Seq = 1; Seq <= N; Seq++) {
      X2 = 0.0;
      for (i = 0; i < Q; i++) {
         Z = unif01_StripB (gen, r, s);

         /* Generate sDivL blocks of L bits */
         for (j = 0; j < sDivL; j++) {
            co = 0;
            Y = Z & MASK;
            while (Y > 0) {        /* Count the number of 1 bits */
               Y &= Y - 1;         /* Clear lowest 1 bit */
               ++co;
            }
            X2 += (co - L2)*(co - L2);
            Z >>= L;
         }
      }
      X2 *= 4.0 / L;
      statcoll_AddObs (res->sVal1, X2);
   }
}


/*-------------------------------------------------------------------------*/

void sstring_HammingWeight2 (unif01_Gen * gen, sres_Basic * res,
   long N, long n, int r, int s, long L)
{
   const long K = n / L;
   double sum;
   double V[1];                   /* Number of Chi2 degrees of freedom */
   char chaine[LEN1 + 1] = "";
   char str[LEN2 + 1] = "";
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sstring_HammingWeight2 test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataLongHead (gen, TestName, N, n, r, s, L);
   util_Assert (r + s <= 32, "sstring_HammingWeight2:   r + s > 32");
   util_Assert (L <= n, "sstring_HammingWeight2:   L > n");
   util_Assert (L >= 2, "sstring_HammingWeight2:   L < 2");

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateBasic ();
   }
   sres_InitBasic (res, N, "sstring_HammingWeight2");
   strncpy (chaine, "sVal1:   a chi-square with ", (size_t) LEN1);
   sprintf (str, "%ld", K);
   strncat (chaine, str, (size_t) LEN2);
   strncat (chaine, " degrees of freedom", (size_t) LEN1);
   statcoll_SetDesc (res->sVal1, chaine);

   if (L >= s)
      HammingWeight2_L (gen, res, N, r, s, L, K);
   else
      HammingWeight2_S (gen, res, N, r, s, L, K);

   V[0] = K;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
                      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sum = N * statcoll_Average (res->sVal1);
   res->sVal2[gofw_Sum] = sum;
   res->pVal2[gofw_Sum] = fbar_ChiSquare2 (N*K, 12, sum);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 2, 1);
   if (swrite_Basic) {
      swrite_AddStrChi (str, LEN2, K);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTestb (N, res->sVal2[gofw_Sum], res->pVal2[gofw_Sum], K);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeleteBasic (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void HammingWeight_L (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r, int s, long L)
/*
 * Generate all the n*L bits for the HammingWeight test in the case L > s. 
 * For the last number generated in a block of L bits, we keep its first
 * LMods bits and discard the other bits.
 */
{
   const int LDivs = L / s;         /* A block uses LDivs numbers ... */
   const int LMods = L % s;         /* + 1 if LMods > 0 */
   int co, j;
   long i, Seq;
   unsigned long Z;
   double X2;

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = res->jmin; i <= res->jmax; i++)
         res->Count[i] = 0;

      for (i = 0; i < n; i++) {
         /* Generate a block of L bits */
         co = 0;
         for (j = 0; j < LDivs; j++) {
            Z = unif01_StripB (gen, r, s);
            while (Z > 0) {               /* Count the number of 1 bits */
               Z &= Z - 1;                /* Clear lowest 1 bit */
               ++co;
            }
         }
         /* The last bits of the block */
         if (LMods > 0) {
            Z = unif01_StripB (gen, r, LMods);
            while (Z > 0) {
               Z &= Z - 1;
               ++co;
            }
         }
         ++res->Count[res->Loc[co]];
      }

      X2 = gofs_Chi2 (res->NbExp, res->Count, res->jmin, res->jmax);
      statcoll_AddObs (res->sVal1, X2);
      if (swrite_Counters)
         tables_WriteTabL (res->Count, res->jmin, res->jmax, 5, 10,
                           "Observed numbers of blocks");
   }
}


/*-------------------------------------------------------------------------*/

static void HammingWeight_S (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r, int s, long L)
/*
 * Generate all the n*L bits for the HammingWeight test in the case L <= s. 
 * A number generates sDivL blocks. If s % L == 0, we use all s bits of the
 * number.
 */
{
   const int sDivL = s / L;         /* A number generates sDivL blocks */
   const int s1 = s - s % L;
   const long Q = n / sDivL;
   const int Q2 = n % sDivL;
   const unsigned long MASK = num_TwoExp[L] - 1.0;
   int co, j;
   long i, Seq;
   unsigned long Z, Y;
   double X2;

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = res->jmin; i <= res->jmax; i++)
         res->Count[i] = 0;

      for (i = 0; i < Q; i++) {
         Z = unif01_StripB (gen, r, s1);

         /* Generate sDivL blocks of L bits */
         for (j = 0; j < sDivL; j++) {
            co = 0;
            Y = Z & MASK;
            while (Y > 0) {        /* Count the number of 1 bits */
               Y &= Y - 1;         /* Clear lowest 1 bit */
               ++co;
            }
            ++res->Count[res->Loc[co]];
            Z >>= L;
         }
      }

      /* The last bits */
      if (Q2 > 0) {
	 Z = unif01_StripB (gen, r, Q2 * L);
         for (j = 0; j < Q2; j++) {
            co = 0;
            Y = Z & MASK;
            while (Y > 0) {        /* Count the number of 1 bits */
               Y &= Y - 1;         /* Clear lowest 1 bit */
               ++co;
            }
            ++res->Count[res->Loc[co]];
            Z >>= L;
         }
      }

      X2 = gofs_Chi2 (res->NbExp, res->Count, res->jmin, res->jmax);
      statcoll_AddObs (res->sVal1, X2);
      if (swrite_Counters)
         tables_WriteTabL (res->Count, res->jmin, res->jmax, 5, 10,
                           "Observed numbers of blocks");
   }
}


/*-------------------------------------------------------------------------*/

void sstring_HammingWeight (unif01_Gen * gen, sres_Chi2 * res,
   long N, long n, int r, int s, long L)
{
   long i;
   double V[1];                   /* Number of Chi2 degrees of freedom */
   char str[LEN1 + 1] = "";
   fmass_INFO Q;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   long jlow, jhigh;
   long NbGroups;                 /* Number of classes */
   char *TestName = "sstring_HammingWeight test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataLongHead (gen, TestName, N, n, r, s, L);
   util_Assert (r + s <= 32, "sstring_HammingWeight:   r + s > 32");
   util_Assert (L >= 2, "sstring_HammingWeight:   L < 2");

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateChi2 ();
   }
   sres_InitChi2 (res, N, L, "sstring_HammingWeight");

   Q = fmass_CreateBinomial (L, 0.5, 0.5);
   for (i = 0; i <= L; i++)
      res->NbExp[i] = n * fmass_BinomialTerm2 (Q, i);
   fmass_DeleteBinomial (Q);

   jlow = 0;
   jhigh = L;
   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, res->Loc, jlow, jhigh, 0);
   gofs_MergeClasses (res->NbExp, res->Loc, &jlow, &jhigh, &NbGroups);
   if (swrite_Classes)
      gofs_WriteClasses (res->NbExp, res->Loc, jlow, jhigh, NbGroups);
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

   if (L >= s)
      HammingWeight_L (gen, res, N, n, r, s, L);
   else
      HammingWeight_S (gen, res, N, n, r, s, L);

   V[0] =  res->degFree;
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_ChiSquare, V,
                      res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetChi2SumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 2, 1);
   if (swrite_Basic) {
      swrite_AddStrChi (str, LEN1, res->degFree);
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2, str);
      swrite_Chi2SumTest (N, res);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeleteChi2 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

#if 0
void sstring_Run0 (unif01_Gen * gen, sres_Basic * res,
   long N, long n, int r, int s)
{
   const long K = n / s;          /* If n % s != 0, a string will contain 
                                     K * s bits instead of n */
   const unsigned long SBIT = 1UL << (s - 1);
   unsigned long jBit;            /* Position of current bit in Z */
   int pBit;                      /* Previous bit */
   long i, Seq;
   long co1;                      /* Counter for number of 1 */
   long cor;                      /* Counter for number of runs */
   unsigned long Z;
   double X, f1;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sstring_Run test";

   Timer = chrono_Create ();
   n = K * s;
   if (swrite_Basic)
      WriteDataPeriod (gen, TestName, N, n, r, s);

   util_Assert (r + s <= 32, "sstring_Run:   r + s > 32");
   /*   util_Assert (100 <= n, "sstring_Run:   n < 100"); */

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateBasic ();
   }
   sres_InitBasic (res, N, "sstring_Run");
   statcoll_SetDesc (res->sVal1, "sVal1:   a standard normal");

   for (Seq = 1; Seq <= N; Seq++) {
      co1 = cor = 0;
      /* Be sure to count the first run with pBit != {0, 1} */
      pBit = 2;
      for (i = 0; i < K; i++) {
         Z = unif01_StripB (gen, r, s);
         jBit = SBIT;

         /* Add the number of 1 bit and number of runs in Z */
         while (jBit > 0) {
	    if (Z & jBit) {                /* bit 1 */
               co1++;
               if (pBit != 1)
                  cor++;
               pBit = 1;
            } else {                       /* bit 0 */
               if (pBit != 0)
                  cor++;
               pBit = 0;
            }
            jBit >>= 1;
         }
      }
      f1 = (double) co1 / (K * s);
      X = (cor - n * 2.0 * f1 * (1.0 - f1)) /
          (2.0 * sqrt ((double) n) * f1 * (1.0 - f1));
      statcoll_AddObs (res->sVal1, X);
   }
   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_Normal,
       (double *) NULL, res->sVal2, res->pVal2);
   res->pVal1->NObs = N;

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2,
         "Normal statistic                      :");
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sres_DeleteBasic(res);
   chrono_Delete (Timer);
}
#endif


/*=========================================================================*/

void sstring_Run (unif01_Gen * gen, sstring_Res3 *res,
   long N, long n, int r, int s)
{
   const unsigned long SBIT = 1UL << (s - 1);
   const double sr = s;
   unsigned long jBit;            /* Position of current bit in Z */
   int pBit;                      /* Previous bit */
   int k, j;
   long Seq;
   double cob;                    /* Counter for number of bits */
   long cor;                      /* Counter for number of 1 runs */
   int len;                       /* Length of current run */
   unsigned long Z;
   double X2, X, temp;
   char str[LEN1 + 1];
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sstring_Run test";
   sres_Basic *NBits;
   sres_Chi2 *NRuns;
   long *Count0,  *Count1;
   double *Prob, *NbExp;
   double Param[1];

   Timer = chrono_Create ();
   k = 1 + num_Log2 (n / gofs_MinExpected);
   if (swrite_Basic)
      WriteDataPeriod (gen, TestName, N, n, r, s);
   util_Assert (r + s <= 32, "sstring_Run:   r + s > 32");
   /*   util_Assert (100 <= n, "sstring_Run:   n < 100"); */

   if (res == NULL) {
      localRes = TRUE;
      res = sstring_CreateRes3 ();
   }
   InitRes3 (res, N, k);
   NBits = res->NBits;
   NRuns = res->NRuns;
   Count0 = res->Count0;
   Count1 = res->Count1;

   statcoll_SetDesc (NBits->sVal1,
       "The N statistic values (a standard normal):");
   sprintf (str, "The N statistic values (a ChiSquare with %1d degrees"
                 " of freedom):", 2*(k - 1));
   statcoll_SetDesc (NRuns->sVal1, str);

   Prob = util_Calloc (1 + (size_t) k, sizeof (double));
   Prob[0] = 1.0;
   for (j = 1; j < k; j++) {
      Prob[j] = Prob[j - 1] / 2.0; 
      NRuns->NbExp[j] = n * Prob[j];
   }
   Prob[k] = Prob[k - 1]; 
   NRuns->NbExp[k] = n * Prob[k];
   util_Assert (NRuns->NbExp[k] >= gofs_MinExpected,
        "sstring_Run:   NRuns->NbExp[k] < gofs_MinExpected");

   if (swrite_Classes)
      gofs_WriteClasses (NRuns->NbExp, NRuns->Loc, 1, k, 0);
   NRuns->jmax = k;
   NRuns->jmin = 1;
   NRuns->degFree = 2*(k - 1);
   if (NRuns->degFree < 1) {
      util_Warning (TRUE, "Chi-square with 0 degree of freedom.");
      if (localRes)
         sstring_DeleteRes3 (res);
      return;
   }

   for (Seq = 1; Seq <= N; Seq++) {
      cob = cor = len = 0;
      for (j = 1; j <= k; j++) {
	 Count0[j] = 0;
	 Count1[j] = 0;
      }

      /* Make sure to count the first run; set pBit != {0, 1} */
      pBit = 2;
      while (cor < n) {
         Z = unif01_StripB (gen, r, s);
         jBit = SBIT;
         cob += sr;
         if (len >= n) {
	    util_Warning (TRUE, "sstring_Run:   all bits are 0 !");
            util_Free (Prob);
            if (localRes)
               sstring_DeleteRes3 (res);
            return;
	 }

         /* Add the number of runs in Z */
         while (jBit > 0) {
	    if (Z & jBit) {                /* bit 1 */
	       if (pBit != 1) {
                  cor++;
                  if (len < k)
		     (Count0[len])++;
                  else
		     (Count0[k])++;
                  len = 1;
	       } else {
		 len++;
	       }
               pBit = 1;
            } else {                       /* bit 0 */
	       if (pBit != 0) {
                  if (len < k)
		     (Count1[len])++;
                  else
		     (Count1[k])++;
                  len = 1;
	       } else {
		 len++;
	       }
               pBit = 0;
            }
            jBit >>= 1;
         }
      }

      X2 = 0.0;
      NbExp = NRuns->NbExp;
      for (j = NRuns->jmin; j <= NRuns->jmax; j++) {
	 temp = Count0[j] - NbExp[j];
         X2 += temp * temp / (NbExp[j] * (1.0 - Prob[j]));
      }
      X = X2;
      X2 = 0.0;
      for (j = NRuns->jmin; j <= NRuns->jmax; j++) {
	 temp = Count1[j] - NbExp[j];
         X2 += temp * temp / (NbExp[j] * (1.0 - Prob[j]));
      }
      statcoll_AddObs (NRuns->sVal1, X2 + X);

      if (swrite_Counters) {
         tables_WriteTabL (Count0, 1, k, 5, 10,
             "Observed number of runs of 0");
         tables_WriteTabL (Count1, 1, k, 5, 10,
             "Observed number of runs of 1");
      }

      X = (cob - 4.0 * n) / sqrt (8.0 * n);
      statcoll_AddObs (NBits->sVal1, X);
   }

   Param[0] = 2*(k - 1);
   gofw_ActiveTests2 (NRuns->sVal1->V, NRuns->pVal1->V, N, wdist_ChiSquare,
      Param, NRuns->sVal2, NRuns->pVal2);
   NRuns->pVal1->NObs = N;
   sres_GetChi2SumStat (NRuns);

   gofw_ActiveTests2 (NBits->sVal1->V, NBits->pVal1->V, N, wdist_Normal,
       (double *) NULL, NBits->sVal2, NBits->pVal2);
   NBits->pVal1->NObs = N;
   sres_GetNormalSumStat (NBits);


   if (swrite_Basic) {
      printf ("\n-----------------------------------------------\n");
      if (N == 1) {
         printf ("Total number of 1 runs:  %1ld\n\n", cor);
         printf ("Number of degrees of freedom          : %4ld\n",
                 NRuns->degFree);
         printf ("Chi2 statistic for number of runs     :");
         gofw_Writep2 (NRuns->sVal2[gofw_Mean], NRuns->pVal2[gofw_Mean]);
      } else {
         printf ("Test results for the number of runs:\n");
         gofw_WriteActiveTests0 (N, NRuns->sVal2, NRuns->pVal2);
         swrite_Chi2SumTest (N, NRuns);
      }
      if (swrite_Collectors)
         statcoll_Write (NRuns->sVal1, 5, 14, 4, 3);

      printf ("\n-----------------------------------------------\n");
      if (N == 1) {
         printf ("Total number of bits:  %.0f\n\n", cob);
         printf ("Normal statistic for number of bits   :");
         gofw_Writep2 (NBits->sVal2[gofw_Mean], NBits->pVal2[gofw_Mean]);
      } else {
         printf ("Test results for the number of bits:\n");
         gofw_WriteActiveTests0 (N, NBits->sVal2, NBits->pVal2);
         swrite_NormalSumTest (N, NBits);
      }
      if (swrite_Collectors)
         statcoll_Write (NBits->sVal1, 5, 14, 4, 3);

      printf ("\n\n");
      swrite_Final (gen, Timer);
   }
   util_Free (Prob);
   if (localRes)
      sstring_DeleteRes3 (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataAutoCor (unif01_Gen *gen, char *Test,
   long N, long n, int r, int s, int d)
{
   swrite_Head (gen, Test, N, n, r);
   printf (",   s = %1d,   d = %1d\n\n", s, d);
}


/*-------------------------------------------------------------------------*/

void sstring_AutoCor (unif01_Gen * gen, sres_Basic * res,
   long N, long n, int r, int s, int d)
{
   const long K = (n - d) / s;
   const long M = d / s + 2;
   unsigned long *Y;              /* Circular buffer for random numbers */
   unsigned long A;               /* Correlation */
   unsigned long Z, s1, s2;
   unsigned long mask1, mask2;    /* Masks of s1, s2 least sig. bits */
   double X;
   long i, Seq;
   int j1, j2;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sstring_AutoCor test";

   Timer = chrono_Create ();
   /* There are a few bits less than n */
   n -= (n - d) % s;
   if (swrite_Basic)
      WriteDataAutoCor (gen, TestName, N, n, r, s, d);

   util_Assert (r + s <= 32, "sstring_AutoCor:   r + s > 32");
   util_Assert (d <= n / 2, "sstring_AutoCor:   d > n/2");
   util_Assert (d > 0, "sstring_AutoCor:   d < 1");

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateBasic ();
   }
   sres_InitBasic (res, N, "sstring_AutoCor");
   Y = util_Calloc ((size_t) M, sizeof (unsigned long));
   statcoll_SetDesc (res->sVal1, "sVal1:   a standard normal");
   s1 = d % s;
   s2 = s - s1;
   mask1 = num_TwoExp[s1] - 1.0;
   mask2 = num_TwoExp[s2] - 1.0;

   for (Seq = 1; Seq <= N; Seq++) {
      /* Fill circular buffer with first random numbers */
      for (i = 0; i < M-1; i++)
         Y[i] = unif01_StripB (gen, r, s);

      A = 0;
      j1 = M - 1;
      j2 = M - 2;
      for (i = 0; i < K; i++) {
         Y[j1] = unif01_StripB (gen, r, s);
         j1 = (j1 + 1) % M;
         Z = ((Y[j1] >> s1) ^ Y[j2]) & mask2;
         while (Z > 0) {          /* Count the number of 1 bits in Z */
            Z &= Z - 1;           /* Clear lowest 1 bit */
            ++A;
         }
         j2 = (j2 + 1) % M;

         Z = ((Y[j2] >> s2) ^ Y[j1]) & mask1;
         while (Z > 0) {
            Z &= Z - 1;
            ++A;
         }
      }

      X = 2.0 * (A - (n - d) / 2.0) / sqrt ((double) (n - d));
      statcoll_AddObs (res->sVal1, X);
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
   util_Free (Y);
   if (localRes)
      sres_DeleteBasic (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataHammingCorr (unif01_Gen *gen, char *TestName,
   long N, long n, int r, int s, int L)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   s = %1d,   L = %1d\n\n\n", s, L);
}


/*-------------------------------------------------------------------------*/

static void HammingCorr_L (unif01_Gen *gen, sstring_Res * res,
   long n, int r, int s, int L)
/*
 * Generate all the n bits for the HammingCorr test in the case L > s. 
 * For the last number generated in a block of L bits, we keep its first
 * LMods bits and discard the other bits.
 */
{
   const int LMods = L % s;
   const int LDivs = L / s;
   int Pre, X;
   int j;
   long k;
   unsigned long Z;

   /* Junk value to avoid a test "if (k == 1)" for every generated
      number; it will not be counted. */
   Pre = L + 1;
   for (k = 1; k <= n; k++) {
      /* Generate a sequence of L bits */
      X = 0;
      for (j = 1; j <= LDivs; j++) {
	 Z = unif01_StripB (gen, r, s);
	 /* Count the number of 1 bits */
	 while (Z > 0) {
	    Z &= Z - 1;                /* Clear lowest 1 bit */
	    ++X;
	 }
      }
      /* The last bits of the sequence */
      if (LMods > 0) {
	 Z = unif01_StripB (gen, r, LMods);
	 while (Z > 0) {
	    Z &= Z - 1;
	    ++X;
	 }
      }
      ++res->Counters[Pre][X];
      Pre = X;
   }
}


/*-------------------------------------------------------------------------*/

static void HammingCorr_S (unif01_Gen *gen, sstring_Res * res,
   long n, int r, int s, int L)
/*
 * Generate all the n bits for the HammingCorr test in the case L <= s. 
 * A number generates sDivL blocks. If s % L == 0, we use all s bits of the
 * number.
 */
{
   const int sDivL = s / L;         /* A number generates sDivL blocks */
   const long Q = n / sDivL;
   const long Q1 = n % sDivL;
   const unsigned long MASK = num_TwoExp[L] - 1.0;
   int Pre, X;
   int j;
   long k;
   unsigned long Z, Y;

   /* Junk value to avoid a test "if (k == 1)" for every generated
      number; it will not be counted. */
   Pre = L + 1;
   for (k = 0; k < Q; k++) {
      Z = unif01_StripB (gen, r, s);

      for (j = 0; j < sDivL; j++) {
	 X = 0;
	 Y = Z & MASK;
	 while (Y > 0) {        /* Count the number of 1 bits */
	    Y &= Y - 1;         /* Clear lowest 1 bit */
	    ++X;
	 }
         ++res->Counters[Pre][X];
         Pre = X;
	 Z >>= L;
      }
   }

   /* The last Q1 blocks */
   if (Q1 > 0) {
      Z = unif01_StripB (gen, r, s);
      for (j = 0; j < Q1; j++) {
	 X = 0;
	 Y = Z & MASK;
	 while (Y > 0) {        /* Count the number of 1 bits */
	    Y &= Y - 1;         /* Clear lowest 1 bit */
	    ++X;
	 }
	 ++res->Counters[Pre][X];
	 Pre = X;
	 Z >>= L;
      }
   }
}


/*-------------------------------------------------------------------------*/

void sstring_HammingCorr (unif01_Gen * gen, sstring_Res * res,
   long N, long n, int r, int s, int L)
{
   int i, j;
   long Seq;
   double Sum;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sstring_HammingCorr test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataHammingCorr (gen, TestName, N, n, r, s, L);

   util_Assert (s <= num_MaxTwoExp, "sstring_HammingCorr:   s too large");
   util_Assert ((unsigned) s <= CHAR_BIT * sizeof (unsigned long),
                "sstring_HammingCorr:   s too large");

   if (res == NULL) {
      localRes = TRUE;
      res = sstring_CreateRes ();
   }
   InitRes (res, N, L, -1, "sstring_HammingCorr");
   statcoll_SetDesc (res->Bas->sVal1, "HammingCorr sVal1:   standard normal");

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 0; i <= L; i++) {
         for (j = 0; j <= L; j++)
            res->Counters[i][j] = 0;
      }
      if (L >= s)
	 HammingCorr_L (gen, res, n, r, s, L);
      else
	 HammingCorr_S (gen, res, n, r, s, L);

      if (swrite_Counters)
         /* Print the matrix of counters */
         tables_WriteMatrixL (res->Counters, 0, L, 0, L, 8,
                              res->Style, "Number of pairs [0..L, 0..L]");

      /* Calculate statistic */
      Sum = 0.0;
      for (i = 0; i <= L; i++) {
         for (j = 0; j <= L; j++)
            Sum += res->Counters[i][j] * (i - L / 2.0) * (j - L / 2.0);
      }
      Sum = Sum * 4.0 / (L * sqrt (n - 1.0));
      statcoll_AddObs (res->Bas->sVal1, Sum);
   }

   gofw_ActiveTests2 (res->Bas->sVal1->V, res->Bas->pVal1->V, N,
      wdist_Normal, (double *) NULL, res->Bas->sVal2, res->Bas->pVal2);
   res->Bas->pVal1->NObs = N;
   sres_GetNormalSumStat (res->Bas);

   if (swrite_Collectors)
      statcoll_Write (res->Bas->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      gofw_WriteActiveTests2 (N, res->Bas->sVal2, res->Bas->pVal2,
         "Normal statistic                      :");
      swrite_NormalSumTest (N, res->Bas);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sstring_DeleteRes (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteDataHammingIndep (unif01_Gen * gen, char *TestName,
   long N, long n, int r, int s, int L, int d)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   s = %1d,   L = %1d,   d = %1d\n\n\n", s, L, d);
}


/*-------------------------------------------------------------------------*/

static void HammingIndep_L (unif01_Gen *gen, sstring_Res * res,
   long n, int r, int s, int L)
/*
 * Generate all the n bits for the HammingIndep test in the case L > s. 
 * For the last number generated in a block of L bits, we keep its first
 * LMods bits and discard the other bits.
 */
{
   int Pre;                         /* Previous value of X */
   int X;
   const int LMods = L % s;
   const int LDivs = L / s;
   int j;
   unsigned long U;
   unsigned long TwonUL;
   unsigned long ic;

   /* For the test with n >= 2^30 */
   TwonUL = 2 * (unsigned long) n;

   Pre = 0;                        /* Eliminate a warning from compiler */
   for (ic = 1; ic <= TwonUL; ic++) {
      /* Generate 1 block of L bits */
      X = 0;
      for (j = 1; j <= LDivs; j++) {
	 U = unif01_StripB (gen, r, s);
	 /* Count the number of 1 bits */
	 while (U > 0) {
	    U &= U - 1;        /* Clear lowest 1 bit */
	    ++X;
	 }
      }
      /* The last bits of the block */
      if (LMods > 0) {
	 U = unif01_StripB (gen, r, LMods);
	 while (U > 0) {
	    U &= U - 1;
	    ++X;
	 }
      }
      /* Non-overlapping pairs; count only when ic % 2 == 0 */
      if (!(ic & 1))
	 ++res->Counters[Pre][X];
      Pre = X;
   }
}


/*-------------------------------------------------------------------------*/

static void HammingIndep_S (unif01_Gen *gen, sstring_Res * res,
   long n, int r, int s, int L)
/*
 * Generate all the n bits for the HammingIndep test in the case L <= s. 
 * A number generates sDivL blocks. If s % L == 0, we use all s bits of the
 * number.
 */
{
   const int sDivL = s / L;         /* A number generates sDivL blocks */
   const unsigned long MASK = num_TwoExp[L] - 1.0;
   int Pre;                         /* Previous value of X */
   int X;
   int j;
   unsigned long Q, Q1, i;
   unsigned long Z, Y;
   unsigned long TwonUL;
   unsigned long bloc = 0;

   /* For the test with n >= 2^30 */
   TwonUL = 2 * (unsigned long) n;
   Q = TwonUL / sDivL;
   Q1 = TwonUL % sDivL;

   Pre = 0;                     /* Eliminate a warning from compiler */
   for (i = 0; i < Q; i++) {
      Z = unif01_StripB (gen, r, s);

      for (j = 0; j < sDivL; j++) {
	 X = 0;
	 Y = Z & MASK;
	 while (Y > 0) {        /* Count the number of 1 bits */
	    Y &= Y - 1;         /* Clear lowest 1 bit */
	    ++X;
	 }
         /* Non-overlapping pairs; count only when bloc % 2 == 0 */
	 if (!(++bloc & 1))
	    ++res->Counters[Pre][X];
         Pre = X;
	 Z >>= L;
      }
   }

   /* The last Q1 blocks */
   if (Q1 > 0)
      Z = unif01_StripB (gen, r, s);
   for (i = 0; i < Q1; i++) {
      X = 0;
      Y = Z & MASK;
      while (Y > 0) {        /* Count the number of 1 bits */
	 Y &= Y - 1;         /* Clear lowest 1 bit */
	 ++X;
      }
      if (!(++bloc & 1))
         ++res->Counters[Pre][X];
      Pre = X;
      Z >>= L;
   }
}


/*-------------------------------------------------------------------------*/

static void CountBlocks (
   sstring_Res *res,
   int L,                     /* Length of blocks (num bits) */
   int d                      /* Rank of rows-columns counted */

)
/*
 * Add the diagonal blocks for the matrix res->Counters, i.e. count the
 * number of values in the 4 corners of the matrix with the center removed.
 * We count only rows-columns of rank >= d starting from the center of
 * the matrix. We eliminate 2d - 1 rows and columns at the center of the
 * matrix when L is even. When L is odd, we eliminate 2d - 2 rows and
 * columns at the center (except when d = 1, where we keep all the rows
 * and columns).
 */
{
   int L2, L1, k, j, i;

   L2 = L1 = L / 2;
   if ((L & 1))
      ++L1;

   for (k = 1; k <= d; k++) {
      res->XD[k][0] = 0;
      res->XD[k][1] = 0;

      /* Block ++ */
      for (i = 0; i <= L1 - k; i++) {
         for (j = 0; j <= L1 - k; j++) {
            res->XD[k][0] += res->Counters[i][j];
         }
      }
      /* Block -- */
      for (i = L2 + k; i <= L; i++) {
         for (j = L2 + k; j <= L; j++) {
            res->XD[k][0] += res->Counters[i][j];
         }
      }
      /* Block +- */
      for (i = 0; i <= L1 - k; i++) {
         for (j = L2 + k; j <= L; j++) {
            res->XD[k][1] += res->Counters[i][j];
         }
      }
      /* Block -+ */
      for (i = L2 + k; i <= L; i++) {
         for (j = 0; j <= L1 - k; j++) {
            res->XD[k][1] += res->Counters[i][j];
         }
      }
   }
}

/*-------------------------------------------------------------------------*/

static void WriteBlocs (
   sstring_Res *res,
   int d                 /* Rank of rows-columns */
   )
/* 
 * Print the sum of counters in the diagonal blocks for different d
 */
{
   int i;
   printf ("--------------------------------------------------\n");

   for (i = 1; i <= d; i++) {
      printf ("The number of blocks ++, -- with d >= %1d is %10ld\n",
              i, res->XD[i][0]);
      printf ("The number of blocks +-, -+ with d >= %1d is %10ld\n\n",
              i, res->XD[i][1]);
   }
   printf ("\n");
}

/*-------------------------------------------------------------------------*/

void sstring_HammingIndep (unif01_Gen * gen, sstring_Res * res,
   long N, long n, int r, int s, int L, int d)
{
   int Liber = 0;           /* Num of degrees of freedom for main test */
   int Liberte;             /* Num of degrees of freedom for block tests */
   int X;
   int i, j;
   long Seq;
   double NbEsp;
   double Var;
   double NbMoyen;
   double X2;
   double Sum;
   const double nLR = n;
   double Z;
   double *Prob;
   double *NbEsp5;
   long *Count5;
   double V[1];
   fmass_INFO Q;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sstring_HammingIndep test";
   char chaine[LEN1 + 1] = "";
   char str[LEN2 + 1];

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataHammingIndep (gen, TestName, N, n, r, s, L, d);

   if (n < 2.0 * gofs_MinExpected) {
      util_Warning (TRUE, "sstring_HammingIndep:   n < 20");
      return;
   }
   /*  util_Assert (n >= 30, "sstring_HammingIndep:   n < 30"); */
   util_Assert (d <= sstring_MAXD, "sstring_HammingIndep:   d > sstring_MAXD");
   util_Assert (((L + 1) / 2) >= d, "sstring_HammingIndep:   d > (L + 1) / 2");
   util_Assert (s <= num_MaxTwoExp, "sstring_HammingIndep:   s too large");
   util_Assert ((unsigned) s <= CHAR_BIT * sizeof (unsigned long),
      "sstring_HammingIndep:   s too large");

   if (res == NULL) {
      localRes = TRUE;
      res = sstring_CreateRes ();
   }
   InitRes (res, N, L, d, "sstring_HammingIndep");

   for (i = 1; i <= d; i++) {
      strncpy (chaine, "HammingIndep Block[", (size_t) LEN1);
      sprintf (str, "%1d", i);
      strncat (chaine, str, (size_t) LEN2);
      strncat (chaine, "]", (size_t) 2);
      statcoll_SetDesc (res->Block[i]->sVal1, chaine);
   }
   strncpy (chaine, "\nCounters with expected numbers >= ", (size_t) LEN1);
   sprintf (str, "%g", gofs_MinExpected);
   strncat (chaine, str, (size_t) LEN2);
   statcoll_SetDesc (res->Bas->sVal1, chaine);

   Prob   = util_Calloc ((size_t) L + 1, sizeof (double));
   NbEsp5 = util_Calloc ((size_t) (L + 1)*(L + 1) + 1, sizeof (double));
   Count5 = util_Calloc ((size_t) (L + 1)*(L + 1) + 1, sizeof (long));

   Q = fmass_CreateBinomial (L, 0.5, 0.5);
   for (i = 0; i <= L; i++)
      Prob[i] = fmass_BinomialTerm2 (Q, i);
   fmass_DeleteBinomial (Q);

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 0; i <= L; i++) {
         for (j = 0; j <= L; j++) {
            res->Counters[i][j] = 0;
         }
      }
      if (L >= s)
	 HammingIndep_L (gen, res, n, r, s, L);
      else
	 HammingIndep_S (gen, res, n, r, s, L);

      /* The cells for which the expected number >= gofs_MinExpected will */
      /* correspond to one class each. Merge all other cells such that    */
      /* their expected number is < gofs_MinExpected into 1 big class: Z  */
      /* will contain the sum of all their expected numbers, and the      */
      /* counter X will contain the sum of all their observed numbers. We */
      /* shall apply a chi-square test with Liber degrees of freedom. We  */
      /* have (Liber + 1) classes. */
      Liber = 0;
      Z = 0.0;
      X = 0;
      for (i = 0; i <= L; i++) {
         for (j = 0; j <= L; j++) {
            res->ZCounters[i][j] = nLR * Prob[i] * Prob[j];
            if (res->ZCounters[i][j] >= gofs_MinExpected) {
               NbEsp5[Liber] = res->ZCounters[i][j];
               Count5[Liber] = res->Counters[i][j];
               ++Liber;
            } else {
               Z += res->ZCounters[i][j];
               X += res->Counters[i][j];
            }
         }
      }
      if (Z >= gofs_MinExpected) {
	/* We have one more class */
	/*  if (swrite_Classes && Seq == 1) {
            printf ("All cells with NbExp < ");
            printf ("%s", str);
            printf (" are merged into one big class with NbExp = %f\n\n", Z);
         } */
         NbEsp5[Liber] = Z;
         Count5[Liber] = X;
      } else if (Liber > 0) {
         /* We add them to the last class instead */
         --Liber;
         NbEsp5[Liber] += Z;
         Count5[Liber] += X;
      }

      /* Everything has been put in a single class; separate all in two
         classes */
      if (Liber == 0) {
	 Z = 0.0;
	 X = 0;
	 for (i = 0; i <= L; i++) {
	    for (j = 0; j <= L / 2; j++) {
	       res->ZCounters[i][j] = nLR * Prob[i] * Prob[j];
	       Z += res->ZCounters[i][j];
	       X += res->Counters[i][j];
	    }
	 }
         NbEsp5[0] = Z;
         Count5[0] = X;
	 Z = 0.0;
	 X = 0;
	 for (i = 0; i <= L; i++) {
	    for (j = 1 + L / 2; j <= L; j++) {
	       res->ZCounters[i][j] = nLR * Prob[i] * Prob[j];
	       Z += res->ZCounters[i][j];
	       X += res->Counters[i][j];
	    }
	 }
         NbEsp5[1] = Z;
         Count5[1] = X;
	 Liber = 1;
      }
      if (Liber > 0) {
         X2 = gofs_Chi2 (NbEsp5, Count5, 0, Liber);
         statcoll_AddObs (res->Bas->sVal1, X2);
      }

      /* Compute the normalized observed number for each pair */
      for (i = 0; i <= L; i++) {
         for (j = 0; j <= L; j++) {
            NbEsp = nLR * Prob[i] * Prob[j];
            Var = NbEsp * (1.0 - Prob[i] * Prob[j]);
            if (Var <= 0.0) {
               /* This case will occur when NbEsp = 0; if so, fail the test */
               /* when res->Counters[i][j] != 0 */
               res->ZCounters[i][j] = (res->Counters[i][j] - NbEsp) * 1.E100;
            } else {
               res->ZCounters[i][j] = (res->Counters[i][j] - NbEsp) /
                  sqrt (Var);
            }
         }
      }
      if (sstring_Counters)
         /* Print the matrix of counters */
         tables_WriteMatrixL (res->Counters, 0, L, 0, L, 8, res->Style,
            "res->Counters, the number of pairs [0..L, 0..L]");
      if (swrite_Counters)
         /* Print the matrix of normalized counters */
         tables_WriteMatrixD (res->ZCounters, 0, L, 0, L, 12, 4,
            res->Style, "res->ZCounters, the normalized counters");

      /* These blocks are sub-matrices symmetrically placed with respect to */
      /* the diagonals in the matrix of the number of pairs [Xi, X(i+1)]. */
      /* For those, we shall apply a chi-square test for the total number */
      /* in the diagonal sub-matrices. */
      CountBlocks (res, L, d);

      if (swrite_Counters)
         WriteBlocs (res, d);

      for (i = 1; i <= d; i++) {
         double NumExp[3];
         long Count[3];
         Sum = 0.0;
         for (j = 0; j <= (L + 1) / 2 - i; j++)
            Sum += Prob[j];
         /* Probability of a sub-matrix block */
         Sum *= Sum;
         /* Average number for each of the 2 blocks */
         NbMoyen = Sum * nLR * 2.0;
         NumExp[0] = NumExp[1] = NbMoyen;
         if (2.0*NbMoyen < gofs_MinExpected)
              printf ("******* sample too small for chi-square for d = %d\n", i);
         NumExp[2] = nLR - 2.0*NbMoyen;
         Count[0] = res->XD[i][0];
         Count[1] = res->XD[i][1];
         Count[2] = n - Count[0] - Count[1];
         X2 = gofs_Chi2 (NumExp, Count, 0, 2);
         statcoll_AddObs (res->Block[i]->sVal1, X2);
      }
   }

   for (i = 1; i <= d; i++) {
      /* Degrees of freedom */
      if ((L & 1) && i == 1)
         Liberte = 1;
      else
         Liberte = 2;

      V[0] = Liberte;
      gofw_ActiveTests2 (res->Block[i]->sVal1->V, res->Block[i]->pVal1->V, N,
         wdist_ChiSquare, V, res->Block[i]->sVal2, res->Block[i]->pVal2);
      res->Block[i]->pVal1->NObs = N;
      Sum = N * statcoll_Average (res->Block[i]->sVal1);
      res->Block[i]->sVal2[gofw_Sum] = Sum;
      res->Block[i]->pVal2[gofw_Sum] = fbar_ChiSquare2 (N*Liberte, 12, Sum);

      if (swrite_Basic) {
         printf ("\nDiagonal blocks with d = %2d", i);
         swrite_AddStrChi (str, LEN2, Liberte);
         gofw_WriteActiveTests2 (N, res->Block[i]->sVal2,
                                 res->Block[i]->pVal2, str);
         swrite_Chi2SumTestb (N, res->Block[i]->sVal2[gofw_Sum],
                                 res->Block[i]->pVal2[gofw_Sum], Liberte);
         if (swrite_Collectors) {
            strncpy (chaine, res->Block[i]->sVal1->Desc, (size_t) LEN1);
            strncat (chaine, ":   a chi2 with ", (size_t) LEN1);
            sprintf (str, "%1d", Liberte);
            strncat (chaine, str, (size_t) LEN2);
            strncat (chaine, " degrees of freedom", (size_t) LEN1);
            statcoll_SetDesc (res->Block[i]->sVal1, chaine);
            statcoll_Write (res->Block[i]->sVal1, 5, 14, 4, 3);
         }
      }
   }

   if (Liber > 0) {
      V[0] = Liber;
      gofw_ActiveTests2 (res->Bas->sVal1->V, res->Bas->pVal1->V, N,
         wdist_ChiSquare, V, res->Bas->sVal2, res->Bas->pVal2);
      res->Bas->pVal1->NObs = N;
      Sum = N * statcoll_Average (res->Bas->sVal1);
      res->Bas->sVal2[gofw_Sum] = Sum;
      res->Bas->pVal2[gofw_Sum] = fbar_ChiSquare2 (N*Liber, 12, Sum);

      if (swrite_Basic) {
         printf ("%s", res->Bas->sVal1->Desc);
         swrite_AddStrChi (str, LEN2, Liber);
         gofw_WriteActiveTests2 (N, res->Bas->sVal2, res->Bas->pVal2, str);
         swrite_Chi2SumTestb (N, res->Bas->sVal2[gofw_Sum],
                                 res->Bas->pVal2[gofw_Sum], Liber);
         if (swrite_Collectors) {
            strncpy (chaine, res->Bas->sVal1->Desc, (size_t) LEN1);
            strncat (chaine, ":   a ChiSquare with ", (size_t) LEN1);
            sprintf (str, "%1d", Liber);
            strncat (chaine, str, (size_t) LEN2);
            strncat (chaine, " degrees of freedom", (size_t) LEN1);
            statcoll_SetDesc (res->Bas->sVal1, chaine);
            statcoll_Write (res->Bas->sVal1, 5, 14, 4, 3);
         }
      }
   } else {
      /* for the module tvaria */
      res->Bas->pVal2[gofw_Mean] = -1.0;
      if (d < 1)
         util_Warning (1,
            "n is too small:   ChiSquare with 0 degree of freedom");
   }

   if (swrite_Basic)
      swrite_Final (gen, Timer);

   util_Free (Prob);
   util_Free (NbEsp5);
   util_Free (Count5);
   if (localRes)
      sstring_DeleteRes (res);
   chrono_Delete (Timer);
}
