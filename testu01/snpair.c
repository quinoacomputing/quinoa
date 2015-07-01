/*************************************************************************\
 *
 * Package:        TestU01
 * File:           snpair.c
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
#include "num.h"
#include "num2.h"

#include "snpair.h"
#include "swrite.h"
#include "unif01.h"

#include "statcoll.h"
#include "fdist.h"
#include "fbar.h"
#include "fmass.h"
#include "gofs.h"
#include "gofw.h"

#include <math.h>
#include <limits.h>
#include <float.h>
#include <stddef.h>
#include <string.h>

#undef DEBUG
#ifdef DEBUG

#include <stdio.h>
/* Prints all the points */
#define TRACEP(n, name) { \
   FILE *f; \
   long i, j; \
   f = util_Fopen (name, "w"); \
   fprintf (f, "------------------------\n"); \
   for (i = 1; i <= n; i++) { \
      for (j = 1; j <= kk; j++) \
	fprintf (f, "%f    ", res->Points[1][i][j]); \
      fprintf (f, "\n"); \
   } \
}
#endif




/*---------------------------- extern variables ---------------------------*/

snpair_Envir snpair_env = {
   20, 30, 30, 1000
};

/* For now, we do not use it. Is used in the t modules */
long snpair_MaxNumPoints = LONG_MAX;

lebool snpair_TimeBB = FALSE;
lebool snpair_mNP2S_Flag = TRUE;




/*---------------------------- static variables ---------------------------*/

typedef struct {

   int L1;                     /* Change coordinate at each L1 recur-  */
                               /* sion level in FindClosePairs         */
   int L2;                     /* Same thing for CheckBoundary         */
   int kk;                     /* = k = Dimension                      */
   int pp;                     /* = p (L_p norm); p = 0 is sup norm    */
   int mm;                     /* = m = (number of kept distances)     */
   int mcd;                    /* Dimension of CloseDist[]             */
   double dlim1;               /* (m \mu2)^{1/k}                       */
   double dlim1p;              /* dlim1^p (when p > 0)                 */
   double dlim;                /* = max (dlim1, CloseDist[m]); search  */
                               /* for points at distance < dlim        */
   double dlimp;               /* = dlim^p (when p > 0)                */
   double pLR;                 /* = p                                  */
   double Invp;                /* 1/p; ( = 1 if p = 0)                 */
   int Maxnp;                  /* Max level of recursion               */
   lebool Torus;              /* TRUE if in Torus; FALSE if in cube   */
   lebool BBFlag;             /* TRUE if BickelBreimann test          */
   wdist_CFUNC FDistBB;        /* BickelBreimann CDF                   */

/* The largest distance for snpair_DistanceCPBitM: the largest number of
   equal bits for all components of a pair of points (all components of
   the pair have at least YLim identical bits) amongst all pairs. */
   int YLim;

} WorkType;


/*------------------ Module variables for G+, G-, H+, H-  -----------------*/
/*
 * To compute the jumps in FDistGPlus, FDistGMinus, FDistHPlus, FDistHMinus
 */


#if 0
static double t0;

#define xinf 1.0e50

static lebool GPlusFlag  = FALSE;
static lebool GMinusFlag = FALSE;
static lebool HPlusFlag  = FALSE;
static lebool HMinusFlag = FALSE;

static double GPlust0  = -xinf;    /* = t0 for GPlus  */
static double GMinust0 = -xinf;    /* = t0 for GMinus */
static double HPlust0  = -xinf;    /* = t0 for HPlus  */
static double HMinust0 = -xinf;    /* = t0 for HMinus */
static double GPlust1  = -xinf;    /* = t1 for GPlus  */
static double GMinust1 = -xinf;    /* = t1 for GMinus */
static double HPlust1  = -xinf;    /* = t1 for HPlus  */
static double HMinust1 = -xinf;    /* = t1 for HMinus */

static double *GPlusJumpX;         /* Position x of the jumps of G+       */
static double *GPlusJumpYBottom;   /* Value y left of the jumps of G+     */
static double *GPlusJumpYTop;      /* Value y right of the jumps of G+    */

static double *GMinusJumpX;        /* Similarly for G- */
static double *GMinusJumpYBottom;
static double *GMinusJumpYTop;

static double *HPlusJumpX;         /* Similarly for H+ */
static double *HPlusJumpYBottom;
static double *HPlusJumpYTop;

static double *HMinusJumpX;        /* Similarly for H- */
static double *HMinusJumpYBottom;
static double *HMinusJumpYTop;

static int GPlusNJump;             /* The number of jumps of G+ */
static int GMinusNJump;            /* The number of jumps of G- */
static int HPlusNJump;             /* The number of jumps of H+ */
static int HMinusNJump;            /* The number of jumps of H- */

#endif

/*------------------ Module variables for Bickel-Breiman test  ------------*/

static double BB2[132];
static double BB3[132];
static double BB4[43];
static double BB5[22];






/*-------------------------------- functions ------------------------------*/

static void CopyPoints (snpair_PointType A[], snpair_PointType B[], long r,
   long s)
/*
 * Copies A[r..s] into B[r..s]
 */
{
   long i;
   for (i = r; i <= s; i++)
      B[i] = A[i];
}


/*=========================================================================*/

void snpair_QuickSort (snpair_PointType A[], long l, long r, int c)
/*
 * Sort points of indices l to r using coordinate c as key. Exchange
 * pointers instead of points.
 */
{
   long j;
   long i;
   double pivot;
   snpair_PointType vec;
   i = l;
   j = r;
   pivot = A[(l + r) / 2][c];
   do {
      while (A[i][c] < pivot)
         ++i;
      while (pivot < A[j][c])
         --j;
      if (i <= j) {
         vec = A[i];
         A[i] = A[j];
         A[j] = vec;
         ++i;
         --j;
      }
   } while (i <= j);
   if (l < j)
      snpair_QuickSort (A, l, j, c);
   if (i < r)
      snpair_QuickSort (A, i, r, c);
}


/*=========================================================================*/

void snpair_DistanceCP (snpair_Res * res, snpair_PointType P1,
   snpair_PointType P2)
/*
 * For ClosePairs, checks if the distance between P1 and P2 is < dlim.
 * If so, updates dlim, dlimp, and adds the new distance in the list
 * of shortest distances.
 */
{
   int i;
   double temp;
   double dist;
   double distp = 0.0;
   WorkType *work = res->work;

   for (i = 1; i <= work->kk; i++) {
      temp = P1[i] - P2[i];
      if (temp < 0.0)
         temp = -temp;
      if (work->Torus && temp > 0.5)
         temp = 1.0 - temp;
      if (work->pp == 0) {
         if (temp > distp)
            distp = temp;
      } else if (work->pp == 1)
         distp += temp;
      else if (work->pp == 2)
         distp += temp * temp;
      else
         distp += pow (temp, work->pLR);
      if (distp >= work->dlimp)
         return;
   }

#define NUM_JUMPS_LIM 50000
/* I put this arbitrary limit on the size of res->CloseDist because for bad
   generators, there will be many pairs of points with 0 distances between
   them, and res->CloseDist will otherwise eat up the whole memory. 
   If the parameters N and m in ClosePairs should be such that 
   N*m > NUM_JUMPS_LIM, then this limit will have to be increased. (RS) */

   if (distp < work->dlimp) {
      if (work->pp <= 1)
         dist = distp;
      else if (work->pp == 2)
         dist = sqrt (distp);
      else
         dist = pow (distp, work->Invp);
      if ((res->NumClose < work->mm ||
           res->CloseDist[res->NumClose] < work->dlim1) &&
          (res->NumClose < NUM_JUMPS_LIM)) {
         ++res->NumClose;         /* Complete the list of close pairs */
         if (res->NumClose >= work->mcd) {
            double *A;
            work->mcd *= 2;
            A = util_Realloc (res->CloseDist, (work->mcd + 1)*sizeof (double));
            if (A == NULL) {
               util_Warning (1, "Cannot realloc res->CloseDist");
            } else
            res->CloseDist = A;
         }
         util_Warning ((res->NumClose >= NUM_JUMPS_LIM) && swrite_Basic,
              "res->NumClose > 50000");
      }

      /* Insert the new distance in the sorted list */
      i = res->NumClose;
      while (i > 1 && dist < res->CloseDist[i - 1]) {
         --i;
         res->CloseDist[i + 1] = res->CloseDist[i];
      }
      res->CloseDist[i] = dist;

      if (res->NumClose == work->mm && res->CloseDist[work->mm] < work->dlim
         && work->dlim1 < work->dlim) {
         work->dlim = res->CloseDist[work->mm];
         if (work->dlim < work->dlim1) {
            work->dlim = work->dlim1;
            work->dlimp = work->dlim1p;
         } else if (work->pp <= 1)
            work->dlimp = work->dlim;
         else if (work->pp == 2)
            work->dlimp = work->dlim * work->dlim;
         else
            work->dlimp = pow (work->dlim, work->pLR);
      }
   }
}


/*=========================================================================*/

void snpair_DistanceBB (snpair_Res * res, snpair_PointType P1,
   snpair_PointType P2)
/*
 * For Bickel-Breiman, checks whether the distance between P1 and P2, to the
 * power p, is less than P1[0] or P2[0]. If so, updates these values.
 */
{
   int i;
   double bound;
   double temp;
   double distp;
   WorkType *work = res->work;

   if (P2[0] > P1[0])
      bound = P2[0];
   else
      bound = P1[0];

   distp = 0.0;
   for (i = 1; i <= work->kk; i++) {
      temp = P1[i] - P2[i];
      if (temp < 0.0)
         temp = -temp;
      if (work->Torus && temp > 0.5)
         temp = 1.0 - temp;
      if (work->pp == 1)
         distp += temp;
      else if (work->pp == 2)
         distp += temp * temp;
      else if (work->pp == 0) {
         if (temp > distp)
            distp = temp;
      } else
         distp += pow (temp, work->pLR);
      if (distp >= bound)
         return;
   }

   if (distp < P1[0])
      P1[0] = distp;
   if (distp < P2[0])
      P2[0] = distp;
}


/*=========================================================================*/

void snpair_VerifPairs0 (snpair_Res * res, snpair_PointType A[], long r,
   long s, int junk1, int junk2)
/*
 * Compute the distance between all pairs of points with indices in the
 * interval [r..s] for the array A; updates the distances for BB. We assume
 * the points are sorted with respect to coordinate c.
 */
{
   long i, j;
   for (i = r; i < s; i++) {
      for (j = i + 1; j <= s; j++) {
         res->Distance (res, A[i], A[j]);
      }
   }
}


/*=========================================================================*/

void snpair_VerifPairs1 (snpair_Res * res, snpair_PointType A[], long r,
   long s, int np, int c)
/*
 * Compute the distance between all pairs of points with indices in the
 * interval [r..s] for the array A; updates dlim and dlimp if necessary.
 * We assume the points are sorted with respect to coordinate c.
 */
{
   long i, j;
   double high;
   WorkType *work = res->work;

   util_Assert (np <= work->Maxnp,
      "Calling snpair_VerifPairs1 with np > Maxnp");
   for (i = r; i <= s; i++) {
      /* util_Assert (r <= s, "Calling snpair_VerifPairs1 with r > s"); */

      /* Consider only points at distance <= dlim from A[i] with respect to
         coordinate c */
      high = A[i][c] + work->dlim;
      j = i + 1;
      while (j <= s && A[j][c] < high) {
         res->Distance (res, A[i], A[j]);
         ++j;
      }
      if (j > s && work->Torus && np <= work->kk) {
         high -= 1.0;
         j = r;
         while (j < i && A[j][c] < high) {
            /* util_Assert (i != j, "Calling distance with i=j in
               snpair_VerifPairs1"); */
            res->Distance (res, A[i], A[j]);
            ++j;
         }
      }
   }
}


/*=========================================================================*/

static void dlimSlice (
   snpair_Res *res, 
   snpair_PointType A[],
   long *r,
   long *imed,
   long *jmed,
   long *s,
   int c,
   lebool Tor
   )
/*
 * Let E1 = A [*r..*imed] and E2 = A [*jmed..*s] be two sets of points sorted
 * with respect to the c-th coordinate.
 * If Torus = FALSE, will reduce E1 to a slice of width dlim (in the direction
 * c) of the leftmost point of E2; and similarly for E2 with respect to the
 * rightmost point of E1 (always with respect to the c-th coordinate).
 * There will be thus a slice on each side of the boundary line.
 * If Torus = TRUE, we consider rather the points to the left of
 * (E1 + 1.0) as being close to those to the right of E2. This procedure
 * tries to decrease imed and to increase jmed.
 */
{
   long i;
   double temp;
   WorkType *work = res->work;

   if (*r > *imed || *jmed > *s)
      return;

#ifdef DEBUG
   printf ("ENTER dlimslice ");
   num_WriteD (work->dlim, 10, 5, 1);
   printf (" %5ld %5ld %5ld %5ld %3d  ", *r, *imed, *jmed, *s, c);
   util_WriteBool (Tor, 5);
   printf ("\n");
#endif

   if (Tor) {
      temp = A[*s][c] - 1.0;
      i = *r;
      while (i <= *imed && A[i][c] - temp < work->dlim)
         ++i;
      *imed = i - 1;
      temp = A[*r][c] + 1.0;
      i = *s;
      while (i >= *jmed && temp - A[i][c] < work->dlim)
         --i;
      *jmed = i + 1;

   } else {
      temp = A[*jmed][c];
      i = *imed;
      while (i >= *r && temp - A[i][c] < work->dlim)
         --i;
      *r = i + 1;
      temp = A[*imed][c];
      i = *jmed;
      while (i <= *s && A[i][c] - temp < work->dlim)
         ++i;
      *s = i - 1;
   }

#ifdef DEBUG
   printf ("EXIT dlimslice            ");
   printf (" %5ld %5ld %5ld %5ld\n", *r, *imed, *jmed, *s);
#endif
}


/*=========================================================================*/

void snpair_MiniProc0 (snpair_Res * res, snpair_PointType T[], long r,
   long s, long u, long v, int junk1, int junk2)
/*
 * Call "res->Distance" for each point of T[r..s] with each point
 * of T[u..v].
 */
{
   long i, j;

   for (i = r; i <= s; i++)
      for (j = u; j <= v; j++)
         res->Distance (res, T[i], T[j]);
}


/*=========================================================================*/

void snpair_MiniProc1 (snpair_Res * res, snpair_PointType T[], long r,
   long s, long u, long v, int np, int c)
/* 
 * Compute the distance between each point of the set E1 = T[r..s] and 
 * those of E2 = T[u..v]. We consider only the pairs for which the  
 * differences between the c-th coordinate is <= dlim.                  
 * Assume that the points of E1 and those of E2 are sorted with respect to
 * coordinate c. If we have the case of snpair_DistanceCP, updates dlim and
 * dlimp whenever a shorter value is found.
 */
{
   long inf;
   long l, k, j, i;
   double high, low;
   WorkType *work = res->work;

#ifdef DEBUG
   printf ("ENTER MiniProc1 ");
   num_WriteD (work->dlim, 10, 5, 1);
   printf (" %5ld %5ld %5ld %5ld    %3d\n", r, s, u, v, np);
   if (v < r)
      printf ("MiniProc1 with v < r !!!\n");
   util_Assert (np <= work->Maxnp, "MiniProc1:  np > Maxnp");
   util_Assert ((s < u) || (v < r),
      "MiniProc1:   Overlap of E1 and E2 dans MiniProc");
#endif

   if (s < r || v < u)
      return;
   inf = u;
   /* sup = v; */

#ifdef DEBUG
   for (i = r; i < s; i++)
      util_Assert (T[i][c] <= T[i + 1][c], "Wrong order in MiniProc1");
   for (i = u; i < v; i++)
      util_Assert (T[i][c] <= T[i + 1][c], "Wrong order in MiniProc1");
#endif


   for (i = r; i <= s; i++) {
      low = T[i][c] - work->dlim;
      high = low + 2.0 * work->dlim;
      /* consider only points at distance <= dlim of T[i] w.r. to coord. c */
      while (inf <= v && T[inf][c] <= low)
         ++inf;
      j = inf;
      while (j <= v && T[j][c] < high) {
         res->Distance (res, T[i], T[j]);
         ++j;
      }
      if (work->Torus) {          /* AND (np <= kk): does not work with this 
                                     cond. */

         /* Search for close points in the torus */
         low += 1.0;
         high -= 1.0;
         k = u;
         l = v;
         while (k <= v && T[k][c] < high) {
            res->Distance (res, T[i], T[k]);
            ++k;
         }
         while (l >= u && T[l][c] > low) {
            res->Distance (res, T[i], T[l]);
            --l;
         }
      }
   }
}


/*=========================================================================*/

void snpair_CheckBoundary (snpair_Res * res, long r, long s, long u, long v,
   int nr, int nrb, int np, int c)
/*
 * Compute the minimal distance between the points of the sets E1 = A[r..s]
 * and E2 = A[u..v].
 * nrb is the recursion level of the calls of snpair_CheckBoundary.
 * A always stands for Points[np], sorted with respect to coordinate c.
 */
{
   long jmed2, imed2;
   long jmed, imed;
   long nextc;
   lebool newc;
   snpair_PointTableType B, A;
   WorkType *work = res->work;

#ifdef DEBUG
   printf ("CheckBoundary:   ");
   printf (" %8ld %8ld %8ld %8ld %3d %3d %3d %3d\n",
      r, s, u, v, nr, nrb, np, c);
#endif

   if (r > s || u > v)
      return;

   util_Assert (np <= work->Maxnp, "np > Maxnp in snpair_CheckBoundary");
   A = res->Points[np];
   newc = ((nrb - 1) % work->L2) == 0;
   if (newc && np < work->Maxnp) {
      B = res->Points[np + 1];
      ++np;
      if (c < work->kk)
         nextc = c + 1;
      else
         nextc = 1;

      /* Copy the remaining points in the table of level np+1, then */
      /* sort with respect to coordinate nextc. */
      CopyPoints (A, B, r, s);
      CopyPoints (A, B, u, v);
      snpair_QuickSort (B, r, s, nextc);
      snpair_QuickSort (B, u, v, nextc);

   } else {
      nextc = c;
      B = A;
   }

   if ((nrb >= work->kk || s - r < snpair_env.Seuil2)
      || v - u < snpair_env.Seuil2) {
      /* Max recursion or small sets of points */
      res->MiniProc (res, B, r, s, u, v, np, nextc);
      return;
   }

   /* We halve each set of points, and we check each half on one side */
   /* with each half of the other side. */
   imed = (r + s) / 2;
   jmed = (u + v) / 2;

   /* Check the halves which are face to face */
   snpair_CheckBoundary (res, r, imed, u, jmed, nr + 1, nrb + 1, np, nextc);
   snpair_CheckBoundary (res, imed + 1, s, jmed + 1, v, nr + 1, nrb + 1, np,
                         nextc);

   /* Check the other (crossed) halves */
   if ((work->Torus && np <= work->kk) && newc) {
      imed2 = imed;
      jmed2 = jmed + 1;
      dlimSlice (res, B, &r, &imed2, &jmed2, &v, nextc, TRUE);
      snpair_CheckBoundary (res, r, imed2, jmed2, v, nr + 1, nrb + 1, np,
         nextc);

      imed2 = imed + 1;
      jmed2 = jmed;
      dlimSlice (res, B, &u, &jmed2, &imed2, &s, nextc, TRUE);
      snpair_CheckBoundary (res, u, jmed2, imed2, s, nr + 1, nrb + 1, np,
         nextc);
   }

   jmed2 = jmed + 1;
   imed2 = imed + 1;
   if (newc)
      dlimSlice (res, B, &r, &imed, &jmed2, &v, nextc, FALSE);
   snpair_CheckBoundary (res, r, imed, jmed + 1, v, nr + 1, nrb + 1, np,
      nextc);
   if (newc)
      dlimSlice (res, B, &u, &jmed, &imed2, &s, nextc, FALSE);
   snpair_CheckBoundary (res, u, jmed, imed + 1, s, nr + 1, nrb + 1, np,
      nextc);
}


/*=========================================================================*/

static void Setdlim (snpair_Res *res, snpair_PointType A[], long r, long s)
/* 
 * Used in snpair_BickelBreiman to update dlim
 */
{
   long i;
   WorkType *work = res->work;

   work->dlimp = 0.0;
   for (i = r; i <= s; i++) {
      if (A[i][0] > work->dlimp)
         work->dlimp = A[i][0];
   }
   if (work->pp == 0 || work->pp == 1)
      work->dlim = work->dlimp;
   else if (work->pp == 2)
      work->dlim = sqrt (work->dlimp);
   else
      work->dlim = pow (work->dlimp, work->Invp);

#ifdef DEBUG
   printf ("Setdlim: ");
   num_WriteD (work->dlim, 10, 5, 1);
   printf ("\n");
#endif
}


/*=========================================================================*/

void snpair_FindClosePairs (snpair_Res * res, long r, long s,
   int nr, int np, int c)
/*
 * Checks whether the minimal distance between the 2 nearest points,
 * amongst those with indices [r..s] in the table of level np, is < dlim.  
 * If so, updates dlim and dlimp. A and B are always Points [np] and
 * Points [np+1].                   
 * Assume that A = Point [np] is sorted with respect to coordinate c.    
 */
{
   long jmed2;
   long imed2;
   long imed;
   long nextc;                    /* Next coordinate to be used */
   snpair_PointTableType B;
   snpair_PointTableType A;
   WorkType *work = res->work;

   /* IF (((nr-1) MOD L) = 0) THEN newc := TRUE ELSE newc := FALSE END; */
#ifdef DEBUG
   printf ("FindClosePairs:   ");
   printf (" %8ld %8ld %3d %3d %3d\n", r, s, nr, np, c);
#endif

   util_Assert (np <= work->Maxnp, "np > Maxnp in snpair_FindClosePairs");
   A = res->Points[np];
   if (s - r < snpair_env.Seuil1) {
      res->VerifPairs (res, A, r, s, np, c);
      /* Here we are finished */
      return;
   }

   /* We divide the points in 2 approximately equal sets E1 and E2; */
   /* then recursions upon the sets E1 and E2. */
   imed = (r + s) / 2;
   if (nr % work->L1 == 0 && np < work->Maxnp && np < work->kk) {
      /****** Condition np < kk is temporary... ********/
      util_Assert (np == 1 + (nr - 1) / work->L1,
         "Bad np in snpair_FindClosePairs");
      /* IF np >= Maxnp THEN VerifPairs (A^, r, s, nr, c); RETURN END; */

      /* We shall now increase np and switch coordinate.  */
      /* Copy the points in a new table for the next level. */
      B = res->Points[np + 1];
      CopyPoints (A, B, r, s);
      if (c < work->kk)
         nextc = c + 1;
      else
         nextc = 1;
      util_Assert (nextc == 1 + (np % work->kk),
         "Bad nextc dans snpair_FindClosePairs");
      snpair_QuickSort (B, r, imed, nextc);
      snpair_QuickSort (B, imed + 1, s, nextc);
      snpair_FindClosePairs (res, r, imed, nr + 1, np + 1, nextc);
      snpair_FindClosePairs (res, imed + 1, s, nr + 1, np + 1, nextc);

   } else {
      snpair_FindClosePairs (res, r, imed, nr + 1, np, c);
      snpair_FindClosePairs (res, imed + 1, s, nr + 1, np, c);
   }

   /* It remains to check the boundary between E1 and E2. */
   if (work->kk == 1) {
      res->Distance (res, A[imed], A[imed + 1]);
      if (work->Torus)
         res->Distance (res, A[r], A[s]);
      return;
   }

   /* Bring m and n closer in order to sqeeze only the points which could */
   /* be at a distance less than dlim from the median. */
   if (work->BBFlag)
      Setdlim (res, A, r, s);
   if (work->Torus && np <= work->kk && (nr - 1) % work->L1 == 0) {
      imed2 = imed;
      jmed2 = imed + 1;
      dlimSlice (res, A, &r, &imed2, &jmed2, &s, c, TRUE);
      snpair_CheckBoundary (res, r, imed2, jmed2, s, nr, 1, np, c);
   }

   jmed2 = imed + 1;
   dlimSlice (res, A, &r, &imed, &jmed2, &s, c, FALSE);
   snpair_CheckBoundary (res, r, imed, jmed2, s, nr, 1, np, c);
}


/*=========================================================================*/

#if 0
#define Epsilon 1.0e-10

static double Probsup (double b, double c, double x)
{
   int jsup;
   int msup;
   int m;
   int j;
   double comb;
   double mLR;
   double jLR;
   double mFact;
   double Sum2;
   double Sum;
   double Previous;

   if (x < 0.0)
      return 0.0;

   Sum = 0.0;
   mFact = 1.0;
   if (x <= 0.0) {
      msup = c * b;
      if (msup > 100) {
         msup = 100;
         util_Warning (TRUE, "Probsup: msup > 100. Reset to 100");
      }
      for (m = 1; m <= msup; m++) {
         mLR = m;
         mFact *= mLR;
         Sum += (pow (b, mLR) / mFact) * (1.0 - mLR / (c * b));
      }
      if (msup >= 0)
         Sum += 1.0;
      return Sum * exp (-b);
   }

   Previous = -1.0;
   m = 1;
   msup = c * b + x;
   if (msup > 100) {
      msup = 100;
      util_Warning (TRUE, "Probsup: msup > 100.  Reset to 100");
   }
   while (m <= msup && Sum - Previous > Epsilon) {
      Previous = Sum;
      Sum2 = 0.0;
      mLR = m;
      mFact *= mLR;
      jsup = x;
      if (jsup > m) {
         jsup = m;
         util_Warning (TRUE, "Probsup: jsup > m.  Reset to m");
      }
      comb = 1.0;
      for (j = 0; j <= jsup; j++) {
         jLR = j;
         Sum2 += comb * pow (jLR - x, jLR) *
            pow (c * b + x - jLR, mLR - jLR - 1.0);
         comb *= (mLR - jLR) / (jLR + 1.0);
      }
      Sum += (c * b + x - mLR) / (mFact * pow (c, mLR)) * Sum2;
      ++m;
   }
   if (msup >= 0)
      Sum += 1.0;
   return Sum * exp (-b);
}


/*=========================================================================*/

static double Probinf (double b, double c, double x)
{
   int msup;
   int m;
   double mLR;
   double mFact;
   double Sum;

   if (x >= 0.0)
      return 1.0;

   msup = c * b + x;
   if (msup > 100) {
      msup = 100;
      util_Warning (TRUE, "Probinf: msup > 100. Reset to 100");
   }
   Sum = 0.0;
   mFact = 1.0;
   for (m = 1; m <= msup; m++) {
      mLR = m;
      mFact *= mLR;
      Sum += exp (-(mLR - x) / c) *
         (pow (mLR - x, mLR - 1.0) / (pow (c, mLR) * mFact));
   }
   if (msup >= 0)
      Sum = -(x * Sum) + exp (x / c);
   return Sum;
}


/*=========================================================================*/

static double FDistGPlus (double Bidon, double c)
/*
 * Obsolete. This discontinuous distribution uses a complicated statistic
 * and did not seem sensitive. We don't use it anymore.

 *  See the reference
 *  P. L'Ecuyer, J.-F. Cordeau, and R. Simard,  "Close-Point Spatial Tests
 *     and their Application to Random Number Generators",
 *     Operations Research, 48, 2 (2000), 308--317
 *
 */
{
   int lSup;
   int l;
   double lLR;
   double lFact;
   double Previous;
   double Sum;

   if ((!GPlusFlag || GPlust0 != t0) || GPlust1 != t1) {
      GPlust0 = t0;
      GPlust1 = t1;
      GPlusFlag = TRUE;
      /* 
         fdist_FindJumps (W, Detail);

      FindJumpsKnown (Bidon, FDistGPlus, GPlust0, 20.0, 0.00001, &GPlusNJump,
      &GPlusJumpX, &GPlusJumpYBottom, &GPlusJumpYTop);*/
   }

   if (c < 0.0)
      return 0.0;
   l = 1;
   Sum = 0.0;
   Previous = -1.0;
   lFact = 1.0;
   lSup = GPlust0 * c;
   while (l <= lSup && Sum - Previous > Epsilon) {
      Previous = Sum;
      lLR = l;
      lFact *= lLR;
      Sum += pow (GPlust0, lLR) / lFact *
                Probsup (GPlust1 - GPlust0, c, GPlust0 * c - lLR);
      ++l;
   }
   Sum += Probsup (GPlust1 - GPlust0, c, GPlust0 * c);
   return Sum * exp (-GPlust0);
}


/*=========================================================================*/

static double FDistGMinus (double Bidon, double c)
/*
 * Obsolete. This discontinuous distribution uses a complicated statistic
 * and did not seem sensitive. We don't use it anymore.

 *  See the reference
 *  P. L'Ecuyer, J.-F. Cordeau, and R. Simard,  "Close-Point Spatial Tests
 *     and their Application to Random Number Generators",
 *     Operations Research, 48, 2 (2000), 308--317
 *
 */
{
   int l;
   double lLR;
   double lFact;
   double Previous;
   double Sum;

   if ((!GMinusFlag || GMinust0 != t0) || GMinust1 != t1) {
      GMinust0 = t0;
      GMinust1 = t1;
      GMinusFlag = TRUE;
      /*
      FindJumpsKnown (Bidon, FDistGMinus, GMinust1, 20.0, 0.00001,
         &GMinusNJump, &GMinusJumpX, &GMinusJumpYBottom, &GMinusJumpYTop);
      */
   }

   if (c < 0.0)
      return 0.0;
   l = 1;
   Sum = 0.0;
   Previous = -1.0;
   lFact = 1.0;
   while (Sum - Previous > Epsilon) {
      Previous = Sum;
      lLR = l;
      lFact *= lLR;
      Sum += pow (GMinust0, lLR) / lFact *
             Probinf (GMinust1 - GMinust0, c, GMinust0 * c - lLR);
      ++l;
   }
   Sum += Probinf (GMinust1 - GMinust0, c, GMinust0 * c);
   return Sum * exp (-GMinust0);
}


/*=========================================================================*/

static double FDistHPlus (double b, double x)
/*
 * Obsolete. This discontinuous distribution uses a complicated statistic
 * and did not seem sensitive. We don't use it anymore.

 *  See the reference
 *  P. L'Ecuyer, J.-F. Cordeau, and R. Simard,  "Close-Point Spatial Tests
 *     and their Application to Random Number Generators",
 *     Operations Research, 48, 2 (2000), 308--317
 *
 */
{
   int msup;
   int jsup;
   int m;
   int j;
   double comb;
   double mLR;
   double jLR;
   double mFact;
   double Sum2;
   double Sum;
   double Previous;

   if ((!HPlusFlag || HPlust0 != t0) || HPlust1 != t1) {
      /* this function has a single jump at x = 0 */
      HPlust0 = t0;
      HPlust1 = t1;
      HPlusFlag = TRUE;/*
      FindJumpsKnown (b, FDistHPlus, 1.0, 0.00001, 1.E-6, &HPlusNJump,
      &HPlusJumpX, &HPlusJumpYBottom, &HPlusJumpYTop);*/
   }

   if (x < 0.0)
      return 0.0;

   Sum = 0.0;
   mFact = 1.0;
   if (x <= 0.0) {
      msup = b;
      if (msup > 100)
         msup = 100;
      for (m = 1; m <= msup; m++) {
         mLR = m;
         mFact *= mLR;
         Sum += pow (b, mLR) / mFact * (1.0 - mLR / b);
      }
      if (msup >= 0)
         Sum += 1.0;
      return Sum * exp (-b);
   }

   Previous = -1.0;
   m = 1;
   msup = b + x;
   if (msup > 100)
      msup = 100;
   while (m <= msup && Sum - Previous > Epsilon) {
      Previous = Sum;
      Sum2 = 0.0;
      mLR = m;
      mFact *= mLR;
      jsup = x;
      if (jsup > m)
         jsup = m;
      comb = 1.0;
      for (j = 0; j <= jsup; j++) {
	 jLR = j;
	 Sum2 += comb*pow (jLR - x, jLR)*pow (b + x - jLR, mLR - jLR - 1.0);
	 comb *= (mLR - jLR) / (jLR + 1.0);
      }
      Sum += (b + x - mLR) / mFact * Sum2;
      ++m;
   }
   if (msup >= 0)
      Sum += 1.0;
   return Sum * exp (-b);
}


/*=========================================================================*/

static double FDistHMinus (double b, double x)
/*
 * Obsolete. This discontinuous distribution uses a complicated statistic
 * and did not seem sensitive. We don't use it anymore.

 *  See the reference
 *  P. L'Ecuyer, J.-F. Cordeau, and R. Simard,  "Close-Point Spatial Tests
 *     and their Application to Random Number Generators",
 *     Operations Research, 48, 2 (2000), 308--317
 *
 */
{
   int msup;
   int m;
   double mLR;
   double mFact;
   double Sum;

   if ((!HMinusFlag || HMinust0 != t0) || HMinust1 != t1) {
      HMinust0 = t0;
      HMinust1 = t1;
      HMinusFlag = TRUE;
      /*
      FindJumpsKnown (b, FDistHMinus, -b, 0.0001, 0.00001, &HMinusNJump,
      &HMinusJumpX, &HMinusJumpYBottom, &HMinusJumpYTop); */
   }

   if (x >= 0.0)
      return 1.0;

   msup = b + x;
   if (msup > 100)
      msup = 100;
   Sum = 0.0;
   mFact = 1.0;
   for (m = 1; m <= msup; m++) {
      mLR = m;
      mFact *= mLR;
      Sum += exp (-(mLR - x)) * (pow (mLR - x, mLR - 1.0) / mFact);
   }

   if (msup >= 0)
      Sum = -(x * Sum) + exp (x);
   return Sum;
}

#endif

/*=========================================================================*/

static void snpair_AllocPoints (snpair_Res *res, long n)
{
   long i;
   WorkType *work = res->work;

   if (n <= 0)
      return;
   /* Allocates Maxnp tables of pointers to the points; one for each level
      of recursion */
   for (i = 1; i <= work->Maxnp; i++)
      res->Points[i] =
         util_Calloc ((size_t) (n + 1), sizeof (snpair_PointType));

   /* Allocates memory for the points; initially, only the first table
      of pointers, i = 1, points to the points. */
   for (i = 0; i <= n; i++)
      res->Points[1][i] = util_Calloc ((size_t) (work->kk + 1),
                                       sizeof (double));

   res->CloseDist = util_Calloc ((size_t) work->mcd + 1, sizeof (double));
}


/*=========================================================================*/

static void snpair_DeletePoints (snpair_Res * res)
/*
 * To clean up after the test.
 */
{
   long i;
   long n = res->n;
   WorkType *work = res->work;

   if (n <= 0)
      return;
   res->CloseDist = util_Free (res->CloseDist);

   for (i = 0; i <= n; i++)
      util_Free (res->Points[1][i]);

   for (i = 1; i <= work->Maxnp; i++)
      res->Points[i] = util_Free (res->Points[i]);
}


/*=========================================================================*/

static void AllocClosePairs (
   snpair_Res *res,           /* Results holder */
   long N,
   long n,
   int m
   )
{
   snpair_AllocPoints (res, n);
   res->Yn = statcoll_Create (m, "Yn: The m jumps of Y");
   res->Y = statcoll_Create (N * m + 100,
      "Y: All the jumps of Y, superposed");
   res->U = statcoll_Create (N * m,
      "U: The jumps of Y transformed into uniforms");
   res->V = statcoll_Create (N * m + 100, "V: A copy of the uniforms");
   res->S = statcoll_Create (N * m + 100, "S: Spacings");
   res->TheWn = statcoll_Create (N, "The N values of the W_n");
   res->TheWni = statcoll_Create (N * m, "The Nm values of the W_{n,i}");
   res->ThepValAD = statcoll_Create (N, "The p-values of A2");
   res->BitMax = statcoll_Create (N, "Largest bit distances");
}


/*=========================================================================*/

static void CleanClosePairs (snpair_Res * res)
{
   int i;
   res->Yn = statcoll_Delete (res->Yn);
   res->Y = statcoll_Delete (res->Y);
   res->U = statcoll_Delete (res->U);
   res->V = statcoll_Delete (res->V);
   res->S = statcoll_Delete (res->S);
   res->TheWn = statcoll_Delete (res->TheWn);
   res->TheWni = statcoll_Delete (res->TheWni);
   res->ThepValAD = statcoll_Delete (res->ThepValAD);
   res->BitMax = statcoll_Delete (res->BitMax);
   snpair_DeletePoints (res);
   for (i = 0; i < snpair_StatType_N; i++) {
      res->sVal[i] = -1.0;
      res->pVal[i] = -1.0;
   }
}


/*=========================================================================*/

static void InitRes (
   snpair_Res *res,           /* Results holder */
   long N,                    /* Number of replications */
   long n,                    /* Number of points */
   int m                      /* Number of closest distances kept */
   )
/* 
 * Initializes the res structure
 */
{
   if (res->CleanFlag)
      CleanClosePairs (res);
   AllocClosePairs (res, N, n, m);
   res->n = n;
   res->CleanFlag = TRUE;
}


/*-------------------------------------------------------------------------*/

snpair_Res *snpair_CreateRes (void)
{
   snpair_Res *res;
   res = util_Malloc (sizeof (snpair_Res));
   memset (res, 0, sizeof (snpair_Res));
   res->work = util_Malloc (sizeof (WorkType));
   res->CleanFlag = FALSE;
   return res;
}


/*-------------------------------------------------------------------------*/

void snpair_DeleteRes (snpair_Res * res)
{
   if (res == NULL)
      return;
   if (res->CleanFlag)
      CleanClosePairs (res);
   res->work = util_Free (res->work);
   util_Free (res);
}


/*=========================================================================*/

static void WriteSeuils (WorkType * work, lebool flag, double mu2,
   double nLR, double kLR)
{
   printf ("\n   Seuil1 = %2d\n   Seuil2 = %2d\n   "
      "Seuil3 = %2d\n   Seuil4 = %2d\n"
      "   L1 = %2d\n   L2 = %2d\n   s1 = ", snpair_env.Seuil1,
      snpair_env.Seuil2, snpair_env.Seuil3, snpair_env.Seuil4, work->L1,
      work->L2);

   /* s1 = n / 2^{kL1} */
   num_WriteD (nLR * pow (2.0, -kLR * work->L1), 9, 2, 2);
   printf ("\n   s2 = ");
   /* s2 = n / 2^{kL2} */
   num_WriteD (nLR * pow (2.0, -kLR * work->L2), 9, 2, 2);
   printf ("\n\n");

   if (flag) {
      printf ("   The minimal distance, to the power k, should be"
         " approximately\n      exponential with mean mu2 = ");
      num_WriteD (mu2, 12, 4, 2);
      printf ("\n\n   dlim1  = ");
      num_WriteD (work->dlim1, 15, 5, 3);
      printf ("\n   dlim1p = ");
      num_WriteD (work->dlim1p, 15, 5, 3);
      printf ("\n\n");
   }
}


/*=========================================================================*/

static void CalcSeuils (WorkType * work, long k, long m, lebool flag,
   double mu2, double nLR, double kLR)
{
   work->L1 = 1 + num_Log2 (nLR / snpair_env.Seuil3) / k;
   work->L2 = 1 + num_Log2 (nLR / snpair_env.Seuil4) / k;
   if (work->L1 < 1)
      work->L1 = 1;
   if (work->L2 < 1)
      work->L2 = 1;
   if (k < 6 && work->L1 < 2)
      work->L1 = 2;
   if (k < 6 && work->L2 < 2)
      work->L2 = 2;
   work->dlim1 = pow (m * mu2, 1.0 / k);
   work->dlim1p = pow (work->dlim1, work->pLR);
   if (swrite_Parameters)
      WriteSeuils (work, flag, mu2, nLR, kLR);
}


/*=========================================================================*/

void snpair_WriteDataCP (unif01_Gen * gen, char *TestName,
   long N, long n, int r, int t, int p, int m, lebool Torus)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",  t = %1d,", t);
   if (p >= 0)
      printf ("  p = %1d,", p);
   printf ("  m = %1d,  Torus = ", m);
   util_WriteBool (Torus, 5);
   printf ("\n\n");
}


/*=========================================================================*/
#define SPACINGS_9
/* 
 * La constante SPACINGS permet l'inclusion des 4 tests de ClosePairs suivants
 * qui sont mis en commentaire sinon: NP-S, NP-PR, mNP1-S, mNP2-S.
 */

void snpair_WriteResultsCP (unif01_Gen * gen, chrono_Chrono * Timer,
   snpair_Res * res, long N, long m)
{
   printf ("\n---------------------------------------\n");
   printf ("Test based on the 2 nearest points (NP):\n\n");

   if (N == 1) {
      printf ("The closest distance                  : ");
      num_WriteD (res->CloseDist[1], 7, 2, 2);
      printf ("\n");
      gofw_Writep1 (res->pVal[snpair_NP]);
   } else {
      printf ("Stat. AD on the N values (NP)         :");
      gofw_Writep2 (res->sVal[snpair_NP], res->pVal[snpair_NP]);
#ifdef SPACINGS
      printf ("Stat. AD after spacings (NP-S)        :");
      gofw_Writep2 (res->sVal[snpair_NPS], res->pVal[snpair_NPS]);
      printf ("Stat. AD after power ratio (NP-PR)    :");
      gofw_Writep2 (res->sVal[snpair_NPPR], res->pVal[snpair_NPPR]);
#endif
   }

   if (m > 1) {
      printf ("\nA2 test based on the spacings between the\n"
              "   successive jump times of process Y_n(t):\n\n");
      printf ("A2 test on the values of A2 (m-NP)    :");
      gofw_Writep2 (res->sVal[snpair_mNP], res->pVal[snpair_mNP]);

      if (N > 1) {
         printf ("Test on the Nm values of W_{n,i}(mNP1):");
         gofw_Writep2 (res->sVal[snpair_mNP1], res->pVal[snpair_mNP1]);
#ifdef SPACINGS
         printf ("Stat. AD after spacings (mNP1-S)      :");
         gofw_Writep2 (res->sVal[snpair_mNP1S], res->pVal[snpair_mNP1S]);
#endif
         printf ("Test on the jump times of Y\n   (superposition of Yn):\n\n");
         printf ("Expected number of jumps of Y = mN    : %7ld\n", m*N);
         printf ("Number of jumps of Y                  ");
         if (res->sVal[snpair_NJumps] >= N*NUM_JUMPS_LIM)
            printf ("> %6.0f     *****\n", res->sVal[snpair_NJumps]);
         else
            printf (": %7.0f\n", res->sVal[snpair_NJumps]);
         gofw_Writep1 (res->pVal[snpair_NJumps]);

         if (res->Y->NObs > 0) {
            printf ("Stat. AD (mNP2)                       :");
            gofw_Writep2 (res->sVal[snpair_mNP2], res->pVal[snpair_mNP2]);
#if 1
            if (snpair_mNP2S_Flag) {
               printf ("Stat. AD after spacings (mNP2-S)      :");
               gofw_Writep2 (res->sVal[snpair_mNP2S], res->pVal[snpair_mNP2S]);
            }
#endif
         }
      }
   }
   swrite_Final (gen, Timer);
}


/*=========================================================================*/

void snpair_ClosePairs (unif01_Gen * gen, snpair_Res * res,
   long N, long n, int r, int k, int p, int m)
/*
 * Looks at the m closest pairs in the torus and the first m jumps of the
 * process Y_n(t). A simplified version of snpair_ClosePairs.  
 */
{
   long j;
   long i;
   long Seq;
   double Wn;
   double t1;
   snpair_PointType T;
   double x;
   double NextJump;
   double Jump;
   double mu2;                    /* Expected minimum distance */
   double A2;
   double Vol;                    /* Volume of unit sphere in k dimension */
   double mLR, nLR, kLR;
   fmass_INFO Mass;
   double pLeft, pRight;
   statcoll_Collector *Q;
   WorkType *work;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "snpair_ClosePairs test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      snpair_WriteDataCP (gen, TestName, N, n, r, k, p, m, TRUE);

   /* util_Assert (k <= snpair_MaxDim, "snpair_ClosePairs: k >
      snpair_MaxDim");
   util_Assert (n <= snpair_MaxNumPoints,
      "snpair_ClosePairs:   n > snpair_MaxNumPoints"); */
   util_Assert (m > 0, "snpair_ClosePairs:   m <= 0");
   util_Assert (m <= snpair_MAXM, "snpair_ClosePairs:   m > snpair_MAXM");
   if (res == NULL) {
      localRes = TRUE;
      res = snpair_CreateRes ();
   }
   work = res->work;
   work->Torus = TRUE;
   work->kk = k;
   work->pp = p;
   work->mm = m;
   kLR = k;
   nLR = n;
   mLR = m;
   work->mcd = 2 * m;
   if (p == 0)
      work->pLR = 1.0;
   else
      work->pLR = p;
   work->Invp = 1.0 / work->pLR;
   if (k < snpair_MAXREC)
      work->Maxnp = k;
   else
      work->Maxnp = snpair_MAXREC;
   work->BBFlag = FALSE;           /* Bickel-Breiman Flag */

   Vol = num2_VolumeSphere ((double) p, k);
   mu2 = 2.0 / (nLR * (nLR - 1.0) * Vol);
   t1 = mLR;

   CalcSeuils (work, k, m, TRUE, mu2, nLR, kLR);
   InitRes (res, N, n, m);
   res->Distance = snpair_DistanceCP;
   res->VerifPairs = snpair_VerifPairs1;
   res->MiniProc = snpair_MiniProc1;

   /* Beginning of test */
   for (Seq = 1; Seq <= N; Seq++) {

      for (i = 1; i <= n; i++) {
         /* Generate n points in dimension k */
         T = res->Points[1][i];
         for (j = 1; j <= k; j++)
            T[j] = unif01_StripD (gen, r);
      }
      res->NumClose = 0;
      work->dlimp = work->dlim = kLR; /* Initial upper bounds */
      snpair_QuickSort (res->Points[1], 1, n, 1);
      snpair_FindClosePairs (res, 1, n, 1, 1, 1);
      Wn = 1.0 - exp (-pow (res->CloseDist[1], kLR) / mu2);
      statcoll_AddObs (res->TheWn, Wn);
      statcoll_Init (res->Yn, m);
      statcoll_Init (res->U, m);
      if (m > 1) {
         /* Calculate the spacings Delta_{n,i} between the jumps of Y_n, */
         /* then the W^*_{n(i)}, which are in principle i.i.d. U(0,1).  */
         Jump = 0.0;
         for (i = 1; i <= m; i++) {
            NextJump = pow (res->CloseDist[i], kLR) / mu2;
            statcoll_AddObs (res->Yn, NextJump);
            x = 1.0 - exp (-(NextJump - Jump));
            statcoll_AddObs (res->U, x);
            statcoll_AddObs (res->TheWni, x);
            Jump = NextJump;
         }

         tables_QuickSortD (res->U->V, 1, m);
         /* res->U should now contain m random var. i.i.d U(0,1), sorted */
         if (swrite_Collectors) {
            statcoll_Write (res->Yn, 5, 14, 4, 3);
            statcoll_Write (res->U, 5, 14, 4, 3);
         }
         A2 = gofs_AndersonDarling (res->U->V, m);
         x = fbar_AndersonDarling (m, A2);
         statcoll_AddObs (res->ThepValAD, x);

         if (N > 1) {
            /* Put in res->Y all the jumps between 0 and t1 */
            for (i = 1; i <= res->NumClose; i++) {
               NextJump = pow (res->CloseDist[i], kLR) / mu2;
               if (NextJump <= t1)
                  statcoll_AddObs (res->Y, NextJump);
            }
         }
      }
   }

   if (N == 1)
      res->pVal[snpair_NP] = 1.0 - Wn;
   else {
      Q = res->TheWn;
      tables_QuickSortD (Q->V, 1, N);
      /* Test NP at level 1 */
      res->sVal[snpair_NP] = gofs_AndersonDarling (Q->V, N);
      res->pVal[snpair_NP] = fbar_AndersonDarling (N, res->sVal[snpair_NP]);

#ifdef SPACINGS
      /* Test NP with the spacings */
      tables_CopyTabD (Q->V, res->V->V, 1, N);
      gofs_DiffD (res->V->V, res->S->V, 1, N, 0.0, 1.0);
      gofs_IterateSpacings (res->V->V, res->S->V, N);
      tables_QuickSortD (res->V->V, 1, N);
      res->sVal[snpair_NPS] = gofs_AndersonDarling (res->V->V, N);
      res->pVal[snpair_NPS] =
         fbar_AndersonDarling (N, res->sVal[snpair_NPS]);
      /* Test NP with power ratio */
      tables_CopyTabD (Q->V, res->V->V, 1, N);
      gofs_PowerRatios (res->V->V, N);
      tables_QuickSortD (res->V->V, 1, N);
      res->sVal[snpair_NPPR] = gofs_AndersonDarling (res->V->V, N);
      res->pVal[snpair_NPPR] =
         fbar_AndersonDarling (N, res->sVal[snpair_NPPR]);
#endif
   }

   if (m > 1) {
      if (N == 1) {
         res->sVal[snpair_mNP] = A2;
         res->pVal[snpair_mNP] = res->ThepValAD->V[1];
      } else {
         tables_CopyTabD (res->ThepValAD->V, res->V->V, 1, N);
         tables_QuickSortD (res->V->V, 1, N);
         res->sVal[snpair_mNP] = gofs_AndersonDarling (res->V->V, N);
         res->pVal[snpair_mNP] =
            fbar_AndersonDarling (N, res->sVal[snpair_mNP]);
         Q = res->TheWni;
         tables_QuickSortD (Q->V, 1, Q->NObs);
         /* Here, Q->NObs = N*m */
         res->sVal[snpair_mNP1] = gofs_AndersonDarling (Q->V, Q->NObs);
         res->pVal[snpair_mNP1] =
            fbar_AndersonDarling (Q->NObs, res->sVal[snpair_mNP1]);
#ifdef SPACINGS
         /* Test NP with the spacings */
         tables_CopyTabD (Q->V, res->V->V, 1, Q->NObs);
         gofs_DiffD (res->V->V, res->S->V, 1, Q->NObs, 0.0, 1.0);
         gofs_IterateSpacings (res->V->V, res->S->V, Q->NObs);
         tables_QuickSortD (res->V->V, 1, Q->NObs);
         res->sVal[snpair_mNP1S] = gofs_AndersonDarling (res->V->V, Q->NObs);
         res->pVal[snpair_mNP1S] =
            fbar_AndersonDarling (Q->NObs, res->sVal[snpair_mNP1S]);
#endif
         /* Superposition process of all the jumps of Y_n in [0, t1].
            Conditionnally on TotJumps = res->Y^.NObs, these jumps should be 
            uniformly distributed in [0, t1]. */

         Q = res->Y;
         for (i = 1; i <= Q->NObs; i++)
            Q->V[i] /= t1;
         if (Q->NObs > 0) {
            tables_QuickSortD (Q->V, 1, Q->NObs);
            /* res->Y must now contain random var. i.i.d U(0,1), sorted */
            res->sVal[snpair_mNP2] = gofs_AndersonDarling (Q->V, Q->NObs);
            res->pVal[snpair_mNP2] =
               fbar_AndersonDarling (Q->NObs, res->sVal[snpair_mNP2]);
	 }

         Mass = fmass_CreatePoisson (N * m);
         pLeft = fdist_Poisson2 (Mass, res->Y->NObs);
         pRight = fbar_Poisson2 (Mass, res->Y->NObs);
         fmass_DeletePoisson (Mass);
         res->sVal[snpair_NJumps] = res->Y->NObs;
         res->pVal[snpair_NJumps] = gofw_pDisc (pLeft, pRight);

#if 1
         /* Test on res->Y with the spacings */
         statcoll_Init (res->V, Q->Dim);
         statcoll_Init (res->S, Q->Dim);
         tables_CopyTabD (Q->V, res->V->V, 1, Q->NObs);
         gofs_DiffD (res->V->V, res->S->V, 1, Q->NObs, 0.0, 1.0);
         gofs_IterateSpacings (res->V->V, res->S->V, Q->NObs);
         tables_QuickSortD (res->V->V, 1, Q->NObs);
         res->sVal[snpair_mNP2S] = gofs_AndersonDarling (res->V->V, Q->NObs);
         res->pVal[snpair_mNP2S] =
            fbar_AndersonDarling (Q->NObs, res->sVal[snpair_mNP2S]);
#endif
      }
   }

   if (swrite_Collectors) {
      if (N > 1)
         statcoll_Write (res->Y, 5, 14, 4, 3);
      statcoll_Write (res->TheWn, 5, 14, 4, 3);
      statcoll_Write (res->TheWni, 5, 14, 4, 3);
      statcoll_Write (res->ThepValAD, 5, 14, 4, 3);
   }

   if (swrite_Basic)
      snpair_WriteResultsCP (gen, Timer, res, N, m);
   if (localRes)
      snpair_DeleteRes (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/
#if 0

void snpair_ReTestY (long N, long n, int m, double tt0, double tt1)
/*
 * Make more (experimental) tests on Y after a call to ClosePairs1
 * (To do?? We should renormalize the Y[i] by multiplying them by nM)

 *****************************************
  This procedure is not used anymore. It makes use of discontinuous
  distribution functions and the associated statistics are very
  complicated and did not seem sensitive.
  ****************************************
 */
{
   long i;
   double Bidon;
   double x;
   double iLR;
   double Fact;
   double HM;
   double HP;
   double GM;
   double GP;
   statcoll_Collector *Q = res->Y;
   fdist_FUNC_JUMPS *GPJumps;     /* All info on the jumps of G+ */
   fdist_FUNC_JUMPS *GMJumps;     /* All info on the jumps of G- */
   fdist_FUNC_JUMPS *HPJumps;     /* All info on the jumps of H+ */
   fdist_FUNC_JUMPS *HMJumps;     /* All info on the jumps of H- */

   /* Calculate statistics G+, G-, H+, H- */
   util_Assert (!swrite_AutoClean,
                "snpair_ReTestY:   swrite_AutoClean must be FALSE");
   util_Assert (res->Y != NULL,
                "snpair_ReTestY:   res->Y is a NULL pointer");

   Fact = N * m;
   for (i = 1; i <= Q->NObs; i++)
      Q->V[i] *= Fact;
   HP = 0.0;
   HM = 0.0;
   i = 1;
   while (Q->V[i] <= tt0 && i <= Q->NObs)
      ++i;
   GP = (i - 1) / tt0;
   while (Q->V[i] <= tt1 && i <= Q->NObs)
      ++i;
   GM = (i - 1) / tt1;
   if (i - 1 - tt1 < 0.0)
      HM = i - 1 - tt1;

   for (i = 1; i <= Q->NObs; i++) {
      if (Q->V[i] <= tt1) {
         x = Q->V[i];
         iLR = i;
         if (iLR - 1.0 - x < HM)
            HM = iLR - 1.0 - x;
         if (iLR - x > HP)
            HP = iLR - x;
         if (Q->V[i] >= tt0) {
            if (iLR / x > GP)
               GP = iLR / x;
            if ((iLR - 1.0) / x < GM)
               GM = (iLR - 1.0) / x;
         }
      }
   }

   res->pVal[snpair_GPlus] = 1.0 - FDistGPlus (Bidon, GP);
   res->sVal[snpair_GPlus] = GP;
   res->pVal[snpair_GMinus] = 1.0 - FDistGMinus (Bidon, GM);
   res->sVal[snpair_GMinus] = GM;
   res->pVal[snpair_HPlus] = 1.0 - FDistHPlus (Bidon, HP);
   res->sVal[snpair_HPlus] = HP;
   res->pVal[snpair_HMinus] = 1.0 - FDistHMinus (Bidon, HM);
   res->sVal[snpair_HMinus] = HM;

   if (swrite_Basic) {
      printf ("Test on the statistic G+              :");
      gofw_Writep2 (GP, res->pVal[snpair_GPlus]);
      printf ("\nTest on the statistic G-              :");
      gofw_Writep2 (GM, res->pVal[snpair_GMinus]);
      printf ("\nTest on the statistic H+              :");
      gofw_Writep2 (HP, res->pVal[snpair_HPlus]);
      printf ("\nTest on the statistic H-              :");
      gofw_Writep2 (HM, res->pVal[snpair_HMinus]);
      printf ("\n");
   }
}

#endif

/*=========================================================================*/

static void InitBBp0k2 (void)
/* 
 * Initialize Bickel-Breiman distribution with p = 0, k = 2
 */
{
   BB2[0] = 0.0;
   BB2[1] = 6.6022859e-5;
   BB2[2] = 2.111e-3;
   BB2[3] = 1.10679e-2;
   BB2[4] = 2.99898e-2;
   BB2[5] = 5.80398e-2;
   BB2[6] = 9.31672e-2;
   BB2[7] = 1.326804e-1;
   BB2[8] = 1.743017e-1;
   BB2[9] = 2.168632e-1;
   BB2[10] = 2.589057e-1;
   BB2[11] = 2.996407e-1;
   BB2[12] = 3.387514e-1;
   BB2[13] = 3.758668e-1;
   BB2[14] = 4.108985e-1;
   BB2[15] = 4.442291e-1;
   BB2[16] = 4.757295e-1;
   BB2[17] = 5.053408e-1;
   BB2[18] = 5.330166e-1;
   BB2[19] = 5.589979e-1;
   BB2[20] = 0.58358;
   BB2[21] = 6.067753e-1;
   BB2[22] = 6.281726e-1;
   BB2[23] = 6.483016e-1;
   BB2[24] = 6.670896e-1;
   BB2[25] = 6.848204e-1;
   BB2[26] = 7.016251e-1;
   BB2[27] = 7.17358e-1;
   BB2[28] = 7.319895e-1;
   BB2[29] = 7.458925e-1;
   BB2[30] = 7.589198e-1;
   BB2[31] = 7.712947e-1;
   BB2[32] = 7.82992e-1;
   BB2[33] = 7.939033e-1;
   BB2[34] = 8.044324e-1;
   BB2[35] = 8.14079e-1;
   BB2[36] = 8.233257e-1;
   BB2[37] = 8.319796e-1;
   BB2[38] = 8.402721e-1;
   BB2[39] = 0.84794;
   BB2[40] = 8.55173e-1;
   BB2[41] = 8.621625e-1;
   BB2[42] = 8.686291e-1;
   BB2[43] = 8.748893e-1;
   BB2[44] = 8.809175e-1;
   BB2[45] = 8.864606e-1;
   BB2[46] = 8.918151e-1;
   BB2[47] = 8.968588e-1;
   BB2[48] = 9.016258e-1;
   BB2[49] = 9.061704e-1;
   BB2[50] = 9.104084e-1;
   BB2[51] = 9.143753e-1;
   BB2[52] = 9.182571e-1;
   BB2[53] = 9.219001e-1;
   BB2[54] = 9.254505e-1;
   BB2[55] = 9.28772e-1;
   BB2[56] = 9.320133e-1;
   BB2[57] = 9.351203e-1;
   BB2[58] = 9.380936e-1;
   BB2[59] = 9.408087e-1;
   BB2[60] = 9.434704e-1;
   BB2[61] = 9.459541e-1;
   BB2[62] = 9.483018e-1;
   BB2[63] = 0.95057;
   BB2[64] = 9.526812e-1;
   BB2[65] = 9.547303e-1;
   BB2[66] = 9.566821e-1;
   BB2[67] = 9.585939e-1;
   BB2[68] = 9.603836e-1;
   BB2[69] = 9.62073e-1;
   BB2[70] = 9.637419e-1;
   BB2[71] = 9.652689e-1;
   BB2[72] = 9.666498e-1;
   BB2[73] = 9.680299e-1;
   BB2[74] = 9.693984e-1;
   BB2[75] = 9.707229e-1;
   BB2[76] = 9.720219e-1;
   BB2[77] = 9.731801e-1;
   BB2[78] = 9.742979e-1;
   BB2[79] = 9.753166e-1;
   BB2[80] = 9.763355e-1;
   BB2[81] = 9.773626e-1;
   BB2[82] = 9.782835e-1;
   BB2[83] = 9.792146e-1;
   BB2[84] = 9.800791e-1;
   BB2[85] = 9.808841e-1;
   BB2[86] = 9.816958e-1;
   BB2[87] = 9.824477e-1;
   BB2[88] = 9.831492e-1;
   BB2[89] = 9.838786e-1;
   BB2[90] = 9.845241e-1;
   BB2[91] = 9.851295e-1;
   BB2[92] = 9.857702e-1;
   BB2[93] = 9.863849e-1;
   BB2[94] = 9.869562e-1;
   BB2[95] = 9.875006e-1;
   BB2[96] = 9.879895e-1;
   BB2[97] = 9.88473e-1;
   BB2[98] = 9.889793e-1;
   BB2[99] = 9.894184e-1;
   BB2[100] = 9.898547e-1;
   BB2[101] = 9.902526e-1;
   BB2[102] = 9.906462e-1;
   BB2[103] = 9.910496e-1;
   BB2[104] = 9.914303e-1;
   BB2[105] = 9.918174e-1;
   BB2[106] = 9.921392e-1;
   BB2[107] = 9.924491e-1;
   BB2[108] = 9.92784e-1;
   BB2[109] = 9.930638e-1;
   BB2[110] = 9.933363e-1;
   BB2[111] = 9.936298e-1;
   BB2[112] = 9.93886e-1;
   BB2[113] = 9.941112e-1;
   BB2[114] = 9.943411e-1;
   BB2[115] = 9.945669e-1;
   BB2[116] = 9.947672e-1;
   BB2[117] = 9.94994e-1;
   BB2[118] = 9.951802e-1;
   BB2[119] = 9.953648e-1;
   BB2[120] = 9.955457e-1;
   BB2[121] = 9.957099e-1;
   BB2[122] = 9.959196e-1;
   BB2[123] = 9.961046e-1;
   BB2[124] = 9.962811e-1;
   BB2[125] = 9.964261e-1;
   BB2[126] = 9.965653e-1;
   BB2[127] = 9.967088e-1;
   BB2[128] = 0.99684;
   BB2[129] = 9.969537e-1;
   BB2[130] = 9.970835e-1;
   BB2[131] = 9.972087e-1;
}


/*-------------------------------------------------------------------------*/

static void InitBBp2k2 (void)
/* 
 * Initialize Bickel-Breiman distribution with p = 2, k = 2
 */
{
   BB4[0] = -8.2912955e-1;
   BB4[1] = -9.4432194e-1;
   BB4[2] = -1.0567132;
   BB4[3] = -1.1679847;
   BB4[4] = -1.2776563;
   BB4[5] = -1.384483;
   BB4[6] = -1.4916059;
   BB4[7] = -1.5956447;
   BB4[8] = -1.6994536;
   BB4[9] = -1.8012517;
   BB4[10] = -1.9014279;
   BB4[11] = -2.0006153;
   BB4[12] = -2.0997178;
   BB4[13] = -2.1987994;
   BB4[14] = -2.2959638;
   BB4[15] = -2.391997;
   BB4[16] = -2.4867876;
   BB4[17] = -2.5815698;
   BB4[18] = -2.6761017;
   BB4[19] = -2.7658218;
   BB4[20] = -2.8582757;
   BB4[21] = -2.9522569;
   BB4[22] = -3.0406141;
   BB4[23] = -3.1311066;
   BB4[24] = -3.2179075;
   BB4[25] = -3.3057192;
   BB4[26] = -3.3933087;
   BB4[27] = -3.4815725;
   BB4[28] = -3.5719191;
   BB4[29] = -3.6592077;
   BB4[30] = -3.7437809;
   BB4[31] = -3.8274559;
   BB4[32] = -3.9149689;
   BB4[33] = -4.000307;
   BB4[34] = -4.0874655;
   BB4[35] = -4.1724253;
   BB4[36] = -4.2619679;
   BB4[37] = -4.3498336;
   BB4[38] = -4.4349335;
   BB4[39] = -4.5214761;
   BB4[40] = -4.607099;
   BB4[41] = -4.6921565;
   BB4[42] = -4.7799781;

   BB5[0] = -4.5909e-3;
   BB5[1] = -3.666e-4;
   BB5[2] = 7.508e-5;
   BB5[3] = 2.15483e-3;
   BB5[4] = 1.115755e-2;
   BB5[5] = 3.033271e-2;
   BB5[6] = 5.881422e-2;
   BB5[7] = 9.422896e-2;
   BB5[8] = 1.3423286e-1;
   BB5[9] = 1.7618124e-1;
   BB5[10] = 2.1865118e-1;
   BB5[11] = 2.6082507e-1;
   BB5[12] = 3.0215075e-1;
   BB5[13] = 3.4140313e-1;
   BB5[14] = 3.7898955e-1;
   BB5[15] = 4.1454877e-1;
   BB5[16] = 4.4830003e-1;
   BB5[17] = 4.7980029e-1;
   BB5[18] = 5.093375e-1;
   BB5[19] = 5.3717465e-1;
   BB5[20] = 5.6357091e-1;
   BB5[21] = 5.8817876e-1;
}


/*-------------------------------------------------------------------------*/

static void InitBBp0k15 (void)
/* 
 * Initialize Bickel-Breiman distribution with p = 0, k = 15
 */
{
   BB3[0] = 0.0;
   BB3[1] = 1.6778e-4;
   BB3[2] = 2.6967455e-3;
   BB3[3] = 1.28187e-2;
   BB3[4] = 3.25519e-2;
   BB3[5] = 0.06001;
   BB3[6] = 9.28778e-2;
   BB3[7] = 1.292254e-1;
   BB3[8] = 1.674211e-1;
   BB3[9] = 2.066797e-1;
   BB3[10] = 2.439418e-1;
   BB3[11] = 2.805974e-1;
   BB3[12] = 3.156376e-1;
   BB3[13] = 3.487236e-1;
   BB3[14] = 3.804003e-1;
   BB3[15] = 4.103833e-1;
   BB3[16] = 4.394161e-1;
   BB3[17] = 4.673735e-1;
   BB3[18] = 4.935018e-1;
   BB3[19] = 5.181638e-1;
   BB3[20] = 5.403617e-1;
   BB3[21] = 5.609553e-1;
   BB3[22] = 5.813387e-1;
   BB3[23] = 6.003938e-1;
   BB3[24] = 6.188892e-1;
   BB3[25] = 6.353537e-1;
   BB3[26] = 6.509678e-1;
   BB3[27] = 6.658608e-1;
   BB3[28] = 6.797704e-1;
   BB3[29] = 6.93198e-1;
   BB3[30] = 7.059891e-1;
   BB3[31] = 7.185495e-1;
   BB3[32] = 7.306549e-1;
   BB3[33] = 7.413516e-1;
   BB3[34] = 7.517439e-1;
   BB3[35] = 7.617069e-1;
   BB3[36] = 7.709877e-1;
   BB3[37] = 7.804121e-1;
   BB3[38] = 7.898275e-1;
   BB3[39] = 7.984594e-1;
   BB3[40] = 8.058125e-1;
   BB3[41] = 8.129582e-1;
   BB3[42] = 8.199078e-1;
   BB3[43] = 8.26626e-1;
   BB3[44] = 8.332602e-1;
   BB3[45] = 8.393936e-1;
   BB3[46] = 8.452292e-1;
   BB3[47] = 8.510694e-1;
   BB3[48] = 8.569731e-1;
   BB3[49] = 8.621826e-1;
   BB3[50] = 8.671328e-1;
   BB3[51] = 8.723293e-1;
   BB3[52] = 8.770461e-1;
   BB3[53] = 8.814338e-1;
   BB3[54] = 8.853624e-1;
   BB3[55] = 8.897322e-1;
   BB3[56] = 8.937578e-1;
   BB3[57] = 8.97254e-1;
   BB3[58] = 9.008031e-1;
   BB3[59] = 9.042233e-1;
   BB3[60] = 9.076829e-1;
   BB3[61] = 9.11218e-1;
   BB3[62] = 9.139078e-1;
   BB3[63] = 9.170002e-1;
   BB3[64] = 9.199191e-1;
   BB3[65] = 9.226127e-1;
   BB3[66] = 9.250731e-1;
   BB3[67] = 9.277341e-1;
   BB3[68] = 9.301693e-1;
   BB3[69] = 9.324761e-1;
   BB3[70] = 9.347405e-1;
   BB3[71] = 9.370394e-1;
   BB3[72] = 9.390614e-1;
   BB3[73] = 9.41111e-1;
   BB3[74] = 9.429319e-1;
   BB3[75] = 9.448513e-1;
   BB3[76] = 9.466235e-1;
   BB3[77] = 9.483763e-1;
   BB3[78] = 9.500882e-1;
   BB3[79] = 9.517579e-1;
   BB3[80] = 9.531616e-1;
   BB3[81] = 9.546471e-1;
   BB3[82] = 9.561263e-1;
   BB3[83] = 9.576014e-1;
   BB3[84] = 9.592471e-1;
   BB3[85] = 9.605977e-1;
   BB3[86] = 9.618122e-1;
   BB3[87] = 9.632723e-1;
   BB3[88] = 9.644877e-1;
   BB3[89] = 9.654043e-1;
   BB3[90] = 9.666469e-1;
   BB3[91] = 9.676583e-1;
   BB3[92] = 9.687529e-1;
   BB3[93] = 9.697718e-1;
   BB3[94] = 9.708359e-1;
   BB3[95] = 9.716986e-1;
   BB3[96] = 9.726066e-1;
   BB3[97] = 9.734057e-1;
   BB3[98] = 9.743224e-1;
   BB3[99] = 9.751716e-1;
   BB3[100] = 9.759489e-1;
   BB3[101] = 9.766958e-1;
   BB3[102] = 9.774256e-1;
   BB3[103] = 9.783317e-1;
   BB3[104] = 9.789422e-1;
   BB3[105] = 9.795293e-1;
   BB3[106] = 9.801187e-1;
   BB3[107] = 9.807522e-1;
   BB3[108] = 9.812972e-1;
   BB3[109] = 9.818664e-1;
   BB3[110] = 9.825167e-1;
   BB3[111] = 9.831091e-1;
   BB3[112] = 9.835873e-1;
   BB3[113] = 9.840919e-1;
   BB3[114] = 9.845122e-1;
   BB3[115] = 9.850374e-1;
   BB3[116] = 9.854874e-1;
   BB3[117] = 9.859857e-1;
   BB3[118] = 9.865129e-1;
   BB3[119] = 9.869294e-1;
   BB3[120] = 9.873618e-1;
   BB3[121] = 9.877482e-1;
   BB3[122] = 9.880475e-1;
}


/*-------------------------------------------------------------------------*/

static double FDistBBp0k2 (double junk[], double x)
/*
 * Bickel-Breiman distribution obtained by simulation with 
 *    N = 1000000,  n = 1000,  r = 0,  k =  2,  p = 0,  Torus =  TRUE
 * 
 * We first interpolated the empirical distribution on the points xs = j/100
 * (integer j) by building a parabola using a least-square fit with all
 * the points in [xs - 0.005, xs + 0.005], and then by computing ys(xs) on
 * the parabola, in order to reduce the noise. 
 * We use a Newton cubic interpolation with the 4 points closest to x to 
 * compute the distribution y(x).
 */
{
   static lebool BBp0k2Flag = FALSE;
   int j;
   double q;
   double y;
   if (x >= 6.0)
      return 1.0;
   if (x >= 1.3)
      return 1.0 - exp (-5.94558e-1 - 3.99672 * x);
   if (x <= 0.014)
      return 0.0;
   if (x <= 0.02) {
      return -2.66337e-3 + x * (5.12234e-1 + x * (-32.8023 + 701.167 * x));
   }

   if (FALSE == BBp0k2Flag) {
      InitBBp0k2 ();
      BBp0k2Flag = TRUE;
   }

   j = 100.0 * x + 2;             /* x is in [P(j-2), P(j-1)] */
   q = 100.0 * x - j;

   /* Newton backward cubic interpolation */
   y = BB2[j - 1] + (BB2[j - 1] - BB2[j - 2]) * q + (((BB2[j - 3]
            - 2.0 * BB2[j - 2]) + BB2[j - 1]) * q * (q + 1.0)) / 2.0
      + ((((-BB2[j - 4] + 3.0 * BB2[j - 3]) -
            3.0 * BB2[j - 2]) + BB2[j - 1]) * q * (q + 1.0) * (q +
         2.0)) / 6.0;

   return y;
}


/*-------------------------------------------------------------------------*/

static double FDistBBp0k15 (double junk[], double x)
/*
 * Bickel-Breiman distribution obtained by simulation with 
 *    N = 100000,  n = 1000,  r = 0,  k =  15,  p = 0,  Torus =  TRUE
 * 
 * We first interpolated the empirical distribution on the points xs = j/100
 * (integer j) by building a parabola using a least-square fit with all
 * the points in [xs - 0.005, xs + 0.005], and then by computing ys(xs) on
 * the parabola, in order to reduce the noise. 
 * We use a Newton cubic interpolation with the 4 points closest to x to 
 * compute the distribution y(x).
 */
{
   static lebool BBp0k15Flag = FALSE;
   int j;
   double q;
   double y;
   if (x <= 0.015)
      return 0.0;
   if (x <= 0.02)
      return (6.1123 * x - 0.18384) * x + 1.3984e-3;
   if (x >= 6.0)
      return 1.0;
   if (x >= 1.2)
      return 1.0 - exp (-3.15786 * x - 5.41639e-1);

   if (FALSE == BBp0k15Flag) {
      InitBBp0k15 ();
      BBp0k15Flag = TRUE;
   }

   j = 100.0 * x + 2;             /* x is in [P(j-2), P(j-1)] */
   q = 100.0 * x - j;

   /* Newton backward cubic interpolation */
   y = (BB3[j - 1] - BB3[j - 2]) * q + BB3[j - 1] +
      (BB3[j - 3] - 2.0 * BB3[j - 2] + BB3[j - 1]) * q * (q + 1.0) / 2.0
      + (-BB3[j - 4] + 3.0 * BB3[j - 3] - 3.0 * BB3[j - 2]
      + BB3[j - 1]) * q * (q + 1.0) * (q + 2.0) / 6.0;

   return y;
}


/*-------------------------------------------------------------------------*/

static double FDistBBp2k2 (double junk[], double x)
/*
 * Bickel-Breiman distribution obtained by simulation with 
 *    N = 1000000,  n = 1000,  r = 0,  k =  2,  p = 2,  Torus =  TRUE
 */
{
   static lebool BBp2k2Flag = FALSE;
   int j;
   double q;
   double y;
   if (x < 0.016)
      return 0.0;
   if (x >= 6.0)
      return 1.0;
   if (x >= 1.0)
      return 1.0 - exp ((0.1408724 * x - 4.485674) * x - 0.264116);

   if (FALSE == BBp2k2Flag) {
      InitBBp2k2 ();
      BBp2k2Flag = TRUE;
   }

   if (x >= 0.2) {
      /* Newton quadratic interpolation based on the points 0.02*j */
      /* in the interval [0.2, 1.0] */
      j = x * 50.0;
      q = x * 50.0 - j;
      y = BB4[j - 10] + q * (BB4[j - 9] - BB4[j - 10]) +
         (q * (q - 1.0) * ((BB4[j - 8] - 2.0 * BB4[j - 9]) + BB4[j -
               10])) / 2.0;

      return 1.0 - exp (y);
   }

   /* Newton backward cubic interpolation based on the points */
   /* 0.01*j in the interval [0, 0.2] */
   j = 100.0 * x + 2;             /* x is in [P(j-2), P(j-1)] */
   q = 100.0 * x - j;

   y = (BB5[j] - BB5[j - 1]) * q + BB5[j]
      + (((BB5[j - 2] - 2.0 * BB5[j - 1]) + BB5[j])
      * q * (q + 1.0)) / 2.0
      + ((((-BB5[j - 3] + 3.0 * BB5[j - 2]) - 3.0 * BB5[j - 1])
         + BB5[j]) * q * (q + 1.0) * (q + 2.0)) / 6.0;

   return y;
}


/*-------------------------------------------------------------------------*/

void snpair_WriteDataBB (unif01_Gen * gen, char *TestName,
   long N, long n, int r, int k, int p, lebool Torus, int L1, int L2)
{
   double z;

   swrite_Head (gen, TestName, N, n, r);
   printf (",  k = %1d,  p = %1d,   Torus = ", k, p);
   util_WriteBool (Torus, 5);
   printf ("\n");

   if (swrite_Parameters) {
      printf ("\n   Seuil1 = %5d\n   Seuil2 = %5d\n   Seuil3 = %5d\n"
         "   Seuil4 = %5d\n   L1 = %2d\n   L2 = %2d\n",
         snpair_env.Seuil1, snpair_env.Seuil2, snpair_env.Seuil3,
         snpair_env.Seuil4, L1, L2);

      z = n * pow (2.0, -L1 * (double) k);
      printf ("   s1 = ");        /* n / 2^{k L1} = "); */
      num_WriteD (z, 9, 2, 2);
      printf ("\n   s2 = ");      /* n / 2^{k L2} = "); */
      z = n * pow (2.0, -L2 * (double) k);
      num_WriteD (z, 9, 2, 2);
   }
   printf ("\n\n\n");
}


/*-------------------------------------------------------------------------*/

void snpair_WriteResultsBB (unif01_Gen * gen, chrono_Chrono * Timer,
   snpair_Res * res, long N)
{

   printf ("-----------------------------------------------\n");
   if (N == 1) {
      printf ("Value of the BB statistic             :");
      gofw_Writep2 (res->sVal[snpair_BB], res->pVal[snpair_BB]);
   } else {
      printf ("AD Statistic on the N p-values of BB  :");
      gofw_Writep2 (res->sVal[snpair_BB], res->pVal[snpair_BB]);
   }
   swrite_Final (gen, Timer);
}


/*=========================================================================*/

void snpair_BickelBreiman (unif01_Gen * gen, snpair_Res * res,
   long N, long n, int r, int k, int p, lebool Torus)
{
   int j;
   long i;
   long Seq;
   snpair_PointType T;
   double mu1;                    /* -n * Vol */
   double ksurp;                  /* k / p */
   double Wni;
   double SumBB;
   double Vol;                    /* Volume of unit sphere in k dimension */
   double x, nLR, kLR;
   WorkType *work;
   lebool localRes = FALSE;
   chrono_Chrono *Timer, *Time1;
   char *TestName = "snpair_BickelBreiman test";

   Timer = chrono_Create ();

   if (res == NULL) {
      localRes = TRUE;
      res = snpair_CreateRes ();
   }
   work = res->work;
   work->Torus = Torus;
   work->kk = k;
   kLR = k;
   nLR = n;
   work->pp = p;
   work->mm = 1;
   work->mcd = 2;
   if (p == 0)
      work->pLR = 1.0;
   else
      work->pLR = p;
   work->Invp = 1.0 / work->pLR;
   ksurp = kLR / work->pLR;

   work->L1 = 1 + num_Log2 (nLR / snpair_env.Seuil3) / (sqrt (kLR));
   if (work->L1 < 2)
      work->L1 = 2;
   work->L2 = 1 + num_Log2 (nLR / snpair_env.Seuil4) / (sqrt (kLR));
   if (work->L2 < 2)
      work->L2 = 2;
   if (k < snpair_MAXREC)
      work->Maxnp = k;
   else
      work->Maxnp = snpair_MAXREC;
   Vol = num2_VolumeSphere ((double) p, k);
   mu1 = -nLR * Vol;
   work->BBFlag = TRUE;
   if (swrite_Basic)
      snpair_WriteDataBB (gen, TestName, N, n, r, k, p, Torus,
         work->L1, work->L2);

   /*  util_Assert (n <= snpair_MaxNumPoints,
       "snpair_BickelBreiman:   n is too large"); */
   util_Assert (p == 2 || p == 0,
      "snpair_BickelBreiman implemented only for p = 2 and p = 0");
   util_Assert (k == 2 || k == 15,
      "snpair_BickelBreiman implemented only for k = 2 and k = 15");
   util_Assert (p != 2 || k != 15,
      "snpair_BickelBreiman:   case p = 2, k = 15  not implemented");
   if (p == 0) {
      if (k == 2)
         work->FDistBB = FDistBBp0k2;
      else
         work->FDistBB = FDistBBp0k15;
   } else
      work->FDistBB = FDistBBp2k2;

   InitRes (res, N, n, 1);
   res->Distance = snpair_DistanceBB;
   res->VerifPairs = snpair_VerifPairs0;
   res->MiniProc = snpair_MiniProc1;
   statcoll_SetDesc (res->ThepValAD, "The N p-values of BickelBreiman");

   /* Test begins */
   for (Seq = 1; Seq <= N; Seq++) {

      for (i = 1; i <= n; i++) {
         /* Generate n points in dimension k */
         T = res->Points[1][i];
         /* Initialize nearest distance */
         T[0] = kLR;
         for (j = 1; j <= k; j++)
            T[j] = unif01_StripD (gen, r);
      }

      /* Find the closest points */
      work->dlim = kLR;            /* Initial upper bounds */
      work->dlimp = work->dlim;
      if (snpair_TimeBB)
         Time1 = chrono_Create ();
      snpair_QuickSort (res->Points[1], 1, n, 1);
      snpair_FindClosePairs (res, 1, n, 1, 1, 1);

      /* For each point, coordinate 0 now contains the distance to */
      /* the nearest point raised to power p (for p > 0) */
      snpair_QuickSort (res->Points[1], 1, n, 0);

      /* Compute the BB statistic, etc...  */
      SumBB = 0.0;
      for (i = 1; i <= n; i++) {
         Wni = 1.0 - exp (mu1 * pow (res->Points[1][i][0], ksurp));
         x = Wni - i / nLR;
         SumBB += x * x;
      }
      if (snpair_TimeBB) {
         printf ("   Time to compute the BB statistic:  ");
         chrono_Write (Time1, chrono_sec);
         printf ("\n");
         chrono_Delete (Time1);
      }
      statcoll_AddObs (res->ThepValAD,
         1.0 - work->FDistBB ((double *) NULL, SumBB));
   }

   if (swrite_Collectors)
      statcoll_Write (res->ThepValAD, 5, 14, 4, 3);

   if (N == 1) {
      res->sVal[snpair_BB] = SumBB;
      res->pVal[snpair_BB] = res->ThepValAD->V[1];
   } else {
      tables_QuickSortD (res->ThepValAD->V, 1, N);
      res->sVal[snpair_BB] = gofs_AndersonDarling (res->ThepValAD->V, N);
      res->pVal[snpair_BB] = fbar_AndersonDarling (N, res->sVal[snpair_BB]);
   }

   if (swrite_Basic)
      snpair_WriteResultsBB (gen, Timer, res, N);

   if (localRes)
      snpair_DeleteRes (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void snpair_DistanceCPBitM (snpair_Res * res, snpair_PointType P1,
   snpair_PointType P2)
/*
 * Similar to snpair_DistanceCP, but for snpair_ClosePairsBitMatch. We take
 * at most two groups of SizeUL bits for each coordinate of two points and
 * find how many equal bits (= Y) they have before the first different bit,
 * starting with the most significant. That is the distance in 1 dimension.
 * We do that for each coordinate and the minimum of these is the distance
 * between the two points (all components of the pair have at least Y
 * identical bits).
 */
{
   const int NBitsUL = CHAR_BIT * sizeof (unsigned long);
   const double Mul = num_TwoExp[NBitsUL];
   unsigned long x1, x2, z;
   int i, j;
   int Y = INT_MAX;               /* Distance between the 2 points */
   WorkType *work = res->work;

   for (i = 1; i <= work->kk; i++) {
      /* Take the first NBitsUL bits of each coordinates of the 2 points */
      x1 = Mul * P1[i];
      x2 = Mul * P2[i];
      /* Find the position - 1 of the first (left) bit where they differ */
      z = x1 ^ x2;
      j = 0;
      if (z) {
         while (z < 2 * z) {
            j++;
            z <<= 1;
            if (j >= Y)
               continue;
         }
      } else {
         /* The first NBitsUL bits are equal, consider the NBitsUL next bits 
          */
         x1 = Mul * (Mul * P1[i] - x1);
         x2 = Mul * (Mul * P2[i] - x2);
         z = x1 ^ x2;
         if (z) {
            j = NBitsUL;
            while (z < 2 * z) {
               j++;
               z <<= 1;
               if (j >= Y)
                  continue;
            }
         } else {
            j = 2 * NBitsUL;
         }
      }
      if (j < Y)
         Y = j;
      if (Y <= work->YLim)
         /* We want the maximum (amongst all pairs of points) of the */
         /* minimum Y over all coordinates of a pair. This pair cannot */
         /* give a larger YLim. */
         return;
   }

   /* A larger YLim has been found. From it, we define an inverse distance */
   /* so that the largest YLim gives the smallest new distance. This is */
   /* necessary if we want to use the fast but complicated algorithm for */
   /* finding the nearest pair. */
   if (Y > work->YLim) {
      work->YLim = Y;
      if (work->YLim <= num_MaxTwoExp)
         work->dlim = 1.0 / num_TwoExp[work->YLim];
      else
         work->dlim = pow (2.0, -(double) work->YLim);
      res->CloseDist[1] = work->dlim;
   }
}


/*-------------------------------------------------------------------------*/

static void WriteDataBM (unif01_Gen * gen, char *TestName,
   long N, long n, int r, int k)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",  t = %1d\n\n", k);
}


/*=========================================================================*/

void snpair_ClosePairsBitMatch (unif01_Gen * gen, snpair_Res * res,
   long N, long n, int r, int k)
/*
 * Similar to ClosePairs, but uses the BitMatch distance.
 */
{
   long Seq;
   double z1, nLR;
   snpair_PointType T;
   int m;
   int MaxY;                      /* Max of all bit distances */
   int j;
   long i;
   double pLeft, pRight;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "snpair_ClosePairsBitMatch test";
   WorkType *work;

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataBM (gen, TestName, N, n, r, k);

   /*  util_Assert (n <= snpair_MaxNumPoints,
       "snpair_ClosePairsBitMatch:   n > snpair_MaxNumPoints"); */
   util_Assert (n > 1, "snpair_ClosePairsBitMatch:   n < 2");
   if (res == NULL) {
      localRes = TRUE;
      res = snpair_CreateRes ();
   }
   work = res->work;
   work->Torus = FALSE;
   work->kk = k;
   work->mm = m = 1;
   work->mcd = 2 * m;
   nLR = n;
   work->Invp = work->pLR = work->pp = 1;

   if (k < snpair_MAXREC)
      work->Maxnp = k;
   else
      work->Maxnp = snpair_MAXREC;
   work->BBFlag = FALSE;           /* Bickel-Breiman Flag */

   CalcSeuils (work, k, m, FALSE, 0.0, nLR, (double) k);

   InitRes (res, N, n, m);
   res->Distance = snpair_DistanceCPBitM;
   res->VerifPairs = snpair_VerifPairs1;
   res->MiniProc = snpair_MiniProc1;

   MaxY = 0;

   /* Beginning of test */
   for (Seq = 1; Seq <= N; Seq++) {

      for (i = 1; i <= n; i++) {
         /* Generate n points in dimension k */
         T = res->Points[1][i];
         for (j = 1; j <= k; j++)
            T[j] = unif01_StripD (gen, r);
      }
      res->NumClose = 0;
      work->YLim = 0;              /* Initial lower bound */
      work->dlim = 1.0;            /* Initial upper bound */
      snpair_QuickSort (res->Points[1], 1, n, 1);
      snpair_FindClosePairs (res, 1, n, 1, 1, 1);

#if 0
      /* Check by computing distances between all pairs; very slow. */
      printf ("%12d", work->YLim);
      if (Seq % 5 == 0)
         printf ("\n");
      swrite_Collectors = TRUE;
      work->YLim = 0;
      snpair_VerifPairs0 (res->Points[1], 1, n, 0, 0);
#endif

      statcoll_AddObs (res->BitMax, (double) work->YLim);
      MaxY = util_Max (work->YLim, MaxY);
   }

   if (swrite_Collectors)
      statcoll_Write (res->BitMax, 5, 14, 4, 3);

   /* z1 = Probability [Min {k geometric (0.5)} >= MaxY] */
   if (k * (MaxY + 1) <= num_MaxTwoExp)
      z1 = 1.0 / num_TwoExp[k * (MaxY + 1)];
   else
      z1 = pow (2.0, -(double) k * (MaxY + 1));

   /* There are n*(n - 1)/2 pairs of points and we replicate that basic test 
      N times, so we compute pLeft = the Probability [Max {N*n*(n - 1)/2 of
      above random var.} <= MaxY] */
   if (z1 > DBL_EPSILON) {
      pLeft = 1.0 - z1;
      z1 = log (pLeft) * N * n * (n - 1) / 2;
      pLeft = exp (z1);
      pRight = 1.0 - pLeft;
   } else {
      /* Use approximation log (1 - z) = -z to avoid loss of precision */
      pRight = z1 * N * n * (n - 1) / 2;
      pLeft = 1.0 - pRight;
   }
   res->pVal[snpair_BM] = gofw_pDisc (pLeft, pRight);
   res->sVal[snpair_BM] = MaxY;

   if (swrite_Basic) {
      printf ("\n-----------------------------------------------\n");
      printf ("Max of all bit distances              :");
      gofw_Writep2 ((double) MaxY, res->pVal[snpair_BM]);
      swrite_Final (gen, Timer);
   }

   if (localRes)
      snpair_DeleteRes (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/
