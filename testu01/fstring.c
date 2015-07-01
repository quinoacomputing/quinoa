/*************************************************************************\
 *
 * Package:        TestU01
 * File:           fstring.c
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
#include "gofs.h"

#include "fstring.h"
#include "fcho.h"
#include "ffam.h"
#include "fres.h"
#include "ftab.h"
#include "sstring.h"
#include "unif01.h"

#include <string.h>


long fstring_Maxn = 1024*1024*4;
long fstring_MaxL = 1024*1024*1;


/* Type of structure to choose the parameter of the sample size */
typedef enum {
   CHO_CHO,
   CHO_STRING
} cho_Type;



/*------------------------------ Functions --------------------------------*/

static void InitRes1 (
   ffam_Fam *fam,
   fstring_Res1 *res,         /* Results holder */
   int N,                     /* Number of replications */
   int Nr,
   int j1, int j2, int jstep,
   char *name1,
   char *name2
)
/* 
 * Initializes the fstring_Res1 structure
 */
{
   fres_InitCont (fam, res->BLen, N, Nr, j1, j2, jstep, name1);
   fres_InitDisc (fam, res->GLen, Nr, j1, j2, jstep, name2);
}


/*-------------------------------------------------------------------------*/

fstring_Res1 * fstring_CreateRes1 (void)
{
   fstring_Res1 *res;
   res = util_Malloc (sizeof (fstring_Res1));
   res->BLen = fres_CreateCont ();
   res->GLen = fres_CreateDisc ();
   return res;
}


/*-------------------------------------------------------------------------*/

void fstring_DeleteRes1 (fstring_Res1 *res)
{
   if (res == NULL)
      return;
   fres_DeleteCont (res->BLen);
   fres_DeleteDisc (res->GLen);
   util_Free (res);
}


/*=========================================================================*/

static void InitRes2 (
   ffam_Fam *fam,
   fstring_Res2 *res,         /* Results holder */
   int N,                     /* Number of replications */
   int Nr,
   int j1, int j2, int jstep,
   char *name1,
   char *name2
)
/* 
 * Initializes the fstring_Res2 structure
 */
{
   fres_InitCont (fam, res->NBits, N, Nr, j1, j2, jstep, name1);
   fres_InitCont (fam, res->NRuns, N, Nr, j1, j2, jstep, name2);
}


/*-------------------------------------------------------------------------*/

fstring_Res2 * fstring_CreateRes2 (void)
{
   fstring_Res2 *res;
   res = util_Malloc (sizeof (fstring_Res2));
   res->NBits = fres_CreateCont ();
   res->NRuns = fres_CreateCont ();
   return res;
}


/*-------------------------------------------------------------------------*/

void fstring_DeleteRes2 (fstring_Res2 *res)
{
   if (res == NULL)
      return;
   fres_DeleteCont (res->NBits);
   fres_DeleteCont (res->NRuns);
   util_Free (res);
}


/*=========================================================================*/

static void PrintHead (char *test, ffam_Fam *fam,
   long N, long n, int r, int s, long L, int d, 
   int Nr, int j1, int j2, int jstep)
{
   printf
   ("\n\n================================================================\n");
   printf ("Family:  %s\n\n", fam->name);
   printf ("Test:    %s\n", test);
   printf ("   N  = %ld,", N);
   if (n)
      printf ("   n = %ld,", n);
   printf ("   r = %d,   s = %d", r, s);
   if (L)
      printf (",   L = %ld", L);
   if (d)
      printf (",   d = %d", d);
   printf ("\n   Nr = %d,   j1 = %d,   j2 = %d,   jstep = %d\n\n",
      Nr, j1, j2, jstep);
}


/*=========================================================================*/

static int ChooseParam (int prec, void *cho, cho_Type ver,
   long *pn, int *pr, int *ps, long *pL, long LMin, int i, int j)
/*
 * Set the values of the parameters for the test. If a parameter is < 0,
 * will call a choose function to set it. Otherwise, will accept it as is.
 * Returns 0 if parameters are ok for the test, returns -1 if the test
 * should not be done for these parameters.
 */
{
   fcho_Cho *cho1 = cho;
   fcho_Cho2 *cho2 = cho;
   fcho_Cho *chon;
   fcho_Cho *choL;

   switch (ver) {
   case CHO_STRING:
      util_Assert (cho2, "fstring:   cho2 is NULL");
      chon = cho2->Chon;
      choL = cho2->Chop2;
      util_Assert (*pn < 0 || *pL < 0, 
         "fstring:   Either n or L must be < 0" );
      break;
   case CHO_CHO:
      chon = cho1;
      break;
   default:
      util_Error ("in fstring, ChooseParam:  no such case");
   }

   if (*pn < 0) {
      util_Assert (chon, "fstring:   n < 0 and chon is NULL");
      *pn = chon->Choose (chon->param, i, j);

      if (*pn <= 3.0 * gofs_MinExpected) {
	 printf ("n is too small\n\n");
	 return -1;
      }
      if (*pn > fstring_Maxn) {
	 printf ("n > %2ld\n\n", fstring_Maxn);
	 return -1;
      }
   }

   *ps = fcho_Chooses (*pr, *ps, prec);
   if (*ps <= 0)
      return -1;

   if (*pL < 0) {
      util_Assert (choL, "fstring:   L < 0 and choL is NULL");
      *pL = choL->Choose (choL->param, i, j);

      if (*pL <= LMin) {
	 printf ("L is too small\n\n");
	 return -1;
      }
      if (*pL > fstring_MaxL) {
	 printf ("L > %2ld\n\n", fstring_MaxL);
	 return -1;
      }
   }

   return 0;
}


/*=========================================================================*/

static void TabPeriod (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, s;
   long N, n, junk = 0;
   const long *Par = par1;
   fres_Cont *fres = res1;
   sres_Chi2 *sres;

   N = Par[0];
   n = Par[1];
   r = Par[2];
   s = Par[3];

   if (ChooseParam (fam->Resol[irow], cho, CHO_CHO, &n, &r, &s, &junk, 0, i, j))
      return;

   sres = sres_CreateChi2 ();
   sstring_PeriodsInStrings (fam->Gen[irow], sres, N, n, r, s);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteChi2 (sres);
}


/*-------------------------------------------------------------------------*/

void fstring_Period1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
   long N, int r, int s,
   int Nr, int j1, int j2, int jstep)
{
   long Par[5] = { 0 };
   lebool localRes;

   Par[0] = N;
   Par[1] = -1;
   Par[2] = r;
   Par[3] = s;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fstring_Period1", fam, N, 0, r, s, 0, 0, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fstring_Period1");
   ftab_MakeTables (fam, res, cho, Par, TabPeriod, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*=========================================================================*/

static void TabRun (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, s;
   long N, n, junk = 0;
   const long *Par = par1;
   fstring_Res2 *fres = res1;
   sstring_Res3 *sres;

   N = Par[0];
   n = Par[1];
   r = Par[2];
   s = Par[3];

   if (ChooseParam (fam->Resol[irow], cho, CHO_CHO, &n, &r, &s, &junk, 0, i, j))
      return;

   sres = sstring_CreateRes3 ();
   sstring_Run (fam->Gen[irow], sres, N, n, r, s);
   fres_FillTableEntryC (fres->NRuns, sres->NRuns->pVal2, N, irow, icol);
   fres_FillTableEntryC (fres->NBits, sres->NBits->pVal2, N, irow, icol);
   sstring_DeleteRes3 (sres);
}


/*-------------------------------------------------------------------------*/

void fstring_Run1 (ffam_Fam *fam, fstring_Res2 *res, fcho_Cho *cho,
   long N, int r, int s,
   int Nr, int j1, int j2, int jstep)
{
   long Par[5] = {0};
   lebool localRes;

   Par[0] = N;
   Par[1] = -1;
   Par[2] = r;
   Par[3] = s;
   if (res == NULL) {
      localRes = TRUE;
      res = fstring_CreateRes2 ();
   } else
      localRes = FALSE;

   PrintHead ("fstring_Run1", fam, N, 0, r, s, 0, 0, Nr, j1, j2, jstep);
   InitRes2 (fam, res, N, Nr, j1, j2, jstep, "fstring_Run1, Number of Bits",
       "fstring_Run1, Number of Runs");
   ftab_MakeTables (fam, res, cho, Par, TabRun, Nr, j1, j2, jstep);
   fres_PrintCont (res->NRuns);
   fres_PrintCont (res->NBits);
   if (localRes)
      fstring_DeleteRes2 (res);
}


/*=========================================================================*/

static void TabAutoCor (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, s, d;
   long N, n, junk = 0;
   const long *Par = par1;
   fres_Cont *fres = res1;
   sres_Basic *sres;

   N = Par[0];
   n = Par[1];
   r = Par[2];
   s = Par[3];
   d = Par[5];

   if (ChooseParam (fam->Resol[irow], cho, CHO_CHO, &n, &r, &s, &junk, 0, i, j))
      return;

   sres = sres_CreateBasic  ();
   sstring_AutoCor (fam->Gen[irow], sres, N, n, r, s, d);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteBasic (sres);
}


/*-------------------------------------------------------------------------*/

void fstring_AutoCor1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
   long N, int r, int s, int d,
   int Nr, int j1, int j2, int jstep)
{
   long Par[6] = {0};
   lebool localRes;

   Par[0] = N;
   Par[1] = -1;
   Par[2] = r;
   Par[3] = s;
   Par[5] = d;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fstring_AutoCor1", fam, N, 0, r, s, 0, d, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fstring_AutoCor1");
   ftab_MakeTables (fam, res, cho, Par, TabAutoCor, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*=========================================================================*/

static void TabLongHead (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, s;
   long N, n, L;
   const long *Par = par1;
   fstring_Res1 *fres = res1;
   sstring_Res2 *sres;

   N = Par[0];
   n = Par[1];
   r = Par[2];
   s = Par[3];
   L = Par[4];

   if (ChooseParam (fam->Resol[irow], cho, CHO_STRING,
          &n, &r, &s, &L, 1050, i, j))
      return;
   if (L < 1000 + s)
      return;

   sres = sstring_CreateRes2 ();
   sstring_LongestHeadRun (fam->Gen[irow], sres, N, n, r, s, L);
   fres_FillTableEntryC (fres->BLen, sres->Chi->pVal2, N, irow, icol);
   fres_FillTableEntryD (fres->GLen, sres->Disc->pLeft, sres->Disc->pRight,
                         sres->Disc->pVal2, irow, icol);
   sstring_DeleteRes2 (sres);
}


/*-------------------------------------------------------------------------*/

void fstring_LongHead1 (ffam_Fam *fam, fstring_Res1 *res, fcho_Cho2 *cho,
   long N, long n, int r, int s, long L,
   int Nr, int j1, int j2, int jstep)
{
   long Par[5] = {0};
   lebool localRes;

   Par[0] = N;
   Par[1] = n;
   Par[2] = r;
   Par[3] = s;
   Par[4] = L;
   if (res == NULL) {
      localRes = TRUE;
      res = fstring_CreateRes1 ();
   } else
      localRes = FALSE;

   PrintHead ("fstring_LongHead1", fam, N, n, r, s, L, 0, Nr, j1, j2, jstep);
   InitRes1 (fam, res, N, Nr, j1, j2, jstep,
      "fstring_LongHead1, n block lengths",
      "fstring_LongHead1, 1 global length");
   ftab_MakeTables (fam, res, cho, Par, TabLongHead, Nr, j1, j2, jstep);
   fres_PrintCont (res->BLen);
   ftab_PrintTable (res->GLen->PVal2);
   if (localRes)
      fstring_DeleteRes1 (res);
}


/*=========================================================================*/

static void TabHamWeight2 (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, s, ver;
   long N, n, L;
   const long *Par = par1;
   fres_Cont *fres = res1;
   sres_Chi2 *sres1;
   sres_Basic *sres2;

   N = Par[0];
   n = Par[1];
   r = Par[2];
   s = Par[3];
   L = Par[4];
   ver = Par[5];

   if (ChooseParam (fam->Resol[irow], cho, CHO_STRING, &n, &r, &s, &L, 0, i, j))
      return;
   if ((ver == 2) && (L > n)) {
      printf ("L > n\n\n");
      return;
   }
   if ((ver == 1) && (n <= 2.0 * gofs_MinExpected)) {
      printf ("n <= 2 gofs_MinExpected\n\n");
      return;
   }
   if (ver == 2) {
      sres2 = sres_CreateBasic  ();
      sstring_HammingWeight2 (fam->Gen[irow], sres2, N, n, r, s, L);
      fres_FillTableEntryC (fres, sres2->pVal2, N, irow, icol);
      sres_DeleteBasic (sres2);
   } else {
      sres1 = sres_CreateChi2  ();
      sstring_HammingWeight (fam->Gen[irow], sres1, N, n, r, s, L);
      fres_FillTableEntryC (fres, sres1->pVal2, N, irow, icol);
      sres_DeleteChi2 (sres1);
   }
}


/*-------------------------------------------------------------------------*/

static void Ver_HamWeight (ffam_Fam *fam, fres_Cont *res,
   fcho_Cho2 *cho, long N, long n, int r, int s, long L,
   int Nr, int j1, int j2, int jstep, int ver)
{
   long Par[6] = {0};
   lebool localRes;
   char Name[60];

   Par[0] = N;
   Par[1] = n;
   Par[2] = r;
   Par[3] = s;
   Par[4] = L;
   Par[5] = ver;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   if (ver == 2)
     strcpy (Name, "fstring_HamWeight2");
   else
     strcpy (Name, "fstring_HamWeight1");

   PrintHead (Name, fam, N, n, r, s, L, 0, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, Name);
   ftab_MakeTables (fam, res, cho, Par, TabHamWeight2, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*-------------------------------------------------------------------------*/

void fstring_HamWeight2 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
   long N, long n, int r, int s, long L,
   int Nr, int j1, int j2, int jstep)
{
   Ver_HamWeight (fam, res, cho, N, n, r, s, L, Nr, j1, j2, jstep, 2);
}


void fstring_HamWeight1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
   long N, long n, int r, int s, long L,
   int Nr, int j1, int j2, int jstep)
{
   Ver_HamWeight (fam, res, cho, N, n, r, s, L, Nr, j1, j2, jstep, 1);
}


/*=========================================================================*/

static void TabHamCorr (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, s;
   long N, n, L;
   const long *Par = par1;
   fres_Cont *fres = res1;
   sstring_Res *sres;

   N = Par[0];
   n = Par[1];
   r = Par[2];
   s = Par[3];
   L = Par[4];

   if (ChooseParam (fam->Resol[irow], cho, CHO_STRING, &n, &r, &s, &L, 0, i, j))
      return;

   sres = sstring_CreateRes  ();
   sstring_HammingCorr (fam->Gen[irow], sres, N, n, r, s, L);
   fres_FillTableEntryC (fres, sres->Bas->pVal2, N, irow, icol);
   sstring_DeleteRes (sres);
}


/*-------------------------------------------------------------------------*/

void fstring_HamCorr1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
   long N, long n, int r, int s, long L,
   int Nr, int j1, int j2, int jstep)
{
   long Par[5];
   lebool localRes;

   Par[0] = N;
   Par[1] = n;
   Par[2] = r;
   Par[3] = s;
   Par[4] = L;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fstring_HamCorr1", fam, N, n, r, s, L, 0, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fstring_HamCorr1");
   ftab_MakeTables (fam, res, cho, Par, TabHamCorr, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*=========================================================================*/

static void TabHamIndep (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, s;
   long N, n, L;
   const long *Par = par1;
   fres_Cont *fres = res1;
   sstring_Res *sres;

   N = Par[0];
   n = Par[1];
   r = Par[2];
   s = Par[3];
   L = Par[4];

   if (ChooseParam (fam->Resol[irow], cho, CHO_STRING, &n, &r, &s, &L, 0, i, j))
      return;

   sres = sstring_CreateRes  ();
   sstring_HammingIndep (fam->Gen[irow], sres, N, n, r, s, L, -1);
   fres_FillTableEntryC (fres, sres->Bas->pVal2, N, irow, icol);
   sstring_DeleteRes (sres);
}


/*-------------------------------------------------------------------------*/

void fstring_HamIndep1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
   long N, long n, int r, int s, long L,
   int Nr, int j1, int j2, int jstep)
{
   long Par[5];
   lebool localRes;

   Par[0] = N;
   Par[1] = n;
   Par[2] = r;
   Par[3] = s;
   Par[4] = L;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fstring_HamIndep1", fam, N, n, r, s, L, 0, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fstring_HamIndep1");
   ftab_MakeTables (fam, res, cho, Par, TabHamIndep, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*=========================================================================*/
