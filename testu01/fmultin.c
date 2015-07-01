/*************************************************************************\
 *
 * Package:        TestU01
 * File:           fmultin.c
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
#include "num.h"
#include "util.h"

#include "fmultin.h"
#include "ffam.h"
#include "fres.h"
#include "fcho.h"
#include "ftab.h"
#include "smultin.h"

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>


long fmultin_Maxn = 16777216;

#define LEN 127
#define LEN2 50
#define LEN3 10

typedef enum {
   A_MULT,
   A_BITS,
   A_PERM
} TestType;

typedef enum {
   B_DT,
   B_2HT,
   B_2L,
   B_T
} KType;

typedef struct {
   long N;                    /* Number of replications */
   int t;                     /* Dimension */
   double R;
   char *name;
   KType kcho;
} Multin_Param;

static const double EPSILON = 1.0E-14;

/* The stable parameters values used by default */
static smultin_Param ParamDefault = {
   2,                             /* NbDelta */
   {-1, 1},                       /* ValDelta */
   smultin_GenerCellSerial,       /* GenerCell */
   -1                             /* bmax */
};


/* The approximation used in the distribution of the CollisionOver test */
#define NMESS 5
static char COAMessages[NMESS][LEN3 + 1];




/*------------------------------- Functions ------------------------------*/


static void InitCOAStrings (void)
/*
 * Type of approximation used in CollisionOver test.
 */
{
   strncpy (COAMessages[smultin_CollNormal], "   N    ", (size_t) LEN3);
   strncpy (COAMessages[smultin_CollPoissonSparse], "   C    ", (size_t) LEN3);
   strncpy (COAMessages[smultin_CollPoissonDense], "   V    ", (size_t) LEN3);
   strncpy (COAMessages[smultin_CollNotInit], "  ---   ", (size_t) LEN3);
}


/*========================================================================*/


static void InitRes (
   ffam_Fam *fam,
   fmultin_Res *res,
   smultin_Param *par,
   int N,
   int Nr,
   int j1, int j2, int jstep,
   char *name,
   lebool Over
)
/* 
 * Initializes the fmultin_Res structure
 */
{
   int s, i;
   char str[LEN + 1] = { 0 };
   char str2[LEN2 + 1] = { 0 };
   smultin_Param *parold = res->Par;
   Nr = util_Min (Nr, fam->Ng);

   for (s = par->NbDelta; s < parold->NbDelta; s++)
      fres_DeleteCont (res->PowDiv[s]);
   for (s = parold->NbDelta; s < par->NbDelta; s++)
      res->PowDiv[s] = fres_CreateCont ();

   for (s = 0; s < par->NbDelta; s++) {
      if (fabs (par->ValDelta[s] + 1.0) < EPSILON) {
         strncpy (str, name, (size_t) LEN);
         if (Over) {
            strncat (str, ": CollisionOver test", (size_t) LEN2);
            ftab_DeleteTable (res->COApprox);
            res->COApprox = ftab_CreateTable (Nr, j1, j2, jstep,
               "Approximation used for distribution of CollisionOver",
               ftab_String, NMESS);
            for (i = 0; i < Nr; i++)
               res->COApprox->LSize[i] = fam->LSize[i];
            InitCOAStrings ();
            for (i = 0; i < NMESS; i++) {
               res->COApprox->Strings[i] = COAMessages[i];
	    }
         } else {
            strncat (str, ": Collision test", (size_t) LEN2);
         }
         fres_InitPoisson (fam, res->Coll, Nr, j1, j2, jstep, str);
         strncpy (str, name, (size_t) LEN);
         strncat (str, ": empty cells", (size_t) LEN2);
         fres_InitPoisson (fam, res->Empty, Nr, j1, j2, jstep, str);
         for (i = 1; i <= par->bmax; i++) {
            strncpy (str, name, (size_t) LEN);
            strncat (str, ": cells with at least ", (size_t) LEN2);
            sprintf (str2, "%1d", i);
            strncat (str, str2, (size_t) 3);
            strncat (str, " balls", (size_t) LEN2);
            fres_InitPoisson (fam, res->Balls[i], Nr, j1, j2, jstep, str);
         }
      }

      strncpy (str, name, (size_t) LEN);
      strncat (str, ": ValDelta = ", (size_t) LEN2);
      sprintf (str2, "%6.3f,", par->ValDelta[s]);
      strncat (str, str2, (size_t) LEN2);
      strncat (str, " test", (size_t) LEN2);
      fres_InitCont (fam, res->PowDiv[s], N, Nr, j1, j2, jstep, str);
   }
   /*
   ftab_DeleteTable (res->CellRatio);
   res->CellRatio = ftab_CreateTable (Nr, j1, j2, jstep,
      "Actual number of cells / Chosen number of cells", ftab_Real, 0);
   */
}


/*-------------------------------------------------------------------------*/

fmultin_Res *fmultin_CreateRes (smultin_Param * par)
{
   fmultin_Res *res;
   int s, j;

   res = util_Malloc (sizeof (fmultin_Res));
   if (par == NULL)
      par = &ParamDefault;
   res->Par = par;

   for (s = 0; s < par->NbDelta; s++)
      res->PowDiv[s] = fres_CreateCont ();

   /* For the collision test */
   res->Coll = fres_CreatePoisson ();
   res->Empty = fres_CreatePoisson ();
   for (j = 1; j <= par->bmax; j++)
      res->Balls[j] = fres_CreatePoisson ();
   res->COApprox = ftab_CreateTable (1, 1, 1, 1,
   "Approximation used for distribution of CollisionOver", ftab_String, 4);
   /*
   res->CellRatio = ftab_CreateTable (1, 1, 1, 1,
      "Actual number of cells / Chosen number of cells", ftab_Real, 0);
   */
   return res;
}


/*-------------------------------------------------------------------------*/

void fmultin_DeleteRes (fmultin_Res * res)
{
   int s, j;
   if (res == NULL)
      return;
   for (s = 0; s < res->Par->NbDelta; s++)
      fres_DeleteCont (res->PowDiv[s]);
   fres_DeletePoisson (res->Coll);
   fres_DeletePoisson (res->Empty);
   for (j = 1; j <= res->Par->bmax; j++)
      fres_DeletePoisson (res->Balls[j]);
   ftab_DeleteTable (res->COApprox);
   /* ftab_DeleteTable (res->CellRatio); */
   util_Free (res);
}


/*========================================================================*/

static void PrintHead (char *name, TestType test, ffam_Fam * fam,
   smultin_Param * spar, long *Par, int Nr, int j1, int j2, int jstep)
{
   int j;

   printf
      ("\n\n================================================================\n");
   printf ("Family:  %s\n\n", fam->name);
   printf ("Test:    %s\n", name);
   printf ("   N  = %ld,   r = %ld", Par[0], Par[1]);

   switch (test) {
   case A_BITS:
      printf (",   s = %ld,   Sparse = ", Par[2]);
      break;
   case A_MULT:
      printf (",   t = %ld,   Sparse = ", Par[3]);
      break;
   case A_PERM:
      printf (",   Sparse = ");
      break;
   default:
      util_Error ("in fmultin, PrintHead:  no such case");
   }
   util_WriteBool (Par[4], -5);
   printf ("\n   Nr = %d,   j1 = %d,   j2 = %d,   jstep = %d\n\n",
      Nr, j1, j2, jstep);

   if (spar) {
      printf ("   NbDelta = %d,   ValDelta = { ", spar->NbDelta);
      for (j = 0; j < spar->NbDelta; j++) {
	 printf ("%5.3g", spar->ValDelta[j]);
	 if (j < spar->NbDelta - 1)
	    printf (", ");
	 else
	    printf (" }\n\n ");
      }
   }
}


/*=========================================================================*/

static void PrintRes (fmultin_Res * res, lebool Over)
{
   int s, j;
   if (res == NULL)
      return;
   for (s = 0; s < res->Par->NbDelta; s++) {
      if (fabs (res->Par->ValDelta[s] + 1.0) > EPSILON)
         fres_PrintCont (res->PowDiv[s]);
      if (fabs (res->Par->ValDelta[s] + 1.0) < EPSILON) {
         fres_PrintPoisson (res->Coll, FALSE, FALSE);
         if (res->Par->bmax >= 0)
            fres_PrintPoisson (res->Empty, FALSE, TRUE);
         for (j = 2; j <= res->Par->bmax; j++)
            fres_PrintPoisson (res->Balls[j], FALSE, FALSE);
         /* ftab_PrintTable (res->CellRatio); */
         if (Over)
            ftab_PrintTable (res->COApprox);
      }
   }
}


/*=========================================================================*/

static fcho_Cho *CreateKcho (long N, int t, double R, KType kcho)
{
   fcho_Cho *cho;
   Multin_Param *Par;

   cho = util_Malloc (sizeof (fcho_Cho));
   Par = util_Malloc (sizeof (Multin_Param));
   Par->N = N;
   Par->t = t;
   Par->R = R;
   Par->kcho = kcho;
   cho->param = Par;
   cho->name = util_Calloc (2, sizeof (char));
   Par->name = cho->name;
   return cho;
}

/*-------------------------------------------------------------------------*/

static double CheckK1 (void *vpar, double K, long n)
{
   Multin_Param *Par = vpar;
   KType kcho = Par->kcho;
   double d, x;
   int h, L, t;

   if (n / K < 1.0 / num_TwoExp[30])
      return -1.0;
   if (K > smultin_env.Maxk) {
      printf ("K > smultin_env->Maxk\n\n");
      return -1.0;
   }
   switch (kcho) {
   case B_DT:
      d = pow (K, 1.0 / Par->t);
      strcpy (Par->name, "d");
      if (d > LONG_MAX)
         return -1.0;
      else
         return d;
      break;

   case B_2HT:
      h = 0.5 + num_Log2(K) / Par->t;
      d = num_TwoExp[h];
      strcpy (Par->name, "d");
      if (d > LONG_MAX)
         return -1.0;
      else
         return d;
      break;

   case B_2L:
      L = 0.5 + num_Log2(K);
      strcpy (Par->name, "L");
      return L;
      break;

   case B_T:
      strcpy (Par->name, "t");
      /* Find smallest t such that t! > k */
      x = 2.0;
      t = 2;
      while (x < K) {
         t++;
         x *= t;
      }
      /* Choose t or (t-1) that gives t! closest to K */
      if ((x - K)/K > (K - x/t)/K)
	 t--;
      return t;      
      break;

   default:
      util_Error ("in fmultin, CheckK1:  no such case");
   }
   return -1.0;
}


/*=========================================================================*/

static void WriteEC (void *vpar, long junk1, long junk2)
{
   Multin_Param *Par = vpar;
   switch (Par->kcho) {
   case B_DT:
      printf ("Choose  EC_DT with EC = ");
      break;
   case B_2HT:
      printf ("Choose  EC_2HT with EC = ");
      break;
   case B_2L:
      printf ("Choose  EC_2L with EC = ");
      break;
   case B_T:
      printf ("Choose  EC_T with EC = ");
      break;
   default:
      util_Error ("in fmultin, WriteEC:  no such case");
   }
   num_WriteD (Par->R, 8, 2, 2);
   printf ("\n\n");
}

/*-------------------------------------------------------------------------*/

static double ChooseEC (void *vpar, long junk, long n)
{
   Multin_Param *Par = vpar;
   double K;
   WriteEC (vpar, 0, 0);
   K = ((double) n * n * Par->N) / (2.0 * Par->R);
   return CheckK1 (vpar, K, n);
}


/*=========================================================================*/

void fmultin_DeleteEC (fcho_Cho * cho)
{
   if (NULL == cho)
      return;
   cho->name = util_Free (cho->name);
   cho->param = util_Free (cho->param);
   util_Free (cho);
}

/*-------------------------------------------------------------------------*/

fcho_Cho *fmultin_CreateEC_DT (long N, int t, double EC)
{
   fcho_Cho *cho;
   cho = CreateKcho (N, t, EC, B_DT);
   cho->Write = WriteEC;
   cho->Choose = ChooseEC;
   return cho;
}

/*-------------------------------------------------------------------------*/

fcho_Cho * fmultin_CreateEC_2HT (long N, int t, double EC)
{
   fcho_Cho *cho;
   cho = CreateKcho (N, t, EC, B_2HT);
   cho->Write = WriteEC;
   cho->Choose = ChooseEC;
   return cho;
}

/*-------------------------------------------------------------------------*/

fcho_Cho * fmultin_CreateEC_2L (long N, double EC)
{
   fcho_Cho *cho;
   cho = CreateKcho (N, -1, EC, B_2L);
   cho->Write = WriteEC;
   cho->Choose = ChooseEC;
   return cho;
}

/*-------------------------------------------------------------------------*/

fcho_Cho * fmultin_CreateEC_T (long N, double EC)
{
   fcho_Cho *cho;
   cho = CreateKcho (N, -1, EC, B_T);
   cho->Write = WriteEC;
   cho->Choose = ChooseEC;
   return cho;
}


/*=========================================================================*/

static void WriteDens (void *vpar, long junk1, long junk2)
{
   Multin_Param *Par = vpar;
   switch (Par->kcho) {
   case B_DT:
      printf ("Choose  Dens_DT with density = ");
      break;
   case B_2HT:
      printf ("Choose  Dens_2HT with density = ");
      break;
   case B_2L:
      printf ("Choose  Dens_2L with density = ");
      break;
   case B_T:
      printf ("Choose  Dens_T with density = ");
      break;
   default:
      util_Error ("in fmultin, WriteDens:  no such case");
   }
   if (Par->R > 0.999999)
      num_WriteD (Par->R, 8, 2, 2);
   else {
      printf (" 1 /");
      num_WriteD (1.0 / Par->R, 8, 2, 2);
   }
   printf ("\n\n");
}

/*-------------------------------------------------------------------------*/

static double ChooseDens (void *vpar, long junk, long n)
{
   Multin_Param *Par = vpar;
   double K;
   WriteDens (vpar, 0, 0);
   K = n / Par->R;
   return CheckK1 (vpar, K, n);
}


/*=========================================================================*/

void fmultin_DeleteDens (fcho_Cho * cho)
{
   if (NULL == cho)
      return;
   cho->name = util_Free (cho->name);
   cho->param = util_Free (cho->param);
   util_Free (cho);
}

/*-------------------------------------------------------------------------*/

fcho_Cho *fmultin_CreateDens_DT (int t, double R)
{
   fcho_Cho *cho;
   cho = CreateKcho (0, t, R, B_DT);
   cho->Write = WriteDens;
   cho->Choose = ChooseDens;
   return cho;
}

/*-------------------------------------------------------------------------*/

fcho_Cho * fmultin_CreateDens_2HT (int t, double R)
{
   fcho_Cho *cho;
   cho = CreateKcho (0, t, R, B_2HT);
   cho->Write = WriteDens;
   cho->Choose = ChooseDens;
   return cho;
}

/*-------------------------------------------------------------------------*/

fcho_Cho * fmultin_CreateDens_2L (double R)
{
   fcho_Cho *cho;
   cho = CreateKcho (0, -1, R, B_2L);
   cho->Write = WriteDens;
   cho->Choose = ChooseDens;
   return cho;
}

/*-------------------------------------------------------------------------*/

fcho_Cho * fmultin_CreateDens_T (double R)
{
   fcho_Cho *cho;
   cho = CreateKcho (0, -1, R, B_T);
   cho->Write = WriteDens;
   cho->Choose = ChooseDens;
   return cho;
}


/*=========================================================================*/

static void WritePer (void *vpar, long junk1, long junk2)
{
   Multin_Param *Par = vpar;
   switch (Par->kcho) {
   case B_DT:
      printf ("Choose  Per_DT with R = ");
      break;
   case B_2HT:
      printf ("Choose  Per_2HT with R = ");
      break;
   case B_2L:
      printf ("Choose  Per_2L with R = ");
      break;
   case B_T:
      printf ("Choose  Per_T with R = ");
      break;
   default:
      util_Error ("in fmultin, WritePer:  no such case");
   }
   if (Par->R > 0.999999)
      num_WriteD (Par->R, 8, 2, 2);
   else {
      printf (" 1 /");
      num_WriteD (1.0 / Par->R, 8, 2, 2);
   }
   printf ("\n\n");
}

/*-------------------------------------------------------------------------*/

static double ChoosePer (void *vpar, long lsize, long n)
{
   Multin_Param *Par = vpar;
   double K;
   WritePer (vpar, 0, 0);
   K = Par->R * pow (2.0, (double) lsize);
   return CheckK1 (vpar, K, n);
}


/*=========================================================================*/

void fmultin_DeletePer (fcho_Cho * cho)
{
   if (NULL == cho)
      return;
   cho->name = util_Free (cho->name);
   cho->param = util_Free (cho->param);
   util_Free (cho);
}

/*-------------------------------------------------------------------------*/

fcho_Cho *fmultin_CreatePer_DT (int t, double R)
{
   fcho_Cho *cho;
   cho = CreateKcho (0, t, R, B_DT);
   cho->Write = WritePer;
   cho->Choose = ChoosePer;
   return cho;
}

/*-------------------------------------------------------------------------*/

fcho_Cho * fmultin_CreatePer_2HT (int t, double R)
{
   fcho_Cho *cho;
   cho = CreateKcho (0, t, R, B_2HT);
   cho->Write = WritePer;
   cho->Choose = ChoosePer;
   return cho;
}

/*-------------------------------------------------------------------------*/

fcho_Cho * fmultin_CreatePer_2L (double R)
{
   fcho_Cho *cho;
   cho = CreateKcho (0, -1, R, B_2L);
   cho->Write = WritePer;
   cho->Choose = ChoosePer;
   return cho;
}

/*-------------------------------------------------------------------------*/

fcho_Cho * fmultin_CreatePer_T (double R)
{
   fcho_Cho *cho;
   cho = CreateKcho (0, -1, R, B_T);
   cho->Write = WritePer;
   cho->Choose = ChoosePer;
   return cho;
}


/*=========================================================================*/

static void FillTables (fmultin_Res *fres, smultin_Res *sres, long N,
    int irow, int icol, lebool Over)
{
  int i, j;
   for (i = 0; i < fres->Par->NbDelta; i++) {
      fres_FillTableEntryC (fres->PowDiv[i], sres->pVal2[i], N, irow, icol);
      if (fabs (fres->Par->ValDelta[i] + 1.0) < EPSILON) {
         fres_FillTableEntryPoisson (fres->Coll, sres->Mu[i],
            sres->NbCollisions, -1.0, -1.0, sres->pColl, irow, icol);
         fres_FillTableEntryPoisson (fres->Empty, sres->EsEmpty,
            (double) sres->NbCells[0], -1.0, -1.0, sres->pEmpty, irow, icol);
         for (j = 1; j <= fres->Par->bmax; j++)
            fres_FillTableEntryPoisson (fres->Balls[j], sres->EsCells[j],
            (double) sres->WbCells[j], -1.0, -1.0, sres->pWb[j], irow, icol);
         if (Over)
            fres->COApprox->Mat[irow][icol] = sres->CollApprox;
      }
   }
}


/*=========================================================================*/

static void TabMultin (ffam_Fam *fam, void *vres, void *vcho,
   void *vpar, int i, int j, int irow, int icol)
{
   int r, t;
   long N, n, d;
   lebool Sparse, Over;
   smultin_Res *sres;
   fmultin_Res *fres = vres;
   fcho_Cho2 *cho = vcho;
   fcho_Cho *chon;
   fcho_Cho *chop2;
   long *Par = vpar;
   Multin_Param *kpar;
   TestType Test;

   N = Par[0];
   r = Par[1];
   d = Par[2];
   t = Par[3];
   Sparse = Par[4];
   Over = Par[5];
   Test = Par[6];

   util_Assert (cho, "fmultin:   cho is NULL");
   chon = cho->Chon;
   chop2 = cho->Chop2;
   util_Assert (chon, "fmultin:   cho->Chon is NULL");
   util_Assert (chop2, "fmultin:   cho->Chop2 is NULL");
   kpar = chop2->param;
   if (Test == A_PERM) {
      util_Assert (kpar->kcho == B_T,
      "cho->Chop2:  wrong function for choosing number of cells");
   } else if (Test == A_MULT) {
      util_Assert ((kpar->kcho == B_DT) || (kpar->kcho == B_2HT),
      "cho->Chop2:  wrong function for choosing number of cells");
   }

   n = fcho_ChooseParamL (chon, 5, fmultin_Maxn, i, j);
   if (n < 0)
      return;
   if (d < 0) {
      strncpy (chop2->name, "d", 1);
      d = fcho_ChooseParamL (chop2, 2, LONG_MAX, i, n);
      if (d < 0)
         return;
      /* The number of bits used must be <= resolution of generator */
      if (r + num_Log2 ((double) d) + 0.5 > fam->Resol[irow]) {
	 printf ("Resolution of generator too small\n\n");
	 return;
      }
   } else {
      strncpy (chop2->name, "t", 1);
      t = fcho_ChooseParamL (chop2, 2, 18, i, n);
      if (t < 0)
         return;
   }
   if (Over && t < 2) {
      printf ("t < 2\n\n");
      return;
   }

   sres = smultin_CreateRes (fres->Par);
   if (Over)
      smultin_MultinomialOver (fam->Gen[irow], fres->Par, sres,
         N, n, r, d, t, Sparse);
   else
      smultin_Multinomial (fam->Gen[irow], fres->Par, sres,
         N, n, r, d, t, Sparse);

   FillTables (fres, sres, N, irow, icol, Over);
   smultin_DeleteRes (sres);
}


/*========================================================================*/

static int Chooses (int r, int s, int prec, int L)
{
   if (r + s > prec)
      s = prec - r;
   if (s <= 0) {
      printf ("r >= Resolution of generator\n\n");
      return -1;
   }
   if (L >= s) {
      while (L % s)
	 s--;
   } else {
      while (s % L)
	 s--;
   }
   return s;
}


/*-------------------------------------------------------------------------*/

static void TabSerialBits (ffam_Fam * fam, void *vres, void *vcho,
   void *vpar, int i, int j, int irow, int icol)
{
   int r, s, s1, L;
   long N, n;
   lebool Sparse, Over;
   smultin_Res *sres;
   fmultin_Res *fres = vres;
   fcho_Cho2 *cho = vcho;
   fcho_Cho *chon;
   fcho_Cho *chop2;
   long *Par = vpar;
   Multin_Param *kpar;

   N = Par[0];
   r = Par[1];
   s = Par[2];
   L = Par[3];
   Sparse = Par[4];
   Over = Par[5];

   util_Assert (cho, "fmultin:   cho is NULL");
   chon = cho->Chon;
   chop2 = cho->Chop2;
   util_Assert (chon, "fmultin:   cho->Chon is NULL");
   util_Assert (chop2, "fmultin:   cho->Chop2 is NULL");
   kpar = chop2->param;
   util_Assert (kpar->kcho == B_2L,
      "cho->Chop2:  wrong function for choosing number of cells");

   n = fcho_ChooseParamL (chon, 5, fmultin_Maxn, i, j);
   if (n < 0)
      return;

   strncpy (chop2->name, "L", 1);
   L = fcho_ChooseParamL (chop2, 1, 53, i, n);
   if (L < 0)
      return;

   if (Over)
      s1 = fcho_Chooses (r, s, fam->Resol[irow]);
   else
      s1 = Chooses (r, s, fam->Resol[irow], L);
   if (s1 <= 0)
      return;

   sres = smultin_CreateRes (fres->Par);
   if (Over)
      smultin_MultinomialBitsOver (fam->Gen[irow], fres->Par, sres,
         N, n, r, s1, L, Sparse);
   else
      smultin_MultinomialBits (fam->Gen[irow], fres->Par, sres,
         N, n, r, s1, L, Sparse);

   FillTables (fres, sres, N, irow, icol, Over);
   smultin_DeleteRes (sres);
}


/*========================================================================*/

void fmultin_Serial1 (ffam_Fam *fam, smultin_Param *spar,
   fmultin_Res *res, fcho_Cho2 *cho, long N, int r, int t,
   lebool Sparse, int Nr, int j1, int j2, int jstep)
{
   long Par[7];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = -1;
   Par[3] = t;
   Par[4] = Sparse;
   Par[5] = FALSE;
   Par[6] = A_MULT;

   if (spar == NULL)
      spar = &smultin_ParamDefault;
   if ((spar->GenerCell != smultin_GenerCellSerial) 
        && (spar->GenerCell != smultin_GenerCellSerial2)) {
      spar->GenerCell = smultin_GenerCellSerial;
      util_Warning (TRUE,
   "fmultin_Serial1:   changing spar->GenerCell to smultin_GenerCellSerial");
   }
   if (res == NULL) {
      localRes = TRUE;
      res = fmultin_CreateRes (spar);
   } else
      localRes = FALSE;

   PrintHead ("fmultin_Serial1", A_MULT, fam, spar, Par, Nr, j1, j2, jstep);
   InitRes (fam, res, spar, N, Nr, j1, j2, jstep, "fmultin_Serial1", FALSE);
   ftab_MakeTables (fam, res, cho, Par, TabMultin, Nr, j1, j2, jstep);
   PrintRes (res, FALSE);
   if (localRes)
      fmultin_DeleteRes (res);
}


/*========================================================================*/

void fmultin_SerialOver1 (ffam_Fam *fam, smultin_Param *spar,
   fmultin_Res *res, fcho_Cho2 *cho, long N, int r, int t,
   lebool Sparse, int Nr, int j1, int j2, int jstep)
{
   long Par[7];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = -1;
   Par[3] = t;
   Par[4] = Sparse;
   Par[5] = TRUE;
   Par[6] = A_MULT;

   if (spar == NULL)
      spar = &smultin_ParamDefault;
   if ((spar->GenerCell != smultin_GenerCellSerial) 
         && (spar->GenerCell != smultin_GenerCellSerial2)) {
      spar->GenerCell = smultin_GenerCellSerial;
      util_Warning (TRUE,
   "fmultin_SerialOver1:   changing spar->GenerCell to smultin_GenerCellSerial");
   }
   if (res == NULL) {
      localRes = TRUE;
      res = fmultin_CreateRes (spar);
   } else
      localRes = FALSE;

   PrintHead ("fmultin_SerialOver1", A_MULT, fam, spar, Par, Nr,
      j1, j2, jstep);
   InitRes (fam, res, spar, N, Nr, j1, j2, jstep, "fmultin_SerialOver1",
      TRUE);
   ftab_MakeTables (fam, res, cho, Par, TabMultin, Nr, j1, j2, jstep);
   PrintRes (res, TRUE);
   if (localRes)
      fmultin_DeleteRes (res);
}


/*========================================================================*/

void fmultin_SerialBits1 (ffam_Fam *fam, smultin_Param *spar,
   fmultin_Res *res, fcho_Cho2 *cho, long N, int r, int s, lebool Sparse,
   int Nr, int j1, int j2, int jstep)
{
   long Par[6];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = s;
   Par[3] = 0;
   Par[4] = Sparse;
   Par[5] = FALSE;

   if (spar == NULL)
      spar = &smultin_ParamDefault;
   if ((spar->GenerCell != smultin_GenerCellSerial) 
            && (spar->GenerCell != smultin_GenerCellSerial2)) {
      spar->GenerCell = smultin_GenerCellSerial;
      util_Warning (TRUE,
   "fmultin_SerialBits1:   changing spar->GenerCell to smultin_GenerCellSerial");
   }
   if (res == NULL) {
      localRes = TRUE;
      res = fmultin_CreateRes (spar);
   } else
      localRes = FALSE;

   PrintHead ("fmultin_SerialBits1", A_BITS, fam, spar, Par, Nr,
      j1, j2, jstep);
   InitRes (fam, res, spar, N, Nr, j1, j2, jstep, "fmultin_SerialBits1",
      FALSE);
   ftab_MakeTables (fam, res, cho, Par, TabSerialBits, Nr, j1, j2, jstep);
   PrintRes (res, FALSE);
   if (localRes)
      fmultin_DeleteRes (res);
}


/*========================================================================*/

void fmultin_SerialBitsOver1 (ffam_Fam *fam, smultin_Param *spar,
   fmultin_Res *res, fcho_Cho2 *cho, long N, int r, int s, lebool Sparse,
   int Nr, int j1, int j2, int jstep)
{
   long Par[6];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = s;
   Par[3] = 0;
   Par[4] = Sparse;
   Par[5] = TRUE;

   if (spar == NULL)
      spar = &smultin_ParamDefault;
   if ((spar->GenerCell != smultin_GenerCellSerial) 
            && (spar->GenerCell != smultin_GenerCellSerial2)) {
      spar->GenerCell = smultin_GenerCellSerial;
      util_Warning (TRUE,
   "fmultin_SerialBitsOver1:   changing spar->GenerCell to smultin_GenerCellSerial");
   }
   if (res == NULL) {
      localRes = TRUE;
      res = fmultin_CreateRes (spar);
   } else
      localRes = FALSE;

   PrintHead ("fmultin_SerialBitsOver1", A_BITS, fam, spar, Par, Nr,
      j1, j2, jstep);
   InitRes (fam, res, spar, N, Nr, j1, j2, jstep, "fmultin_SerialBitsOver1",
      TRUE);
   ftab_MakeTables (fam, res, cho, Par, TabSerialBits, Nr, j1, j2, jstep);
   PrintRes (res, TRUE);
   if (localRes)
      fmultin_DeleteRes (res);
}


/*========================================================================*/

void fmultin_Permut1 (ffam_Fam *fam, smultin_Param *spar,
   fmultin_Res *res, fcho_Cho2 *cho, long N, int r,
   lebool Sparse, int Nr, int j1, int j2, int jstep)
{
   long Par[7];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = 1;
   Par[3] = -1;
   Par[4] = Sparse;
   Par[5] = FALSE;
   Par[6] = A_PERM;

   if (spar == NULL)
      spar = &smultin_ParamDefault;
   if ((spar->GenerCell != smultin_GenerCellPermut)) {
      spar->GenerCell = smultin_GenerCellPermut;
      util_Warning (TRUE,
   "fmultin_Permut1:   changing GenerCell to smultin_GenerCellPermut");
   }
   if (res == NULL) {
      localRes = TRUE;
      res = fmultin_CreateRes (spar);
   } else
      localRes = FALSE;

   PrintHead ("fmultin_Permut1", A_PERM, fam, spar, Par, Nr, j1, j2, jstep);
   InitRes (fam, res, spar, N, Nr, j1, j2, jstep, "fmultin_Permut1", FALSE);
   ftab_MakeTables (fam, res, cho, Par, TabMultin, Nr, j1, j2, jstep);
   PrintRes (res, FALSE);
   if (localRes)
      fmultin_DeleteRes (res);
}


/*========================================================================*/
