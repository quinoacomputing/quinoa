/*************************************************************************\
 *
 * Package:        TestU01
 * File:           fres.c
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
#include "bitset.h"
#include "fres.h"
#include "ftab.h"

#include <string.h>

#define LEN 100




/*=========================================================================*/

void fres_InitCont (ffam_Fam *fam, fres_Cont *res, int N,
   int Nr, int f1, int f2, int fstep, char *nam)
{
   int i, j;
   char str[LEN + 1] = {0};
   size_t len1;
   char *p;

   res->name = util_Realloc (res->name, 1 + strlen (nam) * sizeof (char));
   strcpy (res->name, nam);
   Nr = util_Min (Nr, fam->Ng);
   res->Active = 0;

   for (j = 0; j < gofw_NTestTypes; j++) {
      if ((gofw_Mean == j) ||
            (N > 1 && (bitset_TestBit (gofw_ActiveTests, j)))) {
         strncpy (str, nam, (size_t) LEN);
         len1 = strlen (str);
	 strncat (str, ": ", 3);
         p = strstr (res->PVal[j]->Desc, "p-value");
         if (p)
	    strncat (str, p, (size_t) LEN - len1);
         ftab_DeleteTable (res->PVal[j]);
	 res->PVal[j] = ftab_CreateTable (Nr, f1, f2, fstep, str,
                        ftab_pVal2, 0);
	 ftab_InitMatrix (res->PVal[j], -1.0);
         bitset_SetBit (res->Active, j);
         for (i = 0; i < Nr; i++)
	    res->PVal[j]->LSize[i] = fam->LSize[i];
      }
   }
   if (N > 1)
      bitset_ClearBit (res->Active, gofw_Mean);   
}


/*-------------------------------------------------------------------------*/

fres_Cont * fres_CreateCont (void)
{
   fres_Cont *res;
   char str[LEN + 1];
   gofw_TestType j;
   size_t m;

   res = util_Malloc (sizeof (fres_Cont));
   res->name = util_Calloc (1, sizeof (char));

   m = strlen ("p-value for statistic ");
   for (j = 0; j < gofw_NTestTypes; j++) {
      if ((gofw_Mean == j) || (bitset_TestBit (gofw_ActiveTests, j))) {
	 strncpy (str, "p-value for ", (size_t) LEN);
	 if (gofw_Mean != j)
            strncat (str, gofw_TestNames[j], (size_t) LEN - m);
	 strncat (str, " statistic", (size_t) LEN - m);
	 res->PVal[j] = ftab_CreateTable (1, 0, 1, 1, str, ftab_pVal2, 0);
      }
   }
   return res;
}


/*-------------------------------------------------------------------------*/

void fres_DeleteCont (fres_Cont *res)
{
   gofw_TestType j;

   if (res == NULL)
      return;
   res->name = util_Free (res->name);

   for (j = 0; j < gofw_NTestTypes; j++) {
      if ((gofw_Mean == j) || (bitset_TestBit (gofw_ActiveTests, j))) {
         ftab_DeleteTable (res->PVal[j]);
         res->PVal[j] = NULL;
      }
   }
   util_Free (res);
}


/*=========================================================================*/

void fres_PrintCont (fres_Cont *res)
{
   gofw_TestType j;

   for (j = 0; j <= gofw_Mean; j++) {
      if (bitset_TestBit (res->Active, j))
         ftab_PrintTable (res->PVal[j]);
   }
}


/*=========================================================================*/

void fres_FillTableEntryC (fres_Cont *fres, gofw_TestArray pval,
   int N, int i, int j)
/*
 * Writes the results of one test in the tables.
 */
{
   gofw_TestType k;

   if (N == 1) {
      fres->PVal[gofw_Mean]->Mat[i][j] = pval[gofw_Mean];

   } else {
      for (k = 0; k <= gofw_Mean; k++) {
         if (bitset_TestBit (gofw_ActiveTests, k)) {
            fres->PVal[k]->Mat[i][j] = pval[k];
         }
      }
   }
}


/*=========================================================================*/

void fres_InitDisc (ffam_Fam *fam, fres_Disc *res,
   int Nr, int f1, int f2, int fstep, char *nam)
{
   char str[LEN + 1] = {0};
   char str2[LEN + 1] = {0};
   size_t len1;
   int i;

   res->name = util_Realloc (res->name, 1 + strlen (nam) * sizeof (char));
   strcpy (res->name, nam);
   strncpy (str, nam, (size_t) LEN);
   len1 = strlen (nam);
 
   Nr = util_Min (Nr, fam->Ng);

   ftab_DeleteTable (res->PVal2);
   ftab_DeleteTable (res->PRight);
   ftab_DeleteTable (res->PLeft);

   strncpy (str2, nam, (size_t) LEN);
   strncat (str2, ", Left p-value", (size_t) LEN - len1);
   res->PLeft = ftab_CreateTable (Nr, f1, f2, fstep, str2, ftab_pVal1, 0);

   strncpy (str2, nam, (size_t) LEN);
   strncat (str2, ", Right p-value", (size_t) LEN - len1);
   res->PRight = ftab_CreateTable (Nr, f1, f2, fstep, str2, ftab_pVal1, 0);

   strncpy (str2, nam, (size_t) LEN);
   strncat (str2, ", p-value for discrete statistic", (size_t) LEN - len1);
   res->PVal2 = ftab_CreateTable (Nr, f1, f2, fstep, str2, ftab_pVal2, 0);

   ftab_InitMatrix (res->PLeft, -1.0);
   ftab_InitMatrix (res->PRight, -1.0);
   ftab_InitMatrix (res->PVal2, -1.0);

   for (i = 0; i < Nr; i++) {
      res->PLeft->LSize[i] = fam->LSize[i];
      res->PRight->LSize[i] = fam->LSize[i];
      res->PVal2->LSize[i] = fam->LSize[i];
   }
}


/*-------------------------------------------------------------------------*/

fres_Disc * fres_CreateDisc (void)
{
   fres_Disc *res;

   res = util_Malloc (sizeof (fres_Disc));
   res->name = util_Calloc (1, sizeof (char));

   res->PLeft = ftab_CreateTable (1, 0, 1, 1, "", ftab_pVal1, 0);
   res->PRight = ftab_CreateTable (1, 0, 1, 1, "", ftab_pVal1, 0);
   res->PVal2 = ftab_CreateTable (1, 0, 1, 1, "", ftab_pVal2, 0);
   return res;
}


/*-------------------------------------------------------------------------*/

void fres_DeleteDisc (fres_Disc *res)
{
   if (res == NULL)
      return;
   res->name = util_Free (res->name);
   ftab_DeleteTable (res->PVal2);
   ftab_DeleteTable (res->PRight);
   ftab_DeleteTable (res->PLeft);
   util_Free (res);
}


/*=========================================================================*/

void fres_PrintDisc (fres_Disc *res, lebool LR)
{
   if (LR) {
      ftab_PrintTable (res->PLeft);
      ftab_PrintTable (res->PRight);
   }
   ftab_PrintTable (res->PVal2);
}


/*=========================================================================*/

void fres_FillTableEntryD (fres_Disc *fres,
   double pLeft, double pRight, double pVal2, int i, int j)
/*
 * Writes the results of one test in the tables.
 */
{
   fres->PLeft->Mat[i][j] = pLeft;
   fres->PRight->Mat[i][j] = pRight;
   fres->PVal2->Mat[i][j] = pVal2;
}


/*=========================================================================*/

void fres_InitPoisson (ffam_Fam *fam, fres_Poisson *res,
   int Nr, int f1, int f2, int fstep, char *nam)
{
   char str[LEN + 1] = {0};
   char str2[LEN + 1] = {0};
   size_t len1;
   int i;

   res->name = util_Realloc (res->name, 1 + strlen (nam) * sizeof (char));
   strcpy (res->name, nam);
   strncpy (str, nam, (size_t) LEN);
   len1 = strlen (nam);
 
   Nr = util_Min (Nr, fam->Ng);

   ftab_DeleteTable (res->Obs);
   ftab_DeleteTable (res->Exp);
   ftab_DeleteTable (res->PVal2);
   ftab_DeleteTable (res->PRight);
   ftab_DeleteTable (res->PLeft);

   strncpy (str2, nam, (size_t) LEN);
   strncat (str2, ", Expected numbers", (size_t) LEN - len1);
   res->Exp = ftab_CreateTable (Nr, f1, f2, fstep, str2, ftab_Real, 0);

   strncpy (str2, nam, (size_t) LEN);
   strncat (str2, ", Observed numbers", (size_t) LEN - len1);
   res->Obs = ftab_CreateTable (Nr, f1, f2, fstep, str2, ftab_Integer, 0);

   strncpy (str2, nam, (size_t) LEN);
   strncat (str2, ", Left p-value", (size_t) LEN - len1);
   res->PLeft = ftab_CreateTable (Nr, f1, f2, fstep, str2, ftab_pVal1, 0);

   strncpy (str2, nam, (size_t) LEN);
   strncat (str2, ", Right p-value", (size_t) LEN - len1);
   res->PRight = ftab_CreateTable (Nr, f1, f2, fstep, str2, ftab_pVal1, 0);

   strncpy (str2, nam, (size_t) LEN);
   strncat (str2, ", p-value for discrete statistic", (size_t) LEN - len1);
   res->PVal2 = ftab_CreateTable (Nr, f1, f2, fstep, str2, ftab_pVal2, 0);

   ftab_InitMatrix (res->Exp, -1.0);
   ftab_InitMatrix (res->Obs, -1.0);
   ftab_InitMatrix (res->PLeft, -1.0);
   ftab_InitMatrix (res->PRight, -1.0);
   ftab_InitMatrix (res->PVal2, -1.0);

   for (i = 0; i < Nr; i++) {
      res->PLeft->LSize[i] = fam->LSize[i];
      res->PRight->LSize[i] = fam->LSize[i];
      res->PVal2->LSize[i] = fam->LSize[i];
      res->Exp->LSize[i] = fam->LSize[i];
      res->Obs->LSize[i] = fam->LSize[i];
   }
}


/*-------------------------------------------------------------------------*/

fres_Poisson * fres_CreatePoisson (void)
{
   fres_Poisson *res;

   res = util_Malloc (sizeof (fres_Poisson));
   res->name = util_Calloc (1, sizeof (char));

   res->Obs = ftab_CreateTable (1, 0, 1, 1, "", ftab_pVal1, 0);
   res->Exp = ftab_CreateTable (1, 0, 1, 1, "", ftab_pVal1, 0);
   res->PLeft = ftab_CreateTable (1, 0, 1, 1, "", ftab_pVal1, 0);
   res->PRight = ftab_CreateTable (1, 0, 1, 1, "", ftab_pVal1, 0);
   res->PVal2 = ftab_CreateTable (1, 0, 1, 1, "", ftab_pVal2, 0);
   return res;
}


/*-------------------------------------------------------------------------*/

void fres_DeletePoisson (fres_Poisson *res)
{
   if (res == NULL)
      return;
   res->name = util_Free (res->name);
   ftab_DeleteTable (res->PVal2);
   ftab_DeleteTable (res->PRight);
   ftab_DeleteTable (res->PLeft);
   ftab_DeleteTable (res->Obs);
   ftab_DeleteTable (res->Exp);
   util_Free (res);
}


/*=========================================================================*/

void fres_PrintPoisson (fres_Poisson *res, lebool LR, lebool Ratio)
{
   ftab_PrintTable2 (res->Exp, res->Obs, Ratio);
   if (LR) {
      ftab_PrintTable (res->PLeft);
      ftab_PrintTable (res->PRight);
   }
   ftab_PrintTable (res->PVal2);
}


/*=========================================================================*/

void fres_FillTableEntryPoisson (fres_Poisson *fres, double Exp, double Obs,
   double pLeft, double pRight, double pVal2, int i, int j)
/*
 * Writes the results of one test in the tables.
 */
{
   fres->Obs->Mat[i][j] = Obs;
   fres->Exp->Mat[i][j] = Exp;
   fres->PLeft->Mat[i][j] = pLeft;
   fres->PRight->Mat[i][j] = pRight;
   fres->PVal2->Mat[i][j] = pVal2;
}


/*=========================================================================*/
