/*************************************************************************\
 *
 * Package:        TestU01
 * File:           ftab.c
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
#include "gofw.h"

#include "ftab.h"
#include "ffam.h"
#include "swrite.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>


#define MAXLEN 100                /* Max number of chars in Desc[] */



/*---------------------------- extern variables ---------------------------*/

ftab_StyleType ftab_Style = ftab_Plain;

double ftab_Suspectp = 0.01;

int ftab_SuspectLog2p = 6;




/*---------------------------- module variables ---------------------------*/

static double SuspectLog2pval;





/*-------------------------------- Functions ------------------------------*/


void ftab_SetDesc (ftab_Table *T, char *Desc)
{
   size_t len;
   util_Assert (T != NULL, "ftab_SetDesc:  ftab_Table is a NULL pointer");
   len = strlen (Desc);
   if (len > MAXLEN) {
      len = MAXLEN;
      util_Warning (1, "ftab_Table->Desc truncated");
   }
   if (T->Desc != NULL)
      T->Desc = util_Free (T->Desc);
   T->Desc = util_Calloc (len + 1, sizeof (char));
   strncpy (T->Desc, Desc, (size_t) len);
   T->Desc[len] = '\0';
}


/*=========================================================================*/

ftab_Table *ftab_CreateTable (int Nr, int j1, int j2, int jstep,
   char *Desc, ftab_FormType Form, int Ns)
{
   ftab_Table *T;
   T = util_Malloc (sizeof (ftab_Table));
   memset (T, 0, sizeof (ftab_Table));
   T->Nr = Nr;
   T->j1 = j1;
   T->j2 = j2;
   T->jstep = jstep;
   T->Nc = 1 + (j2 - j1)/jstep;
   T->Mat = tables_CreateMatrixD (T->Nr, T->Nc);
   T->LSize = util_Calloc ((size_t) T->Nr, sizeof (int));
   T->Desc = NULL;
   ftab_SetDesc (T, Desc);
   T->Form = Form;
   if (Form == ftab_String) {
      T->Strings = util_Calloc ((size_t) Ns, sizeof (char *));
      T->Ns = Ns;
   } else
      T->Strings = NULL;
   return T;
}


/*=========================================================================*/

void ftab_DeleteTable (ftab_Table * T)
{
   if (T == NULL)
      return;
   tables_DeleteMatrixD (&T->Mat);
   T->LSize = util_Free (T->LSize);
   T->Desc = util_Free (T->Desc);
   if (T->Form == ftab_String)
      T->Strings = util_Free (T->Strings);
   util_Free (T);
}


/*=========================================================================*/

void ftab_InitMatrix (ftab_Table * T, double x)
{
   int i, j;

   for (i = 0; i < T->Nr; i++)
      for (j = 0; j < T->Nc; j++)
         T->Mat[i][j] = x;
}


/*=========================================================================*/

void ftab_MakeTables (ffam_Fam *fam, void *res, void *cho, void *par,
   ftab_CalcType Calc, int Nr, int f1, int f2, int fstep)
{
   int i, j;   /* Row and column of matrices for results of one test */
   int f;
   chrono_Chrono *Timer;
   unif01_Gen *gen;

   SuspectLog2pval = 1.0 / (num_TwoExp[ftab_SuspectLog2p] - 1.0);

   Timer = chrono_Create ();

   Nr = util_Min (Nr, fam->Ng);
   for (i = 0; i < Nr; i++) {
      if (swrite_Basic) {
         printf ("CPU cumulative time: ");
         chrono_Write (Timer, chrono_hms);
         printf ("\n\n============================================="
                 "==============\n\nLSize = i = %2d\n\n", fam->LSize[i]);
      }
      if ((gen = fam->Gen[i])) {
         f = f1;
         j = 0;
         while (f <= f2) {
            Calc (fam, res, cho, par, fam->LSize[i], f, i, j);
            f += fstep;
            j++;
         }
      }
   }
   if (swrite_Basic) {
      printf ("Total CPU time: ");
      chrono_Write (Timer, chrono_hms);
      printf
         ("\n\n======================================================\n");
   }
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void PrintTexName (char *nam)
/*
 * Make sure that any _ char in Latex format name is printed as \_
 */
{
   char *p, *name = nam;
   size_t len;

   if (NULL == nam)
      return;
   len = strlen (name) + 1;
   name = util_Calloc (len, sizeof (char));
   strncpy (name, nam, (size_t) len);

   while ((p = strchr(name, '_'))) {
      *p = '\0';
      printf ("%s", name);
      printf ("\\_");
      name = p + 1;
   }
   printf ("%s", name);
}


/*=========================================================================*/

static void PrintLog2 (double d)
/*
 * Prints the logarithm (rounded) of d in base 2, when d is outside the
 * interval [SuspectLog2pval, 1 - SuspectLog2pval]; otherwise prints 
 * nothing.
 */
{
   int s;

   if (d <= gofw_Epsilonp) {
      printf ("    inf    ");
   } else if (d <= SuspectLog2pval) {
      s = 0.5 - num_Log2 (d);
      printf ("     %2d    ", s);
   } else if (d >= 1.0 - gofw_Epsilonp1) {
      printf ("   -inf    ");
   } else if (d >= 1.0 - SuspectLog2pval) {
      s = 0.5 - num_Log2 (1.0 - d);
      if (s > 9)
         printf ("    ");
      else
         printf ("     ");
      printf ("-%1d    ", s);
   } else
      printf ("           ");
}


/*-------------------------------------------------------------------------*/

static void PrintLog2Tex (double d)
/*
 * Similar to PrintLog2, but prints in Latex style.
 */
{
   int s;
   if (d <= gofw_Epsilonp) {
      printf (" & $\\infty$  ");
   } else if (d <= SuspectLog2pval) {
      s = 0.5 - num_Log2 (d);
      printf (" &  %3d   ", s);
   } else if (d >= 1.0 - gofw_Epsilonp1) {
      printf (" & $-\\infty$ ");
   } else if (d >= 1.0 - SuspectLog2pval) {
      s = 0.5 - num_Log2 (1.0 - d);
      if (s > 9)
         printf (" &  $-");
      else
         printf (" &   $-");
      printf ("%1d $ ", s);
   } else
      printf (" &        ");
}


/*=========================================================================*/

static void PrintLog10 (double d)
/*
 * Prints the logarithm (rounded) of d in base 10, when d is outside the
 * interval [ftab_Suspectp, 1 - ftab_Suspectp]; otherwise prints 
 * nothing.
 */
{
   int s;
   if (d <= gofw_Epsilonp) {
      printf ("    inf   ");
   } else if (d <= ftab_Suspectp) {
      s = 0.5 - log10 (d);
      printf ("     %2d    ", s);
   } else if (d >= 1.0 - gofw_Epsilonp1) {
      printf ("   -inf   ");
   } else if (d >= 1.0 - ftab_Suspectp) {
      s = 0.5 - log10 (1.0 - d);
      if (s > 9)
         printf ("    ");
      else
         printf ("     ");
      printf ("-%1d    ", s);
   } else {
      printf ("           ");
   }
}


/*-------------------------------------------------------------------------*/

static void PrintLog10Tex (double d)
/*
 * Similar to PrintLog10, but prints in LaTex style.
 */
{
   int s;
   if (d <= gofw_Epsilonp) {
      printf (" &  $\\infty$  ");
   } else if (d <= ftab_Suspectp) {
      s = 0.5 - log10 (d);
      printf (" &  %3d   ", s);
   } else if (d >= 1.0 - gofw_Epsilonp1) {
      printf (" & $-\\infty$ ");
   } else if (d >= 1.0 - ftab_Suspectp) {
      s = 0.5 - log10 (1.0 - d);
      if (s > 9)
         printf (" &  $-");
      else
         printf (" &   $-");
      printf ("%1d $ ", s);
   } else {
      printf (" &        ");
   }
}


/*=========================================================================*/

static void PrintVal (ftab_Table * T, double d, ftab_FormType Form)
/*
 * Prints the value d according to format Form.
 */
{
   int s;
   /* All Table tables are initialized to -1; thus the test was not done for 
      this pair (e, f) if d = -1. */
   if (d < -0.9) {
      printf ("      ---  ");
   } else if (Form == ftab_String) {
      printf ("   ");
      s = 0.5 + d;
      printf ("%s", T->Strings[s]);
   } else if (Form == ftab_Integer) {
      printf ("   ");
      if (d <= LONG_MAX)
         printf ("%8ld", (long) d);
      else
         num_WriteD (d, 8, 0, 0);
   } else if (Form == ftab_Real) {
      printf ("   ");
      num_WriteD (d, 8, 2, 2);
   } else if (Form == ftab_pLog2) {
      PrintLog2 (d);
   } else if (Form == ftab_pLog10) {
      PrintLog10 (d);
   } else if (d < gofw_Epsilonp) {
      printf ("      eps  ");
   } else if (d < ftab_Suspectp) {
      printf ("   ");
      num_WriteD (d, 8, 2, 2);
   } else if (d > 1.0 - gofw_Epsilonp1 && Form == ftab_pVal2) {
      printf ("     -eps1  ");
   } else if (d > 1.0 - ftab_Suspectp && Form == ftab_pVal2) {
      printf ("   ");
      num_WriteD (d - 1.0, 8, 2, 2);
   } else if (Form == ftab_NotInit) {
      util_Error ("ftab_PrintTable:   Form is not initialized");
   } else {
      printf ("           ");
   }
}


/*=========================================================================*/

static void PrintValTex (ftab_Table * T, double d, ftab_FormType Form)
/*
 * Similar to PrintVal, but prints in LaTex style.
 */
{
   int s;
   if (d < -0.9) {
      printf (" &   ---   ");
   } else if (Form == ftab_String) {
      printf (" & ");
      s = d + 0.5;
      printf ("%s", T->Strings[s]);
   } else if (Form == ftab_Integer) {
      printf (" & ");
      if (d <= LONG_MAX)
         printf ("%8ld", (long) d);
      else
         num_WriteD (d, 8, 0, 0);
   } else if (Form == ftab_Real) {
      printf (" & ");
      num_WriteD (d, 8, 2, 2);
   } else if (Form == ftab_pLog10) {
      PrintLog10Tex (d);
   } else if (Form == ftab_pLog2) {
      PrintLog2Tex (d);
   } else if (d < gofw_Epsilonp) {
      printf (" &   \\eps  ");
   } else if (d < ftab_Suspectp) {
      printf (" & ");
      num_WriteD (d, 8, 2, 2);
   } else if (d > 1.0 - gofw_Epsilonp1 && Form == ftab_pVal2) {
      printf (" &  \\epsm  ");
   } else if (d > 1.0 - ftab_Suspectp && Form == ftab_pVal2) {
      printf (" & ");
      num_WriteD (d - 1.0, 8, 2, 2);
   } else if (Form == ftab_NotInit) {
      util_Error ("ftab\\_PrintTable:   Form is not initialized");
   } else {
      printf (" &         ");
   }
}


/*=========================================================================*/

static void PrintTablePlain (ftab_Table * T)
/*
 * Prints table T in plain text style, according to format Form.
 */
{
   int i, j;
   int j1 = T->j1;
   int j2 = T->j2;
   int jstep = T->jstep;
   double d;
   ftab_FormType Form = T->Form;

   printf ("%s", T->Desc);
   printf ("\n\nLSize   j =%2d", j1);
   j = j1 + jstep;
   while (j <= j2) {
      printf ("      j =%2d", j);
      j += jstep;
   }
   printf ("\n------------------------------------------------------\n");

   for (i = 0; i < T->Nr; i++) {
      printf ("%3d", T->LSize[i]);
      for (j = 0; j < T->Nc; j++) {
         d = T->Mat[i][j];
         PrintVal (T, d, Form);
      }
      printf ("\n");
   }
   printf ("\n=======================================================\n");
}


/*=========================================================================*/

static void PrintTableTex (ftab_Table * T)
/*
 * Prints table T in Latex style, according to format Form.
 */
{
   int i, j;
   int j1 = T->j1;
   int j2 = T->j2;
   int jstep = T->jstep;
   ftab_FormType Form = T->Form;

   printf ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
           "\\begin {tabular}{|c|@{\\extracolsep{10pt}}");
   j = j1;
   while (j <= j2) {
      printf ("c");
      j += jstep;
   }
   printf ("|}\n\\multicolumn{%1d", (j2 - j1) / jstep + 2);
   printf ("}{l}{\\makebox[0pt][l]{");
   PrintTexName (T->Desc);
   printf ("}}\\\\\n\\hline\nLSize & $ j=%2d", j1);
   j = j1 + jstep;
   while (j <= j2) {
      printf (" $ & $ j=%2d", j);
      j += jstep;
   }
   printf ("$  \\\\\n\\hline\n");

   for (i = 0; i < T->Nr; i++) {
      printf ("%3d  ", T->LSize[i]);
      for (j = 0; j < T->Nc; j++) {
         PrintValTex (T, T->Mat[i][j], Form);
      }
      printf (" \\\\\n");
   }
   printf ("\\hline\n\\end {tabular} \\\\\n\\medskip\n\n");
}


/*=========================================================================*/

void ftab_PrintTable (ftab_Table * T)
{
   if (NULL == T)
      return;
   if (ftab_Style == ftab_Plain)
      PrintTablePlain (T);
   else
      PrintTableTex (T);
}


/*=========================================================================*/

static void PrintTable2Tex (ftab_Table * T1, ftab_Table * T2, lebool Flag)
/*
 * Prints tables in Latex style, T1 according to format Form1,
 * T2 according to format Form2.
 */
{
   int i, j;
   int j1 = T1->j1;
   int j2 = T1->j2;
   int jstep = T1->jstep;
   double x;
   ftab_FormType Form1 = T1->Form;
   ftab_FormType Form2 = T2->Form;

   printf ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
           "\\begin {tabular}{|c|@{\\extracolsep{10pt}}");
   j = j1;
   while (j <= j2) {
      printf ("rr|");
      j += jstep;
   }
   printf ("}\n\\multicolumn{%1d", 2 * ((j2 - j1) / jstep + 1) + 1);
   printf ("}{l}{\\makebox[0pt][l]{");
   PrintTexName (T1->Desc);
   printf ("---");
   PrintTexName (T2->Desc);
   if (Flag)
      printf (" (RATIO)");
   printf ("}}\\\\\n\\hline\n" " LSize& \\multicolumn{2}{c|}{$  j=%1d $}", j1);
   j = j1 + jstep;
   while (j <= j2) {
      printf (" & \\multicolumn{2}{c|}{$  j=%1d $}", j);
      j += jstep;
   }
   printf ("  \\\\\n\\hline\n");

   for (i = 0; i < T1->Nr; i++) {
      printf ("%3d", T1->LSize[i]);
      for (j = 0; j < T1->Nc; j++) {
         PrintValTex (T1, T1->Mat[i][j], Form1);
         x = T2->Mat[i][j];
         if (!Flag || x < -0.9)
            PrintValTex (T2, x, Form2);
         else {
            x = x / T1->Mat[i][j];
            PrintValTex (T2, x, ftab_Real);
         }
      }
      printf (" \\\\\n");
   }
   printf ("\\hline\n\\end {tabular} \\\\\n\\medskip\n\n");
}


/*=========================================================================*/

static void PrintTable2Plain (ftab_Table * T1, ftab_Table * T2, lebool Flag)
/*
 * Prints tables in plain text style, T1 according to format Form1,
 * T2 according to format Form2.
 */
{
   int i, j;
   int j1 = T1->j1;
   int j2 = T1->j2;
   int jstep = T1->jstep;
   double x;
   ftab_FormType Form1 = T1->Form;
   ftab_FormType Form2 = T2->Form;

   printf ("%s", T1->Desc);
   printf ("---");
   printf ("%s", T2->Desc);
   if (Flag)
      printf (" (RATIO)");
   printf ("\n\n  LSize   j=%1d", j1);
   printf ("       j=%2d", j1);
   j = j1 + jstep;
   while (j <= j2) {
      printf ("       j=%2d", j);
      printf ("       j=%2d", j);
      j += jstep;
   }
   printf ("\n----------------------------------------------------\n");

   for (i = 0; i < T1->Nr; i++) {
      printf ("%3d", T1->LSize[i]);
      for (j = 0; j < T1->Nc; j++) {
         PrintVal (T1, T1->Mat[i][j], Form1);
         x = T2->Mat[i][j];
         if (!Flag || x < -0.9)
            PrintVal (T2, x, Form2);
         else {
            x = x / T1->Mat[i][j];
            PrintVal (T2, x, ftab_Real);
         }
      }
      printf ("\n");
   }
   printf ("\n=======================================================\n");
}


/*=========================================================================*/

void ftab_PrintTable2 (ftab_Table * T1, ftab_Table * T2, lebool Flag)
{
   if (NULL == T1 || NULL == T2)
      return;
   if (ftab_Style == ftab_Plain)
      PrintTable2Plain (T1, T2, Flag);
   else
      PrintTable2Tex (T1, T2, Flag);
}
