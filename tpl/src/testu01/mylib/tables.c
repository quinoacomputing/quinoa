/*************************************************************************\
 *
 * Package:        MyLib
 * File:           tables.c
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


#include "tables.h"
#include "util.h"
#include "mystr.h"
#include "num.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>






tables_StyleType Style = tables_Plain;

static char OuvrantMat = ' ';     /* Matrix delimitors */
static char FermantMat = ' ';

static char OuvrantVec = ' ';     /* Vector delimitors */
static char FermantVec = ' ';

static char SepareVec = ' ';      /* Element separators */
static char SepareElem = ' ';

#define MaxInd 60
static long HacheTab[MaxInd + 1] = {
   8191, 12109, 16381, 24373, 32749, 48871, 65521, 97777, 131071, 195659,
   262139, 393203, 524287, 786407, 1048573, 1572803, 2097143, 2500097,
   3145711, 3600097, 4194301, 5300003, 6291403, 7300003, 8388593, 9500021,
   10500013, 11500003, 12582917, 13500007, 14500001, 15500011, 16777213,
   17500013, 18500017, 19500101, 20500097, 21300101, 22200001,
   23200097, 24200101, 25165807, 28000097, 30000001, 33554393, 39000001,
   45000097, 50331653, 55000013, 61000001, 67108859, 76000091, 85000007,
   94906247, 134217689, 189812501, 268435399, 379625003, 2147483647, -1
};


long **tables_CreateMatrixL (int N, int M)
{
   int i;
   long **T1;
   /* Note: the memory must be allocated in a contiguous way for the matrix
      to be later used properly, without problems. Source: comp.lang.c -
      Answers to Frequently Asked Questions
      http://www.faqs.org/faqs/C-faq/faq/ Questions 6.18, 6.19, 6.20 */

   T1 = (long **) util_Malloc (N * sizeof (long *));
   T1[0] = (long *) util_Malloc (N * M * sizeof (long));
   for (i = 1; i < N; i++)
      T1[i] = T1[0] + i * M;
   return T1;
}


unsigned long **tables_CreateMatrixUL (int N, int M)
{
   int i;
   unsigned long **T3;

   T3 = (unsigned long **) util_Malloc (N * sizeof (unsigned long *));
   T3[0] = (unsigned long *) util_Malloc (N * M * sizeof (unsigned long));
   for (i = 1; i < N; i++)
      T3[i] = T3[0] + i * M;
   return T3;
}


double **tables_CreateMatrixD (int N, int M)
{
   int i;
   double **T2;

   T2 = (double **) util_Malloc (N * sizeof (double *));
   T2[0] = (double *) util_Malloc (N * M * sizeof (double));
   for (i = 1; i < N; i++)
      T2[i] = T2[0] + i * M;
   return T2;
}


void tables_DeleteMatrixL (long ***T)
{
   free ((*T)[0]);
   free (*T);
   *T = NULL;
}

void tables_DeleteMatrixUL (unsigned long ***T)
{
   free ((*T)[0]);
   free (*T);
   *T = NULL;
}

void tables_DeleteMatrixD (double ***T)
{
   free ((*T)[0]);
   free (*T);
   *T = NULL;
}


void tables_CopyTabD (double T1[], double T2[], int n1, int n2)
{
   int i;
   for (i = n1; i <= n2; i++) {
      T2[i] = T1[i];
   }
}

void tables_CopyTabL (long T1[], long T2[], int n1, int n2)
{
   int i;
   for (i = n1; i <= n2; i++) {
      T2[i] = T1[i];
   }
}

void tables_QuickSortD (double T[], int l, int r)
   /* On trie le tableau des observations T[l..r].  */
{
   int j;                         /* Indices dans le tableau Tab.  */
   int i;
   double w;
   double x;
   i = l;
   j = r;
   x = T[(l + r) / 2];
   do {
      while (T[i] < x)
         ++i;
      while (x < T[j])
         --j;
      if (i <= j) {
         w = T[i];
         T[i] = T[j];
         T[j] = w;
         ++i;
         --j;
      }
   } while (i <= j);
   if (l < j)
      tables_QuickSortD (T, l, j);
   if (i < r)
      tables_QuickSortD (T, i, r);
}


void tables_QuickSortL (long T[], int l, int r)
     /* On trie le tableau des observations T[l..r].  */
{
   int j;                         /* Indices dans le tableau Tab.  */
   int i;
   long w;
   long x;
   i = l;
   j = r;
   x = T[(l + r) / 2];
   do {
      while (T[i] < x)
         ++i;
      while (x < T[j])
         --j;
      if (i <= j) {
         w = T[i];
         T[i] = T[j];
         T[j] = w;
         ++i;
         --j;
      }
   } while (i <= j);
   if (l < j)
      tables_QuickSortL (T, l, j);
   if (i < r)
      tables_QuickSortL (T, i, r);
}


/*=======================================================================*/
#ifdef USE_LONGLONG

void tables_QuickSortLL (longlong T[], int l, int r)
{
   int j;
   int i;
   longlong w;
   longlong x;
   i = l;
   j = r;
   x = T[(l + r) / 2];
   do {
      while (T[i] < x)
         ++i;
      while (x < T[j])
         --j;
      if (i <= j) {
         w = T[i];
         T[i] = T[j];
         T[j] = w;
         ++i;
         --j;
      }
   } while (i <= j);
   if (l < j)
      tables_QuickSortLL (T, l, j);
   if (i < r)
      tables_QuickSortLL (T, i, r);
}

void tables_QuickSortULL (ulonglong T[], int l, int r)
{
   int j;
   int i;
   ulonglong w;
   ulonglong x;
   i = l;
   j = r;
   x = T[(l + r) / 2];
   do {
      while (T[i] < x)
         ++i;
      while (x < T[j])
         --j;
      if (i <= j) {
         w = T[i];
         T[i] = T[j];
         T[j] = w;
         ++i;
         --j;
      }
   } while (i <= j);
   if (l < j)
      tables_QuickSortULL (T, l, j);
   if (i < r)
      tables_QuickSortULL (T, i, r);
}

#endif
/*=======================================================================*/

void tables_WriteTabL (long V[], int n1, int n2, int k, int p, char Desc[])
{
   int i;
   printf ("---------------------------------------\n");
   printf ("%s\n", Desc);
   if (k > 1) {
      printf ("Elements  %d  to  %d\n\n", n1, n2);
      for (i = n1; i <= n2; i++) {
         printf ("%*ld ", p, V[i]);
         if (((i + 1 - n1) % k) == 0)
            printf ("\n");
      }
      printf ("\n");
   } else {
      printf ("\n Index        Element\n");
      for (i = n1; i <= n2; i++)
         printf ("%6d   %12ld\n", i, V[i]);
   }
   printf ("\n");
}


void tables_WriteTabD (double V[], int n1, int n2, int k, int p1,
   int p2, int p3, char Desc[])
{
   int i;
   printf ("---------------------------------------\n");
   printf ("%s\n", Desc);
   if (k > 1) {
      printf ("Elements  %d  to  %d\n\n", n1, n2);
      for (i = n1; i <= n2; i++) {
         /* printf ("%*.*G", p1, p2, V[i]); */
         num_WriteD (V[i], p1, p2, p3);
         if (((i + 1 - n1) % k) == 0)
            printf ("\n");
      }
      printf ("\n");
   } else {
      printf ("\n Index            Element\n");
      for (i = n1; i <= n2; i++) {
         printf ("%6d", i);
         num_WriteD (V[i], p1, p2, p3);
         printf ("\n");
      }
   }
   printf ("\n");
}


/*=========================================================================*/
#ifdef USE_LONGLONG

void tables_WriteTabLL (longlong V[], int n1, int n2, int k, int p,
    char Desc[])
{
   int i;
   printf ("---------------------------------------\n");
   printf ("%s\n", Desc);
   if (k > 1) {
      printf ("Elements  %d  to  %d\n\n", n1, n2);
      for (i = n1; i <= n2; i++) {
         printf (" %*" PRIdLEAST64, p, V[i]);
         if (((i + 1 - n1) % k) == 0)
            printf ("\n");
      }
      printf ("\n");
   } else {
      printf ("\n Index        Element\n");
      for (i = n1; i <= n2; i++)
         printf ("%6d     %12" PRIdLEAST64 "\n", i, V[i]);
   }
   printf ("\n");
}

void tables_WriteTabULL (ulonglong V[], int n1, int n2, int k, int p,
    char Desc[])
{
   int i;
   printf ("---------------------------------------\n");
   printf ("%s\n", Desc);
   if (k > 1) {
      printf ("Elements  %d  to  %d\n\n", n1, n2);
      for (i = n1; i <= n2; i++) {
         printf (" %*" PRIuLEAST64, p, V[i]);
         if (((i + 1 - n1) % k) == 0)
            printf ("\n");
      }
      printf ("\n");
   } else {
      printf ("\n Index        Element\n");
      for (i = n1; i <= n2; i++)
         printf ("%6d     %12" PRIuLEAST64 "\n", i, V[i]);
   }
   printf ("\n");
}

#endif
/*=========================================================================*/


static void FixeDelim (tables_StyleType style)
{
   /* Fixe les delimiteurs pour imprimer une matrice selon un format
      approprie */
   Style = style;
   switch (style) {
   case tables_Mathematica:
      OuvrantMat = '{';
      FermantMat = '}';
      OuvrantVec = '{';
      FermantVec = '}';
      SepareVec = ',';
      SepareElem = ',';
      break;
   case tables_Matlab:
      OuvrantMat = '[';
      FermantMat = ']';
      OuvrantVec = ' ';
      FermantVec = ' ';
      SepareVec = ' ';
      SepareElem = ' ';
      break;
   default:
      OuvrantMat = ' ';
      FermantMat = ' ';
      OuvrantVec = ' ';
      FermantVec = ' ';
      SepareVec = ' ';
      SepareElem = ' ';
      break;
   }
}


void tables_WriteMatrixL (long **Mat, int i1, int i2, int j1, int j2,
   int w, tables_StyleType style, char Nom[])
{
   int i;
   int j;

   FixeDelim (style);
   if (strlen (Nom) > 0) {
      printf ("%s = ", Nom);
   }
   printf ("%c\n", OuvrantMat);
   for (i = i1; i <= i2; i++) {
      printf ("%c", OuvrantVec);
      for (j = j1; j <= j2; j++) {
         printf ("%*ld", (int) w, Mat[i][j]);
         if (j < j2)
            printf ("%c", SepareElem);
      }
      printf ("%c", FermantVec);
      if (i < i2)
         printf ("%c\n", SepareVec);
   }
   printf ("%c\n\n", FermantMat);
}


void tables_WriteMatrixD (double **Mat, int i1, int i2, int j1, int j2,
   int w, int p, tables_StyleType style, char Nom[])
{
   int k;
   int m;
   int j;
   int i;
   unsigned int bidon;
   double prec;
   double x;
   int trouve;
   char S[32];

   FixeDelim (style);
   if (strlen (Nom) > 0) {
      printf ("%s = ", Nom);
   }
   prec = pow (10.0, (double) p);
   printf ("%c\n", OuvrantMat);
   for (i = i1; i <= i2; i++) {
      printf ("%c", OuvrantVec);
      for (j = j1; j <= j2; j++) {
         printf (" ");
         switch (style) {
         case tables_Mathematica:
            x = Mat[i][j];
            if (((x != 0.0) && (fabs (x) < 0.1)) || (fabs (x) > prec)) {
               sprintf (S, "%.*G", (int) p, x);
               /* automatique avec %G ... : myst_Subst(S, "e", "E"); */
               mystr_Position ("E", S, 0, &bidon, &trouve);
               if (trouve) {
                  mystr_Subst (S, "E", "*10^(");
                  strncat (S, ")", (size_t) 2);
               }
            } else
               sprintf (S, "%.*f", (int) p, x);
            m = (int) strlen (S);
            for (k = 1; k <= w - m; k++) {
               printf (" ");
            }
            printf ("%s", S);
            break;
         default:
            /* tables_Matlab, Default */
            printf ("%*.*G", (int) w, (int) p, Mat[i][j]);
            break;
         }
         if (j < j2)
            printf ("%c", SepareElem);
      }
      printf ("%c", FermantVec);
      if (i < i2)
         printf ("%c\n", SepareVec);
   }
   printf ("%c\n\n", FermantMat);
}

long tables_HashPrime (long n, double load)
{
   int i;
   double nD;
   util_Assert (n > 0, "tables_HashPrime : n <= 0");
   nD = (double) n;
   i = 1;
   while (i < MaxInd && HacheTab[i] < n)
      ++i;
   while (i < MaxInd && load * (double) (HacheTab[i]) < nD)
      ++i;
   util_Assert (HacheTab[i] > 0, "tables_HashPrime failed");
   return HacheTab[i];
}
