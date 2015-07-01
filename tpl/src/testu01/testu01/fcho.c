/*************************************************************************\
 *
 * Package:        TestU01
 * File:           fcho.c
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
#include "fcho.h"
#include "ftab.h"
#include "swrite.h"

#include <stdio.h>
#include <string.h>
#include <math.h>



#define EPS2 1.0E-10
#define LEN 7




typedef struct {
   double a;
   double b;
   double c;
   fcho_FuncType F;
   char *name;
} Sample_Param;



int fcho_Resolution = 30;


/*-------------------------------- Functions ------------------------------*/



double fcho_Linear (double x)
{
   return x;
}


/*=========================================================================*/

double fcho_LinearInv (double x)
{
   return 1.0/x;
}


/*=========================================================================*/

double fcho_2Pow (double x)
{
   return pow (2.0, x);
}


/*=========================================================================*/

static void WriteSample (void *vpar, long junk, long j)
{
   Sample_Param *param = vpar;
   const double a = param->a;
   const double b = param->b;
   const double c = param->c;

   printf ("Choose  ");
   if (ftab_Style == ftab_Latex)
      printf ("$");
   if (param->name)
      printf ("%s", param->name);

   if (param->F == fcho_2Pow)
      printf (" = 2^{ ");
   else if (param->F == fcho_Linear)
      printf (" = ");
   else
      printf (" = F(");

   if (a > EPS2)
      printf ("%4.2f*i ", a);

   if (fabs (b*j) > EPS2) {
      if (b*j > EPS2)
         printf ("+ ");
      else
         printf ("- ");
      if (fabs (b - 1.0) > EPS2) 
         printf ("%4.2f*%1ld ", fabs (b), labs (j));
      else
         printf ("%1ld ", labs (j));
   }

   if (c > EPS2)
      printf ("+ %4.2f", fabs (c));
   else if (c < -EPS2)
      printf ("- %4.2f", fabs (c));

   if (param->F == fcho_2Pow)
      printf ("}");
   else if (param->F == fcho_Linear)
      ;
   else
      printf (")");

   if (ftab_Style == ftab_Latex)
      printf ("$");  
   printf ("\n\n");
}


/*-------------------------------------------------------------------------*/

static double ChooseSample (void *vpar, long i, long j)
{
   Sample_Param *param = vpar;
   double x, y;
   
   if (swrite_Basic)
      WriteSample (vpar, 0, j);
   x = i * param->a + j * param->b + param->c;
   y = param->F(x);
   return y;
}


/*-------------------------------------------------------------------------*/

fcho_Cho * fcho_CreateSampleSize (double a, double b, double c,
                                  fcho_FuncType F, char *name)
{
   fcho_Cho *cho;
   Sample_Param *param;   
   size_t len;
   char *name0 = "n";

   cho = util_Malloc (sizeof (fcho_Cho));
   param = util_Malloc (sizeof (Sample_Param));
   param->a = a;
   param->b = b;
   param->c = c;
   if (NULL == F)
      param->F = fcho_2Pow;
   else
      param->F = F;

   if (NULL == name)
      name = name0;
   len = strlen (name);
   cho->name = util_Calloc (len + 1, sizeof (char));
   strncpy (cho->name, name, (size_t) len);
   cho->param = param;
   cho->Write = WriteSample;
   cho->Choose = ChooseSample;
   param->name = cho->name;
   return cho;
}


/*-------------------------------------------------------------------------*/

void fcho_DeleteSampleSize (fcho_Cho *cho)
{
   if (NULL == cho)
      return;
   cho->name = util_Free (cho->name);
   cho->param = util_Free (cho->param);
   util_Free (cho);
}


/*=========================================================================*/

long fcho_ChooseParamL (fcho_Cho *cho, long min, long max, long i, long j)
{
   double n;

   util_Assert (cho, "fcho_ChooseParamL:   cho is NULL");
   n = cho->Choose (cho->param, i, j);

   if (n < min) {
      if (cho->name)
         printf ("%s < %ld\n\n", cho->name, min);
      return -1;
   }
   if (n > max) {
      if (cho->name)
         printf ("%s > %ld\n\n", cho->name, max);
      return -1;
   }
   return (long) n;
}


/*=========================================================================*/

int fcho_Chooses (int r, int s, int prec)
{
   int s1;

   if (r + s <= prec)
      return s;

   s1 = prec - r;
   if (s1 <= 0)
      printf ("r >= Resolution of generator\n\n");

   return s1;
}


/*=========================================================================*/

fcho_Cho2 * fcho_CreateCho2 (fcho_Cho *Chon, fcho_Cho *Chop2)
{
   fcho_Cho2 *cho;
   cho = util_Malloc (sizeof (fcho_Cho2));
   memset (cho, 0, sizeof (fcho_Cho2));
   cho->Chon = Chon;
   cho->Chop2 = Chop2;
   return cho;
}

/*-------------------------------------------------------------------------*/

void fcho_DeleteCho2 (fcho_Cho2 *cho)
{
   if (NULL == cho)
      return;
   util_Free (cho);
}

/*=========================================================================*/
