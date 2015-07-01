/*************************************************************************\
 *
 * Package:        TestU01
 * File:           ffsr.c
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

#include "ffsr.h"
#include "ffam.h"
#include "utaus.h"
#include "ulec.h"

#include <stdio.h>

#define MAXCAR  256                /* Max length of a line */
#define SEED    123                /* Default seed for generators */




/* The kinds of predefined families */
typedef enum {
   LFSR1_a,
   LFSR2_a,
   LFSR3_a,
   TausLCG2_a,
   GenN_a
}
GenType;



/*=========================================================================*/

static int ReadOneLFSR1 (char *Line, ffam_Fam * fam, int i)
{
   int status;
   unsigned int k, q, s;

   status = sscanf (Line, "%u %u %u", &k, &q, &s);
   util_Assert (status == 3,
      "ReadOneLFSR1:   Error in reading LFSR1 parameter file");
   if (k <= 32)
      fam->Gen[i] = utaus_CreateTaus (k, q, s,
             SEED & (unsigned long) (num_TwoExp[k] - 1.0));
#ifdef USE_LONGLONG
   else
      fam->Gen[i] = utaus_CreateLongTaus (k, q, s, (ulonglong) SEED);
#endif
   return 0;
}


/*=========================================================================*/

static int ReadOneLFSR2 (FILE *f, char *Line, ffam_Fam * fam, int i)
{
   int status;
   unsigned int k1, q1, s1;
   unsigned int k2, q2, s2;

   status = util_GetLine (f, Line, '#');
   if (status)                 /* if EOF or error */
      return -1;
   status = sscanf (Line, "%u %u %u", &k1, &q1, &s1);
   util_Assert (status == 3,
      "ReadOneLFSR2:   Error in reading LFSR2 parameter file");

   status = util_GetLine (f, Line, '#');
   if (status)                 /* if EOF or error */
      return -1;
   status = sscanf (Line, "%u %u %u", &k2, &q2, &s2);
   util_Assert (status == 3,
      "ReadOneLFSR2:   Error in reading LFSR2 parameter file");


   fam->Gen[i] = utaus_CreateCombTaus2 (k1, k2, q1, q2, s1, s2,
             SEED & (unsigned long) (num_TwoExp[k1] - 1.0),
             SEED & (unsigned long) (num_TwoExp[k2] - 1.0));
   return 0;
}


/*=========================================================================*/

static int ReadOneLFSR3 (FILE *f, char *Line, ffam_Fam * fam, int i)
{
   int status;
   unsigned int k1, q1, s1;
   unsigned int k2, q2, s2;
   unsigned int k3, q3, s3;

   status = util_GetLine (f, Line, '#');
   if (status)                 /* if EOF or error */
      return -1;
   status = sscanf (Line, "%u %u %u", &k1, &q1, &s1);
   util_Assert (status == 3,
      "ReadOneLFSR3:   Error in reading LFSR3 parameter file");

   status = util_GetLine (f, Line, '#');
   if (status)                 /* if EOF or error */
      return -1;
   status = sscanf (Line, "%u %u %u", &k2, &q2, &s2);
   util_Assert (status == 3,
      "ReadOneLFSR3:   Error in reading LFSR3 parameter file");

   status = util_GetLine (f, Line, '#');
   if (status)                 /* if EOF or error */
      return -1;
   status = sscanf (Line, "%u %u %u", &k3, &q3, &s3);
   util_Assert (status == 3,
      "ReadOneLFSR3:   Error in reading LFSR3 parameter file");

   fam->Gen[i] = utaus_CreateCombTaus3 (k1, k2, k3, q1, q2, q3, s1, s2, s3,
             SEED & (unsigned long) (num_TwoExp[k1] - 1.0),
             SEED & (unsigned long) (num_TwoExp[k2] - 1.0),
             SEED & (unsigned long) (num_TwoExp[k3] - 1.0));
   return 0;
}


/*=========================================================================*/

static int ReadOneTausLCG2 (FILE *f, char *Line, ffam_Fam * fam, int i)
{
   int status;
   unsigned int k1, q1, s1;
   unsigned int k2, q2, s2;
   long m, a;

   status = util_GetLine (f, Line, '#');
   if (status)                 /* if EOF or error */
      return -1;
   status = sscanf (Line, "%u %u %u %u %u %u", &k1, &q1, &s1, &k2, &q2, &s2);
   util_Assert (status == 6,
      "ReadOneTausLCG2:   Error in reading TausLCG2 parameter file");

   status = util_GetLine (f, Line, '#');
   if (status)                 /* if EOF or error */
      return -1;
   status = sscanf (Line, "%ld %ld", &m, &a);
   util_Assert (status == 2,
      "ReadOneTausLCG2:   Error in reading TausLCG2 parameter file");

   fam->Gen[i] = ulec_CreateCombTausLCG21 (
      k1, q1, s1, SEED & (unsigned long) (num_TwoExp[k1] - 1.0), 
      k2, q2, s2, SEED & (unsigned long) (num_TwoExp[k2] - 1.0),
      m, a, 0, SEED % m);
   return 0;
}


/*=========================================================================*/

static ffam_Fam *ReadAllGen (char *filename, char *deffile, GenType g,
   int i1, int i2, int istep)
/*
 * Read parameter file. If filename exists, it will be read; otherwise, the
 * default file deffile will be read. 
 * Keeps only generators whose LSize are in [i1, i2], and keep only generators
 * spaced istep apart. If i1 < the smallest LSize in the file, it will be
 * reset to the first LSize in the file; similarly if i2 > the largest LSize
 * in the file, it will be reset to the last LSize in the file.
 */
{
   FILE *f;
   char Line[MAXCAR + 1];
   ffam_Fam *fam;
   int i, j;
   int status, LSize;

   f = ffam_OpenFile (filename, deffile);

   /* Read the family name */
   util_GetLine (f, Line, '#');
   fam = ffam_CreateFam ((i2 - i1 + istep) / istep, Line);

   /* Read the parameters */
   i = j = LSize = 0;
   while (LSize <= i2) {
      status = util_GetLine (f, Line, '#');
      if (status)                 /* if EOF or error */
         break;

      status = sscanf (Line, "%d", &LSize);
      util_Assert (status == 1, "Error in reading LSize of generator");

      /* Consider only generators with (LSize >= i1) and (LSize <= i2) */
      if (LSize < i1)
         continue;
      if (LSize > i2)
         break;
      if (j++ % istep)
         continue;
      if (i >= fam->Ng)
         ffam_ReallocFam (fam, 2 * i);
      fam->LSize[i] = LSize;

      switch (g) {
      case LFSR1_a:
         fam->Resol[i] = util_Min (32, LSize);
         status = ReadOneLFSR1 (Line, fam, i);
         break;
      case LFSR2_a:
         fam->Resol[i] = util_Min (32, LSize);
         status = ReadOneLFSR2 (f, Line, fam, i);
         break;
      case LFSR3_a:
         fam->Resol[i] = util_Min (32, LSize);
         status = ReadOneLFSR3 (f, Line, fam, i);
         break;
      case TausLCG2_a:
         fam->Resol[i] = util_Min (32, LSize);
         status = ReadOneTausLCG2 (f, Line, fam, i);
         break;

      default:
         util_Error ("ReadAllGen:   impossible case");
      }

      if (status)                 /* EOF? */
         break;
      i++;
   }
   util_Assert (i > 0, "ffsr:   no generator");
   ffam_ReallocFam (fam, i);
   return fam;
}


/*=========================================================================*/

ffam_Fam * ffsr_CreateLFSR1 (char *fname, int i1, int i2, int istep)
{
   return ReadAllGen (fname, "LFSR1.par", LFSR1_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void ffsr_DeleteLFSR1 (ffam_Fam * fam)
{
   int i;
   for (i = 0; i < fam->Ng; i++)
      utaus_DeleteGen (fam->Gen[i]);
   ffam_DeleteFam (fam);
}


/*=========================================================================*/

ffam_Fam * ffsr_CreateLFSR2 (char *fname, int i1, int i2, int istep)
{
   return ReadAllGen (fname, "LFSR2.par", LFSR2_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void ffsr_DeleteLFSR2 (ffam_Fam * fam)
{
   ffsr_DeleteLFSR1 (fam);
}


/*=========================================================================*/

ffam_Fam * ffsr_CreateLFSR3 (char *fname, int i1, int i2, int istep)
{
   return ReadAllGen (fname, "LFSR3.par", LFSR3_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void ffsr_DeleteLFSR3 (ffam_Fam * fam)
{
   ffsr_DeleteLFSR1 (fam);
}


/*=========================================================================*/

ffam_Fam * ffsr_CreateTausLCG2 (char *fname, int i1, int i2, int istep)
{
   return ReadAllGen (fname, "TausLCG2.par", TausLCG2_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void ffsr_DeleteTausLCG2 (ffam_Fam * fam)
{
   int i;
   for (i = 0; i < fam->Ng; i++)
      ulec_DeleteCombTausLCG21 (fam->Gen[i]);
   ffam_DeleteFam (fam);
}
