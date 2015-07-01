/*************************************************************************\
 *
 * Package:        TestU01
 * File:           ufile.c
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
#include "util.h"
#include "ufile.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>
#include <limits.h>

/* Length of strings */
#define LEN 200

#define NORM32 2.3283064365386962E-10   /* 1 / 2^32 */

#define ARRAYDIM 1048576             /* = 2^20 */




/*========================== module variables ============================*/

static char S[LEN + 1];

static FILE *f1, *f2;

static int co1 = 0, co2 = 0;      /* Counters */

/* X1 will contain the numbers read from a text input file */
static double *X1 = NULL;

/* X2 will keep the bytes read from a binary input file */
static unsigned char *X2 = NULL;

static unsigned long n1, n2,      /* Current index of tables */
   MaxBin, MaxText,               /* Maximal index in the tables */
   Dim1, Dim2;                    /* Dimension of tables */

static double NBin, NText; /* Number of calls to the generator */



/*========================================================================*/

static void WrReadText (void *junk)
{
   printf (" %.0f  numbers have been read\n", NText);
}

/*-----------------------------------------------------------------------*/

static void FillTextArray (void)
/*
 * Read numbers (double's) until end of file or until Dim1 (dimension of
 * array X1) numbers have been read. 
 */
/*
 * The standard function setvbuf could be used here to increase the
 * speed of reading. Right now, we use default system buffering.
 */
{
   unsigned long i;

   MaxText = Dim1;
   i = 0;
   while ((i < Dim1) && (fscanf (f1, " %lf", (X1 + i)) == 1))
      ++i;

   if (i < MaxText)
      /* The numbers do not fill the whole array: EOF or Error */
      MaxText = i;

   n1 = 0;
}

/*-----------------------------------------------------------------------*/

static double ReadText_U01 (void *junk1, void *junk2)
/*
 * Return the n-th element of the array. If at end of array, call
 * FillTextArray to fill the array once again, and then return the
 * first element.
 */
{
   if (n1 < MaxText) {
      NText += 1.0;
      return X1[n1++];

   } else if (MaxText == Dim1) {
      FillTextArray ();
      NText += 1.0;
      return X1[n1++];

   } else {
      X1 = util_Free (X1);
      util_Fclose (f1);
      sprintf (S, "%.0f numbers have been read.\n", NText);
      strncat (S, "End-of-file detected.\n", (size_t) 25);
      strncat (S, "Not enough numbers in file for these test parameters.",
                 (size_t) 60);
      util_Error (S);
      return -1.0;
   }
}

/*-----------------------------------------------------------------------*/

static unsigned long ReadText_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (ReadText_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

unif01_Gen * ufile_CreateReadText (char *A, long dim)
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];

   util_Assert (dim > 0, "ufile_CreateReadText:   nbuf <= 0.");
   util_Assert (co1 == 0,
      "ufile_CreateReadText:   only 1 generator at a time can be in use");
   co1++;

   gen = util_Malloc (sizeof (unif01_Gen));

   strncpy (name, "ufile_CreateReadText:   ", (size_t) LEN);
   strncat (name, A, (size_t) (LEN - 30));
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   f1 = util_Fopen (A, "r");
   Dim1 = util_Min (ARRAYDIM, dim);
   MaxText = Dim1;
   X1 = util_Calloc ((size_t) Dim1, sizeof (double));
   gen->GetBits = &ReadText_Bits;
   gen->GetU01 = &ReadText_U01;
   gen->Write = &WrReadText;
   gen->param = NULL;
   gen->state = NULL;
   FillTextArray ();
   NText = 0;
   return gen;
}

/*-----------------------------------------------------------------------*/

void ufile_DeleteReadText (unif01_Gen *gen)
{
   X1 = util_Free (X1);
   util_Fclose (f1);
   gen->name = util_Free (gen->name);
   util_Free (gen);
   co1--;
}


/*-----------------------------------------------------------------------*/

void ufile_InitReadText (void)
{
   int j = 0;
   util_Assert (NULL != f1, "ufile_InitReadText:   unable to read from file");
   if (NText > Dim1) {
      j = fseek (f1, 0L, SEEK_SET);
      util_Assert (0 == j, "ufile_InitReadText:   file rewind failed");
      FillTextArray ();
   }
   NText = n1 = 0;
}


/*========================================================================*/

static void FillBinArray (void)
/*
 * Read bits, at most Dim2 bytes.
 */
{
   MaxBin = fread (X2, (size_t) 1, (size_t) Dim2, f2);
   n2 = 0;
}

/*-----------------------------------------------------------------------*/

static unsigned long ReadBin_Bits (void *vpar, void *vsta)
/*
 * Return the n-th element of the array. If at end of the array, call
 * FillBinArray to fill the array once again, and then return the first
 * element. Each number uses 32 bits in big-endian form: the first byte
 * makes the most significant bits, and the fourth one makes the least
 * significant.
 */
{
   unsigned long u;

   if (n2 < MaxBin) {
      u  = (unsigned long) X2[n2++] << 24;
      u |= (unsigned long) X2[n2++] << 16;
      u |= (unsigned long) X2[n2++] << 8;
      u |= (unsigned long) X2[n2++];
      NBin += 1.0;
      return u;

   } else if (MaxBin == Dim2) {
      FillBinArray ();
      return ReadBin_Bits (vpar, vsta);

   } else {
      X2 = util_Free (X2);
      util_Fclose (f2);
      f2 = NULL;
      sprintf (S, "%.0f bits have been read.\n", NBin * 32.0);
      strncat (S, "End-of-file detected.\n", (size_t) 25);
      strncat (S, "Not enough bits in file for these test parameters.",
               (size_t) 53);
      util_Error (S);
      return 0;
   }
}

/*-----------------------------------------------------------------------*/

static double ReadBin_U01 (void *vpar, void *vsta)
{
   return ReadBin_Bits (vpar, vsta) * NORM32;
}

/*-----------------------------------------------------------------------*/

static void WrReadBin (void *junk)
{
   printf (" %.0f  bits have been read.\n", NBin * 32.0);
}

/*-----------------------------------------------------------------------*/

unif01_Gen * ufile_CreateReadBin (char *A, long dim) 
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];

   util_Assert (dim > 0, "ufile_CreateReadBin:   nbuf <= 0.");
   util_Assert (co2 == 0,
      "ufile_CreateReadBin:   only 1 generator at a time can be in use");
   co2++;

   gen = util_Malloc (sizeof (unif01_Gen));

   strncpy (name, "ufile_CreateReadBin:   ", (size_t) LEN);
   strncat (name, A, (size_t) LEN - 30);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   f2 = util_Fopen (A, "rb");

   /* Each random number will be built of 32 bits = 4 bytes */
   Dim2 = util_Min (ARRAYDIM, 4*dim);
   X2 = util_Calloc ((size_t) Dim2, sizeof (unsigned char));
   FillBinArray ();
   NBin = 0;

   gen->GetBits = &ReadBin_Bits;
   gen->GetU01 = &ReadBin_U01;
   gen->Write = &WrReadBin;
   gen->param = NULL;
   gen->state = NULL;
   return gen;
}

/*-----------------------------------------------------------------------*/

void ufile_DeleteReadBin (unif01_Gen *gen)
{
   X2 = util_Free (X2);
   util_Fclose (f2);
   gen->name = util_Free (gen->name);
   util_Free (gen);
   co2--;
}


/*-----------------------------------------------------------------------*/

void ufile_InitReadBin (void)
{
   int j = 0;
   util_Assert (NULL != f2, "ufile_InitReadBin:   unable to read from file");
   if (NBin >= Dim2 / 4) {
      j = fseek (f2, 0L, SEEK_SET);
      util_Assert (0 == j, "ufile_InitReadBin:   file rewind failed");
      FillBinArray ();
   }
   NBin = n2 = 0;
}


/*========================================================================*/

void ufile_Gen2Bin (unif01_Gen *gen, char *fname, double nbits,
   int r, int s)
{
   unsigned long Z;
   unsigned long i, n;
   unsigned char buffer[4];
   FILE *f;
   int k;
   const int KMAX = s / 8;
   int status;

   util_Assert (nbits > 0.0, "ufile_Gen2Bin:   nbits <= 0");
   util_Assert (r >= 0, "ufile_Gen2Bin:   r < 0");
   util_Assert (s % 8 == 0,
                "ufile_Gen2Bin:   s must be in { 8, 16, 24, 32 }");
   util_Assert (nbits / s <= ULONG_MAX,
      "ufile_Gen2Bin:   nbits is too large");
   
   n = 0.5 + nbits / s;
   if (n * (double) s < nbits)
      n++;
   f = util_Fopen (fname, "wb");

   for (i = 0; i < n; i++) {
      Z = unif01_StripB (gen, r, s);
      for (k = KMAX - 1; k >= 0; k--) {
	      buffer[k] = Z & 0xFF;
	      Z >>= 8;
      }
      status = fwrite (buffer, (size_t) 1, (size_t) KMAX, f);
      if (status != KMAX) {
	 perror ("ufile_Gen2Bin:   fwrite");
	 exit (EXIT_FAILURE);
      }
   }

   util_Fclose (f);
}

