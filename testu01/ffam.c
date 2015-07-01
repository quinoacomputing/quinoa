/*************************************************************************\
 *
 * Package:        TestU01
 * File:           ffam.c
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
#include "ffam.h"
#include "unif01.h"

#include <string.h>

#define MAXCAR  256                /* Max length of a line */




/*=========================================================================*/

ffam_Fam * ffam_CreateFam (int Ng, char *name)
{
   ffam_Fam *fam;
   size_t len;

   fam = util_Malloc (sizeof (ffam_Fam));
   fam->Resol = util_Calloc ((size_t) Ng, sizeof (int));
   fam->LSize = util_Calloc ((size_t) Ng, sizeof (int));
   fam->Gen = util_Calloc ((size_t) Ng, sizeof (unif01_Gen *));
   fam->Ng = Ng;
   len = strlen (name);
   fam->name = util_Calloc (1 + len, sizeof (char));
   strncpy (fam->name, name, (size_t) len);
   return fam;
}


/*-------------------------------------------------------------------------*/

void ffam_DeleteFam (ffam_Fam *fam)
{
   if (fam == NULL)
      return;
   util_Free (fam->Resol);
   util_Free (fam->Gen);
   util_Free (fam->LSize);
   util_Free (fam->name);
   util_Free (fam);
}


/*-------------------------------------------------------------------------*/

void ffam_PrintFam (ffam_Fam *fam)
{
   int i;
   if (fam == NULL) {
      util_Warning (TRUE, "ffam_PrintFam:   fam is NULL");
      return;
   }
   printf ("-------------------------------------------------\n");
   printf ("Family:   %s\nNumber of generators:   %d\n\n",
            fam->name, fam->Ng);
   printf ("LSize Resol   Generator\n");
   printf ("-------------------------------------------------\n");
   for (i = 0; i < fam->Ng; i++) {
      printf ("%3d   %3d    %s\n", fam->LSize[i], fam->Resol[i],
         fam->Gen[i]->name);
   }
   printf ("\n\n");
}


/*=========================================================================*/

void ffam_ReallocFam (ffam_Fam * fam, int Ng)
{
   fam->Resol = util_Realloc (fam->Resol, (size_t) Ng * sizeof (int));
   fam->LSize = util_Realloc (fam->LSize, (size_t) Ng * sizeof (int));
   fam->Gen = util_Realloc (fam->Gen, (size_t) Ng * sizeof (unif01_Gen *));
   fam->Ng = Ng;
}


/*=========================================================================*/

FILE *ffam_OpenFile (char *filename, char *deffile)
/*
 * Open the file filename if it exists; otherwise if filename == NULL, open
 * deffile in the param directory; otherwise open filename in the param
 * directory.
 */
{
   FILE *f;
   char path[MAXCAR + 1];         /* Directory of parameter files */

   /* Is the parameter filename in the current directory? */
   if (filename) {
      f = fopen (filename, "r");
      if (f)
         return f;
      else
         printf ("Cannot open file  %s  in current directory."
            " Searching directory param ...\n", filename);
   }

   /* Build directory: "../param/" on Linux, "..\\param\\" on Windows */
   strncpy (path, "..", (size_t) 3);
   strncat (path, DIR_SEPARATOR, (size_t) 3);
   strncat (path, "param", (size_t) 6);
   strncat (path, DIR_SEPARATOR, (size_t) 3);

   /* Is filename NULL? open default file */
   if (filename == NULL)
      strncat (path, deffile, (size_t) MAXCAR - 20);
   else
      strncat (path, filename, (size_t) MAXCAR - 20);

   f = util_Fopen (path, "r");
   return f;
}


/*=========================================================================*/

ffam_Fam * ffam_CreateSingle (unif01_Gen *gen, int prec, int i1, int i2)
{
   int i;
   ffam_Fam *fam;

   fam = ffam_CreateFam (i2 - i1 + 1, gen->name);
   for (i = 0; i < fam->Ng; i++) {
      fam->Gen[i] = gen;
      fam->LSize[i] = i + i1;
      fam->Resol[i] = prec;
   }
   return fam;
}


void ffam_DeleteSingle (ffam_Fam *fam)
{
   ffam_DeleteFam (fam);
}


/*=========================================================================*/
