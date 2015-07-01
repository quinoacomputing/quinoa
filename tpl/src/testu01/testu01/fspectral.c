/*************************************************************************\
 *
 * Package:        TestU01
 * File:           fspectral.c
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


#include "fspectral.h"
#include "ffam.h"
#include "fres.h"
#include "fcho.h"
#include "ftab.h"
#include "sspectral.h"

#include "num.h"

long fspectral_Maxn = 4194304;




/*------------------------------- Functions ------------------------------*/


static void PrintHead (char *name, ffam_Fam *fam, void *par1,
   int Nr, int j1, int j2, int jstep)
{
   int *Par = par1;

   printf
   ("\n\n================================================================\n");
   printf ("Family:  %s\n\n", fam->name);
   printf ("Test:    %s\n", name);
   printf ("   k = %d,   r = %d,   s = %d", Par[0], Par[1], Par[2]);
   printf ("\n   Nr = %d,   j1 = %d,   j2 = %d,   jstep = %d\n\n",
      Nr, j1, j2, jstep);
}


/*=========================================================================*/

static void TabFourier3 (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, s, k;
   long N;
   int *Par = par1;
   fres_Cont *fres = res1;
   sspectral_Res *sres;

   k = Par[0];
   r = Par[1];
   s = Par[2];

   N = fcho_ChooseParamL (cho, (long) num_TwoExp[k], fspectral_Maxn, i, j);
   if (N < 0)
      return;
   s = fcho_Chooses (r, s, fam->Resol[irow]);
   if (s <= 0)
      return;


   sres = sspectral_CreateRes ();
   sspectral_Fourier3 (fam->Gen[irow], sres, N, k, r, s);
   fres_FillTableEntryC (fres, sres->Bas->pVal2, N, irow, icol);
   sspectral_DeleteRes (sres);
}


/*========================================================================*/

void fspectral_Fourier3 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                         int k, int r, int s,
                         int Nr, int j1, int j2, int jstep)
{
   int Par[3];
   lebool localRes;

   Par[0] = k;
   Par[1] = r;
   Par[2] = s;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fspectral_Fourier3", fam, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, 2, Nr, j1, j2, jstep, "fspectral_Fourier3");
   ftab_MakeTables (fam, res, cho, Par, TabFourier3, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}
