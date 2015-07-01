/*************************************************************************\
 *
 * Package:        TestU01
 * File:           uautomata.c
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
#include "num.h"
#include "addstr.h"

#include "uautomata.h"
#include "unif01.h"

#include <stdlib.h>
#include <string.h>

/*
#include <stdio.h>
*/


#define  LEN  200                  /* Max length of strings */


typedef struct {
   int *F;                         /* Rule */
   long L;                         /* Dimension of F */
   int r;                          /* Radius of neighborhoud */
   int k;                          /* Number of cells whose values is used to
                                      build the random numbers */
   int imin, imax;                 /* Minimal and maximal cell numbers used to
      build the random numbers: i.e. the k cells numbered imin, imin + cstep,
      imin + 2*cstep, ... , imax are the ones used. */
   int cstep;                      /* Cell spacing */
   int tstep;                      /* Time spacing */
   int rot;                        /* Shift parameter */
} CA1_param;


typedef struct {
   int *Cell;
   int *OldCell;
   int N;                          /* Number of cells */
   unsigned long *Rand;            /* Random numbers */
   int nk;                         /* Index of Rand */
} CA1_state;


/*-------------------------------------------------------------------------*/

typedef struct {
   int *Cell;
   int *OldCell;
   int m;                          /* Number of cells is m + 2 */
} CA90mp_state;


/*-------------------------------------------------------------------------*/




/*-----------------------------------------------------------------------*/

static void DoOneTimeStep (CA1_param *param, CA1_state *state)
{
   unsigned long ind;
   int i;
   int s;
   int *p;
   p = state->Cell;
   state->Cell = state->OldCell;
   state->OldCell = p;

   for (i = param->r; i < state->N - param->r; i++) {
      ind = 0;
      for (s = -param->r; s <= param->r; s++) {
         ind <<= 1;
	 ind += state->OldCell[i + s];
      }
      state->Cell[i] = param->F[ind];
   }

   /* End of grid */
   for (i = state->N - param->r; i < state->N; i++) {
      ind = 0;
      for (s = -param->r; s <= param->r; s++) {
         ind <<= 1;
	 ind += state->OldCell[(i + s) % state->N];
      }
      state->Cell[i] = param->F[ind];
   }
 
   /* Beginning of grid */
   for (i = 0; i < param->r; i++) {
      ind = 0;
      for (s = -param->r; s <= param->r; s++) {
         ind <<= 1;
	 ind += state->OldCell[(i + s + state->N) % state->N];
      }
      state->Cell[i] = param->F[ind];
   }
}


/*-----------------------------------------------------------------------*/

static void RotateCells (CA1_param *param, CA1_state *state)
{
   int i;
   int *p;
   p = state->Cell;
   state->Cell = state->OldCell;
   state->OldCell = p;

   if (param->rot > 0) {
      for (i = 0; i < state->N - param->rot; i++)
	 state->Cell[i + param->rot] = state->OldCell [i];
      for (i = state->N - param->rot; i < state->N; i++)
	 state->Cell[i + param->rot - state->N] = state->OldCell [i];
 
   } else if (param->rot < 0) {
      for (i = -param->rot; i < state->N; i++)
	 state->Cell[param->rot + i] = state->OldCell [i];
      for (i = 0; i < -param->rot; i++)
	 state->Cell[state->N + param->rot + i] = state->OldCell [i];

   } else {
      p = state->Cell;
      state->Cell = state->OldCell;
      state->OldCell = p;
   }
}


/*-----------------------------------------------------------------------*/

static unsigned long CA1_Bits (void *vpar, void *vsta)
{
   CA1_state *state = vsta;

   if (state->nk <= 0) {
      CA1_param *param = vpar;
      int i, j, s, t;
      for (j = 0; j < param->k; ++j)
         state->Rand[j] = 0;

      /* Build k 32-bits random numbers from the bits of the k chosen cells */
      for (s = 0; s < 32; ++s) {
	 for (t = 0; t < param->tstep; ++t) {
	    DoOneTimeStep (param, state);
            if (param->rot != 0) {
               RotateCells (param, state);
	    }
	 }
         /* Get the bits at the k chosen cells and put them in Rand */
         j = 0;
	 for (i = param->imin; i <= param->imax; i += param->cstep) {
	    state->Rand[j] <<= 1;
            state->Rand[j] |= state->Cell[i];
            ++j;
	 }
      }
      state->nk = param->k;
   }
   return state->Rand[--state->nk];
}


/*-----------------------------------------------------------------------*/

static double CA1_U01 (void *vpar, void *vsta)
{
   return CA1_Bits (vpar, vsta) * unif01_INV32;
}


/*-----------------------------------------------------------------------*/

static void WrCA1 (void *vsta)
{
   CA1_state *state = vsta;
   int j;

   if (unif01_WrLongStateFlag) {
      printf (" S = {\n ");
      for (j = 0; j < state->N; j++)
         printf (" %1d", state->Cell[j]);
      printf ("    }\n");
   } else
      unif01_WrLongStateDef ();
}


/*-----------------------------------------------------------------------*/

unif01_Gen * uautomata_CreateCA1 (int N, int S[ ], int r, int F[ ], 
                                  int k, int ts, int cs, int rot)
{
   unif01_Gen *gen;
   CA1_param *param;
   CA1_state *state;
   size_t leng;
   char name[LEN + 1];
   long i;
   unsigned long Z = 0;
   long L;                 /* Dimension of the rule F */

   i = 4*sizeof(long) - 1;
   util_Assert ((long) r < i, "uautomata_CreateCA1:   r too large");
   util_Assert (ts >= 0, "uautomata_CreateCA1:   ts < 0");
   util_Assert (cs >= 0, "uautomata_CreateCA1:   cs < 0");
   util_Assert (k > 0,  "uautomata_CreateCA1:   k < 1");
   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (CA1_param));
   state = util_Malloc (sizeof (CA1_state));

   strncpy (name, "uautomata_CreateCA1:", (size_t) LEN);
   addstr_Long (name, "   N = ", (long) N);
   addstr_Long (name, ",   r = ", (long) r);
   L = num_TwoExp[2*r + 1];
   addstr_ArrayInt (name, ",   F = ",  (int) L, F);
   for (i = L-1; i >= 0; i--) {
      Z <<= 1;
      Z += F[i];
   }
   addstr_Ulong(name, " = Rule ", Z);   
   addstr_Long (name, ",   k = ", (long) k);
   addstr_Long (name, ",   ts = ", (long) ts);
   addstr_Long (name, ",   cs = ", (long) cs);
   addstr_Long (name, ",   rot = ", (long) rot);
   addstr_ArrayInt (name, ",   S = ",  N, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, (size_t) leng);

   param->r = r;
   param->tstep = ts + 1;
   param->cstep = cs + 1;
   param->rot = rot;
   param->k = k;
   param->F = F;
   param->L = L;
   state->Cell = util_Calloc ((size_t) N, sizeof (int));
   state->OldCell = util_Calloc ((size_t) N, sizeof (int));
   state->Rand = util_Calloc ((size_t) k, sizeof (unsigned long));
   state->N = N;
   state->nk = 0;
   param->imin = N/2 - k/2 * param->cstep;
   param->imax = N/2 + (k - 1)/2 * param->cstep;
   util_Assert (param->imin >= 0, "uautomata_CreateCA1:   k*cs too large");
   util_Assert (param->imax < N, "uautomata_CreateCA1:   k*cs too large");

   for (i = 0; i < N; i++) {
     /*     util_Assert (S[i] == 0 || S[i] == 1,
            "uautomata_CreateCA1:   all S[i] must be in { 0, 1 }."); */
      state->Cell[i] = S[i] & 1;
   }
   
   gen->GetBits = &CA1_Bits;
   gen->GetU01  = &CA1_U01;
   gen->Write   = &WrCA1;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

static unsigned long CA90mp_Bits (void *junk, void *vsta)
{
   CA90mp_state *state = vsta;
   unsigned long Rand = 0;             /* Random number */
   int s, i;
   int *p;

   /* Build a 32-bit random integer from the successive bit of cell m */
   for (s = 0; s < 32; ++s) {
      p = state->Cell;
      state->Cell = state->OldCell;
      state->OldCell = p;

      /* Do one time step */
      for (i = 1; i <= state->m; i++)
	 state->Cell[i] = state->OldCell[i - 1] ^ state->OldCell[i + 1];
      /* Mirror boundary condition */
      state->Cell[state->m + 1] = state->Cell[state->m];

      Rand <<= 1;
      Rand |= state->Cell[state->m];
   }

   return Rand;
}


/*-----------------------------------------------------------------------*/

static double CA90mp_U01 (void *vpar, void *vsta)
{
   return CA90mp_Bits (vpar, vsta) * unif01_INV32;
}


/*-----------------------------------------------------------------------*/

static void WrCA90mp (void *vsta)
{
   CA90mp_state *state = vsta;
   int j;

   if (unif01_WrLongStateFlag) {
      printf (" S = {\n ");
      for (j = 1; j <= state->m; j++)
         printf (" %1d", state->Cell[j]);
      printf ("    }\n");
   } else
      unif01_WrLongStateDef ();
}


/*-----------------------------------------------------------------------*/

unif01_Gen * uautomata_CreateCA90mp (int m, int S[])
{
   unif01_Gen *gen;
   CA90mp_state *state;
   size_t leng;
   char name[LEN + 1];
   int i;

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (CA90mp_state));

   strncpy (name, "uautomata_CreateCA90mp:", (size_t) LEN);
   addstr_Long (name, "   m = ", (long) m);
   addstr_ArrayInt (name, ",   S = ",  m, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, (size_t) leng);

   state->Cell = util_Calloc ((size_t) m + 2, sizeof (int));
   state->OldCell = util_Calloc ((size_t) m + 2, sizeof (int));
   state->m = m;

   for (i = 1; i <= m; i++) {
      util_Assert (S[i] == 0 || S[i] == 1,
        "uautomata_CreateCA90mp:   all S[i] must be in { 0, 1 }.");
      state->Cell[i] = S[i];
   }
   /* Null boundary condition */
   state->Cell[0] = 0;
   state->OldCell[0] = 0;
   gen->GetBits = &CA90mp_Bits;
   gen->GetU01  = &CA90mp_U01;
   gen->Write   = &WrCA90mp;
   gen->state   = state;
   gen->param   = NULL;
   return gen;
}


/**************************************************************************/

void uautomata_DeleteGen (unif01_Gen *gen)
{
   CA1_state *state;
   if (NULL == gen)
      return;
   state = gen->state;
   util_Free (state->Cell);
   util_Free (state->OldCell);
   util_Free (state->Rand);
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


void uautomata_DeleteCA90mp (unif01_Gen *gen)
{
   CA90mp_state *state;
   if (NULL == gen)
      return;
   state = gen->state;
   util_Free (state->Cell);
   util_Free (state->OldCell);
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}

