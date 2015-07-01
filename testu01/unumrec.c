/*************************************************************************\
 *
 * Package:        TestU01
 * File:           unumrec.c
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
#include "addstr.h"

#include "unumrec.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>



/*============================== Constants ===============================*/

#define  LEN  100                 /* Max length of strings */

#define  N_TAB    32
#define  N_TAB8   40              /* N_TAB + 8 */

#define  M1  2147483647
#define  A1  16807
#define  q1  127773
#define  r1  2836
#define  NDiv  67108864           /* = 1 + (M1 - 1) / N_TAB */



/*================================ Types ================================*/

typedef struct {
   double Norm;
} Ran0_param;

typedef struct {
   long S1;
} Ran0_state;

typedef struct {
   long S1, z;
   long Tab[N_TAB + 1];
} Ran1_state;

typedef struct {
   long S1, S2, z;
   long Tab[N_TAB + 1];
} Ran2_state;


/**************************************************************************/

static double Ran0_U01 (void *vpar, void *vsta)
{
   Ran0_param *param = vpar;
   Ran0_state *state = vsta;
   long k;

   k = state->S1 / q1;
   state->S1 = A1 * (state->S1 - k * q1) - k * r1;
   if (state->S1 < 0)
      state->S1 += M1;
   return (state->S1 * param->Norm);
}

static unsigned long Ran0_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Ran0_U01 (vpar, vsta));
}

static void WrRan0 (void *vsta)
{
   Ran0_state *state = vsta;
   printf (" S = %1ld\n", state->S1);
}

unif01_Gen * unumrec_CreateRan0 (long s)
{
   unif01_Gen *gen;
   Ran0_param *param;
   Ran0_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert (s > 0, "unumrec_CreateRan0:   s <= 0");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Ran0_param));
   state = util_Malloc (sizeof (Ran0_state));

   strncpy (name, "unumrec_CreateRan0:", LEN);
   addstr_Long (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S1 = s;
   param->Norm = 1.0 / M1;

   gen->GetBits = &Ran0_Bits;
   gen->GetU01  = &Ran0_U01;
   gen->Write   = &WrRan0;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

static double Ran1_U01 (void *vpar, void *vsta)
{
   Ran0_param *param = vpar;
   Ran1_state *state = vsta;
   long k;
   int j;

   k = state->S1 / q1;
   state->S1 = A1 * (state->S1 - k * q1) - k * r1;
   if (state->S1 < 0)
      state->S1 += M1;

   j = 1 + (int) (state->z / NDiv);
   state->z = state->Tab[j];
   state->Tab[j] = state->S1;

   /* This is slightly different from Press and Teukolsky: they use RNMX. */
   return (state->z * param->Norm);
}

static unsigned long Ran1_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Ran1_U01 (vpar, vsta));
}

static void WrRan1 (void *vsta)
{
   Ran1_state *state = vsta;
   int i;
   if (unif01_WrLongStateFlag) {
      printf (" S = %1ld\n\n", state->S1);
      for (i = 1; i <= N_TAB; i++) {
         printf ("  Tab [%2d] = %12ld\n", i, state->Tab[i]);
      }
   } else
      unif01_WrLongStateDef ();
}

unif01_Gen * unumrec_CreateRan1 (long s)
{
   unif01_Gen *gen;
   Ran0_param *param;
   Ran1_state *state;
   size_t leng;
   char name[LEN + 1];
   int i;
   long k;

   util_Assert (s > 0, "unumrec_CreateRan1:   s <= 0");
   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Ran0_param));
   state = util_Malloc (sizeof (Ran1_state));

   strncpy (name, "unumrec_CreateRan1:", LEN);
   addstr_Long (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->Norm = 1.0 / M1;
   state->S1 = s;

   for (i = N_TAB8; i >= 1; i--) {
      k = state->S1 / q1;
      state->S1 = A1 * (state->S1 - k * q1) - k * r1;
      if (state->S1 < 0)
         state->S1 += M1;
      if (i <= N_TAB)
         state->Tab[i] = state->S1;
   }
   state->z = state->Tab[1];

   gen->GetBits = &Ran1_Bits;
   gen->GetU01  = &Ran1_U01;
   gen->Write   = &WrRan1;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

#undef  M1
#undef  A1
#undef  q1
#undef  r1
#undef  NDiv

#define M1  2147483563
#define M2  2147483399
#define A1  40014
#define A2  40692
#define q1  53668
#define q2  52774
#define r1  12211
#define r2  3791
#define M1m1  2147483562
#define NDiv  67108862        /* = 1 + M1m1 / N_TAB */


static double Ran2_U01 (void *vpar, void *vsta)
{
   Ran0_param *param = vpar;
   Ran2_state *state = vsta;
   long k;
   int j;

   k = state->S1 / q1;
   state->S1 = A1 * (state->S1 - k * q1) - k * r1;
   if (state->S1 < 0)
      state->S1 += M1;

   k = state->S2 / q2;
   state->S2 = A2 * (state->S2 - k * q2) - k * r2;
   if (state->S2 < 0)
      state->S2 += M2;

   j = 1 + (int) (state->z / NDiv);
   state->z = state->Tab[j] - state->S2;
   state->Tab[j] = state->S1;
   if (state->z < 1)
      state->z += M1m1;

   /* This is slightly different from Press and Teukolsky: they use RNMX. */
   return (state->z * param->Norm);
}

static unsigned long Ran2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Ran2_U01 (vpar, vsta));
}

static void WrRan2 (void *vsta)
{
   Ran2_state *state = vsta;
   int i;
   if (unif01_WrLongStateFlag) {
      printf (" S1 = %1ld,   S2 = %1ld\n\n", state->S1, state->S2);
      for (i = 1; i <= N_TAB; i++) {
         printf ("  Tab [%2d] = %12ld\n", i, state->Tab[i]);
      }
   } else
      unif01_WrLongStateDef ();
}

unif01_Gen * unumrec_CreateRan2 (long s)
{
   unif01_Gen *gen;
   Ran0_param *param;
   Ran2_state *state;
   size_t leng;
   char name[LEN + 1];
   int i;
   long k;

   util_Assert (s > 0, "unumrec_CreateRan2:   s <= 0");
   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Ran0_param));
   state = util_Malloc (sizeof (Ran2_state));

   strncpy (name, "unumrec_CreateRan2:", LEN);
   addstr_Long (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->Norm = 1.0 / M1;
   state->S1 = s;
   state->S2 = s;

   for (i = N_TAB8; i >= 1; i--) {
      k = state->S1 / q1;
      state->S1 = A1 * (state->S1 - k * q1) - k * r1;
      if (state->S1 < 0)
         state->S1 += M1;
      if (i <= N_TAB)
         state->Tab[i] = state->S1;
   }
   state->z = state->Tab[1];

   gen->GetBits = &Ran2_Bits;
   gen->GetU01  = &Ran2_U01;
   gen->Write   = &WrRan2;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

void unumrec_DeleteGen (unif01_Gen *gen)
{
   unif01_DeleteGen (gen);
}
