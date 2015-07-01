/*************************************************************************\
 *
 * Package:        TestU01
 * File:           utezuka.c
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

#include "utezuka.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>




/*============================== constants ==============================*/

#define LEN  250                  /* Max length of strings */

#define DeuxExp28      268435456
#define DeuxExp29      536870912
#define DeuxExp31      2147483648U
#define DeuxExp32      4294967296.0

/* For the TezLec91 generator */
#define  Q1  13
#define  Q2  2
#define  S1  12
#define  S2  17
#define  K1mS1  19
#define  K2mS2  12
#define  K1mK2  2
#define  M1  2147483647           /* (2^31) - 1 */
#define  M2  536870911            /* (2^29) - 1 */





/*================================ Types ================================*/

typedef struct {
   double Norm;
} Tezuka_param;

typedef Tezuka_param TezLec91_param;
typedef Tezuka_param Tez95_param;
typedef Tezuka_param TezMRG95_param;

typedef struct {
   unsigned int X1, X2;
} TezLec91_state;

typedef struct {
   unsigned int X1, X2, X3;
} Tez95_state;

typedef struct {
   int j1, k1, j2, k2;
   unsigned int X1[5], X2[7];
} TezMRG95_state;




/*============================== Functions ==============================*/

static void WrTezLec91 (void *vsta)
{
   TezLec91_state *state = vsta;
   printf (" s1 = %1u,   s2 = %1u\n", state->X1, state->X2);
}

static unsigned long TezLec91_Bits (void *junk, void *vsta)
{
   TezLec91_state *state = vsta;
   unsigned int A, B;

   B = (M1 & (state->X1 ^ (state->X1 << Q1))) >> K1mS1;
   A = state->X1 << S1;
   state->X1 = M1 & (A ^ B);

   B = (M2 & (state->X2 ^ (state->X2 << Q2))) >> K2mS2;
   A = state->X2 << S2;
   state->X2 = M2 & (A ^ B);

   return (state->X1 ^ (state->X2 << K1mK2)) << 1;
}

static double TezLec91_U01 (void *vpar, void *vsta)
{
   TezLec91_param *param = vpar;
   return param->Norm * TezLec91_Bits (vpar, vsta);
}


unif01_Gen * utezuka_CreateTezLec91 (unsigned int Y1, unsigned int Y2)
{
   unif01_Gen *gen;
   TezLec91_param *param;
   TezLec91_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert (Y1 < DeuxExp31, "utezuka_CreateTezLec91:   Y1 >= 2^31");
   util_Assert (Y2 < DeuxExp29, "utezuka_CreateTezLec91:   Y2 >= 2^29");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (TezLec91_param));
   state = util_Malloc (sizeof (TezLec91_state));

   strncpy (name, "utezuka_CreateTezLec91:", (size_t) LEN);
   addstr_Uint (name, "   Y1 = ", Y1);
   addstr_Uint (name, ",   Y2 = ", Y2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->X1 = Y1;
   state->X2 = Y2;
   param->Norm = 1.0 / DeuxExp32;

   gen->GetBits = &TezLec91_Bits;
   gen->GetU01  = &TezLec91_U01;
   gen->Write   = &WrTezLec91;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

static void WrTez95 (void *vsta)
{
   Tez95_state *state = vsta;
   printf (" s1 = %1u,   s2 = %1u,   s3 =  %1u\n",
            state->X1, state->X2, state->X3);
}

static unsigned long Tez95_Bits (void *junk, void *vsta)
{
   Tez95_state *state = vsta;
   unsigned int B;

   B = ((state->X1 << 9) ^ state->X1) << 4;
   state->X1 = (state->X1 << 13) ^ (B >> 19);

   B = ((state->X2 << 2) ^ state->X2) << 3;
   state->X2 = (state->X2 << 20) ^ (B >> 12);

   B = ((state->X3 << 6) ^ state->X3) << 1;
   state->X3 = (state->X3 << 17) ^ (B >> 15);

   return state->X1 ^ state->X2 ^ state->X3;
}

static double Tez95_U01 (void *vpar, void *vsta)
{
   Tez95_param *param = vpar;
   return param->Norm * Tez95_Bits (vpar, vsta);
}


unif01_Gen *utezuka_CreateTez95 (unsigned int Y1, unsigned int Y2,
   unsigned int Y3)
{
   unif01_Gen *gen;
   Tez95_param *param;
   Tez95_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned int B;

   util_Assert (Y1 < DeuxExp28, "utezuka_CreateTez95:   Y1 >= 2^28");
   util_Assert (Y2 < DeuxExp29, "utezuka_CreateTez95:   Y2 >= 2^29");
   util_Assert (Y3 < DeuxExp31, "utezuka_CreateTez95:   Y3 >= 2^31");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Tez95_param));
   state = util_Malloc (sizeof (Tez95_state));

   strncpy (name, "utezuka_CreateTez95:", (size_t) LEN);
   addstr_Uint (name, "   Y1 = ", Y1);
   addstr_Uint (name, ",   Y2 = ", Y2);
   addstr_Uint (name, ",   Y3 = ", Y3);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   B = ((Y1 << 9) ^ Y1) << 4;
   state->X1 = (Y1 << 4) ^ (B >> 28);

   B = ((Y2 << 2) ^ Y2) << 3;
   state->X2 = (Y2 << 3) ^ (B >> 29);

   B = ((Y3 << 6) ^ Y3) << 1;
   state->X3 = (Y3 << 1) ^ (B >> 31);

   param->Norm = 1.0 / DeuxExp32;

   gen->GetBits = &Tez95_Bits;
   gen->GetU01  = &Tez95_U01;
   gen->Write   = &WrTez95;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

static void WrTezMRG95 (void *vsta)
{
   TezMRG95_state *state = vsta;
   int k;
   if (unif01_WrLongStateFlag) {
      printf (" S1 = (");
      for (k = 0; k < 5; k++)
         printf ("%12u ", state->X1[k]);
      printf (" )\n\nS2 = (");
      for (k = 0; k < 7; k++) {
         printf ("%12u ", state->X2[k]);
         if (k == 4)
            printf ("\n      ");
      }
      printf (" )\n\n");
   } else
      unif01_WrLongStateDef ();
}

static unsigned long TezMRG95_Bits (void *junk, void *vsta)
{
   TezMRG95_state *state = vsta;
   unsigned int *X1 = state->X1;
   unsigned int *X2 = state->X2;
   unsigned int b1, b2;

   b1 = ((X1[state->k1] << 3) ^ X1[state->k1]) << 1;
   b2 = ((X1[state->j1] << 3) ^ X1[state->j1]) << 1;
   state->X1[state->k1] = (X1[state->k1] << 23) ^ (b1 >> 9) ^
                          (X1[state->j1] << 5) ^ (b2 >> 27);

   b1 = ((X2[state->k2] << 2) ^ X2[state->k2]) << 3;
   b2 = ((X2[state->j2] << 2) ^ X2[state->j2]) << 3;
   state->X2[state->k2] = (X2[state->k2] << 19) ^ (b1 >> 13) ^
                          (X2[state->j2] << 16) ^ (b2 >> 16);

   --state->j1;
   if (state->j1 < 0)
      state->j1 = 4;

   --state->k1;
   if (state->k1 < 0)
      state->k1 = 4;

   --state->j2;
   if (state->j2 < 0)
      state->j2 = 6;

   --state->k2;
   if (state->k2 < 0)
      state->k2 = 6;

   return state->X1[state->k1] ^ state->X2[state->k2];
}

static double TezMRG95_U01 (void *vpar, void *vsta)
{
   TezMRG95_param *param = vpar;
   return param->Norm * TezMRG95_Bits (vpar, vsta);
}


unif01_Gen *utezuka_CreateTezMRG95 (unsigned int Y1[5], unsigned int Y2[7])
{
   unif01_Gen *gen;
   TezMRG95_param *param;
   TezMRG95_state *state;
   size_t leng;
   char name[LEN + 1];
   int k;
   unsigned int b;

   for (k = 0; k < 5; k++)
      util_Assert (Y1[k] < DeuxExp31,
         "utezuka_CreateTezMRG95:   Y1[k] >= 2^31");
   for (k = 0; k < 7; k++)
      util_Assert (Y2[k] < DeuxExp29,
         "utezuka_CreateTezMRG95:   Y2[k] >= 2^29");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (TezMRG95_param));
   state = util_Malloc (sizeof (TezMRG95_state));

   strncpy (name, "utezuka_CreateTezMRG95:", (size_t) LEN);
   addstr_ArrayUint (name, "   Y1 = ", 5, Y1);
   addstr_ArrayUint (name, ",   Y2 = ", 7, Y2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   for (k = 0; k < 5; k++) {
      b = ((Y1[k] << 3) ^ Y1[k]) << 1;
      state->X1[k] = (Y1[k] << 1) ^ (b >> 31);
   }

   for (k = 0; k < 7; k++) {
      b = ((Y2[k] << 2) ^ Y2[k]) << 3;
      state->X2[k] = (Y2[k] << 3) ^ (b >> 29);
   }
   state->j1 = 1;
   state->k1 = 4;
   state->j2 = 4;
   state->k2 = 6;

   param->Norm = 1.0 / DeuxExp32;

   gen->GetBits = &TezMRG95_Bits;
   gen->GetU01  = &TezMRG95_U01;
   gen->Write   = &WrTezMRG95;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

void utezuka_DeleteGen (unif01_Gen *gen)
{
   unif01_DeleteGen (gen);
}
