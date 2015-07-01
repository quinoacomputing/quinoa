/*************************************************************************\
 *
 * Package:        TestU01
 * File:           uweyl.c
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
#include "uweyl.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>


#define LEN 200                   /* Max length of strings */



/*================================ Types ================================*/

typedef struct {
   unsigned long n;
} Weyl_state;

typedef struct {
   double Alpha;
} Weyl_param;

typedef struct {
   double Alpha;
   long M;
} SNWeyl_param;


/*========================= Common functions ============================*/

static void WrWeyl (void *vsta)
{
   Weyl_state *state = vsta;
   printf (" n = %1lu\n", state->n);
}

static unif01_Gen *CreateWeyl_0 (double Alpha, long n0, char name[])
{
   unif01_Gen *gen;
   Weyl_param *param;
   Weyl_state *state;
   size_t leng;

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Weyl_param));
   state = util_Malloc (sizeof (Weyl_state));
   param->Alpha = Alpha;
   state->n = n0;

   addstr_Double (name, "  Alpha = ", Alpha);
   addstr_Long (name, ",   n0 = ", n0);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->Write = &WrWeyl;
   gen->param = param;
   gen->state = state;
   return gen;
}


/*=======================================================================*/

static double Weyl_U01 (void *vpar, void *vsta)
{
   Weyl_param *param = vpar;
   Weyl_state *state = vsta;
   double X;

   state->n++;
   X = state->n * param->Alpha;
   return (X - (long) X);
}

static unsigned long Weyl_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Weyl_U01 (vpar, vsta));
}

unif01_Gen *uweyl_CreateWeyl (double Alpha, long n0)
{
   unif01_Gen *gen;
   char name[LEN + 1];

   util_Assert ((Alpha > 0.0), "uweyl_CreateWeyl:   Alpha <= 0");
   util_Assert ((Alpha < 1.0), "uweyl_CreateWeyl:   Alpha >= 1");

   strncpy (name, "uweyl_CreateWeyl: ", (size_t) LEN);
   gen = CreateWeyl_0 (Alpha, n0, name);
   gen->GetU01 = &Weyl_U01;
   gen->GetBits = &Weyl_Bits;
   return gen;
}


/*=======================================================================*/

static double NWeyl_U01 (void *vpar, void *vsta)
{
   Weyl_param *param = vpar;
   Weyl_state *state = vsta;
   double X;

   state->n++;
   X = state->n * param->Alpha;
   X -= (long) X;
   X *= state->n;
   return (X - (long) X);
}

static unsigned long NWeyl_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * NWeyl_U01 (vpar, vsta));
}

unif01_Gen *uweyl_CreateNWeyl (double Alpha, long n0)
{
   unif01_Gen *gen;
   char name[LEN + 1];

   util_Assert ((Alpha > 0.0), "uweyl_CreateNWeyl:   Alpha <= 0");
   util_Assert ((Alpha < 1.0), "uweyl_CreateNWeyl:   Alpha >= 1");

   strncpy (name, "uweyl_CreateNWeyl (nested): ", (size_t) LEN);
   gen = CreateWeyl_0 (Alpha, n0, name);
   gen->GetU01 = &NWeyl_U01;
   gen->GetBits = &NWeyl_Bits;
   return gen;
}


/*=======================================================================*/

static double SNWeyl_U01 (void *vpar, void *vsta)
{
   /* This function gives very different results on the same machine
      depending on the level of optimization in the compilation.  */
   SNWeyl_param *param = vpar;
   Weyl_state *state = vsta;
   double X, Z;

   state->n++;
   X = state->n * param->Alpha;
   X -= (long) X;
   X *= state->n;
   X -= (long) X;
   X = X * param->M + 0.5;
   Z = X * param->Alpha;
   Z -= (long) Z;
   Z *= X;
   Z -= (long) Z;
   return Z;
}

static unsigned long SNWeyl_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * SNWeyl_U01 (vpar, vsta));
}

unif01_Gen *uweyl_CreateSNWeyl (long M, double Alpha, long n0)
{
   unif01_Gen *gen;
   SNWeyl_param *param;
   Weyl_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert ((Alpha > 0.0), "uweyl_CreateSNWeyl:   Alpha <= 0");
   util_Assert ((Alpha < 1.0), "uweyl_CreateSNWeyl:   Alpha >= 1");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (SNWeyl_param));
   state = util_Malloc (sizeof (Weyl_state));
   param->Alpha = Alpha;
   param->M = M;
   if (n0 < 0)
      n0 = -n0;
   state->n = n0;

   strncpy (name, "uweyl_CreateSNWeyl (shuffled nested):", (size_t) LEN);
   addstr_Long (name, "   M = ", M);
   addstr_Double (name, ",   Alpha = ", Alpha);
   addstr_Long (name, ",   n0 = ", n0);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &SNWeyl_Bits;
   gen->GetU01  = &SNWeyl_U01;
   gen->Write   = &WrWeyl;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/*=======================================================================*/

void uweyl_DeleteGen (unif01_Gen *gen)
{
   unif01_DeleteGen (gen);
}
