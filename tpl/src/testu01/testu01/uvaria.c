/*************************************************************************\
 *
 * Package:        TestU01
 * File:           uvaria.c
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

#include "uvaria.h"
#include "unif01.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>




static int co = 0;                 /* Counter for Ranrot */




/*============================= Constants =================================*/

#define INV31 4.656612875245796E-10  /* 1 / (2^31 - 1) */

#define LEN  200                  /* Max length of strings */





/*================================ Types ================================*/

typedef struct {
   double *X;
   int K;
} ACORN_state;

/*------------------------------------*/

typedef struct {
   long B, a0, a1;
   double BK[33];
   int K;
} Tindo_param;

typedef struct {
   long C[33];
   int Current, L;
} Tindo_state;

/*------------------------------------*/

typedef struct {
   long S;
   unsigned long V;
} CSD_state;

/*------------------------------------*/

typedef struct {
   double A1, B1, A2, B2;
} Rey97_param;

typedef struct {
   unsigned long n;
} Rey97_state;




    

/*============================== Functions ==============================*/

static double ACORN_U01 (void *junk, void *vsta)
{
   ACORN_state *state = vsta;
   int j;
   double temp;

   for (j = 0; j < state->K; j++) {
      temp = state->X[j + 1] + state->X[j];
      if (temp >= 1.0)
         state->X[j + 1] = temp - 1.0;
      else
         state->X[j + 1] = temp;
   }
   return state->X[state->K];
}

/*-----------------------------------------------------------------------*/

static unsigned long ACORN_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (ACORN_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrACORN (void *vsta)
{
   ACORN_state *state = vsta;
   int j;
   if (unif01_WrLongStateFlag) {
      printf (" S = {\n");
      for (j = 0; j < state->K; j++) {
         printf (" %22.16f", state->X[j]);
         if (j < state->K)
            printf (",");
         if (((j + 1) % 3) == 0)
            printf ("\n");
      };
      printf ("\n     }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

unif01_Gen * uvaria_CreateACORN (int k, double S[])
{
   unif01_Gen *gen;
   ACORN_state *state;
   size_t leng;
   char name[LEN + 1];
   int j;

   if (k < 1) {
      util_Error ("uvaria_CreateACORN:   k < 1");
   }
   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (ACORN_state));

   state->X = util_Calloc ((size_t) k + 1, sizeof (double));
   for (j = 0; j < k; j++)
      state->X[j] = S[j];
   state->X[k] = 0.1234567;    /* Arbitrary initial value */

   strcpy (name, "uvaria_CreateACORN:");
   addstr_Int (name, "   k = ", k);
   addstr_ArrayDouble (name, ",   S = ", k, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->K = k;

   gen->GetBits = &ACORN_Bits;
   gen->GetU01  = &ACORN_U01;
   gen->Write   = &WrACORN;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}

/*------------------------------------------------------------------------*/

void uvaria_DeleteACORN (unif01_Gen * gen)
{
   ACORN_state *state;

   if (NULL == gen)
      return;
   state = gen->state;
   util_Free (state->X);
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/**************************************************************************/

static double Tindo_U01 (void *vpar, void *vsta)
{
   Tindo_param *param = vpar;
   Tindo_state *state = vsta;
   long T[33];
   double Sum;
   int j, i;

   Sum = 0.0;
   for (j = 1; j <= param->K; j++) {
      if (state->Current == 0) {
         /* Refresh table C */
         T[1] = (param->a1 * state->C[state->L] + param->a0 * state->C[1] + 1)
               % param->B;
         for (i = 2; i <= state->L; i++)
            T[i] = (param->a1 * state->C[i - 1] + param->a0 * state->C[i] + 1)
               % param->B;
         for (i = 1; i <= state->L; i++)
            state->C[i] = T[i];
      }
      Sum += param->BK[j] * state->C[state->Current + 1];
      state->Current = (state->Current + 1) % state->L;
   }
   while (Sum > 1.0)
      Sum -= 1.0;
   return Sum;
}

/*-----------------------------------------------------------------------*/

static unsigned long Tindo_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (Tindo_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrTindo (void *vsta)
{
   Tindo_state *state = vsta;
   int i;
   if (unif01_WrLongStateFlag) {
      printf (" C = {\n");
      for (i = 1; i <= state->L; i++) {
         printf ("    %10ld", state->C[i]);
         if (i < state->L)
            printf (",");
         if ((i % 3) == 0)
            printf ("\n");
      };
      printf ("\n     }");
      printf ("\n Current = %1d\n", state->Current);
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

unif01_Gen * uvaria_CreateTindo (long b, long Delta, int l, int k)
/*
 * Assumes that ( a0 < b ) , ( a1 < b ) and ( b < 2^15 = 32768 )
 */
{
   unif01_Gen *gen;
   Tindo_param *param;
   Tindo_state *state;
   size_t leng;
   char name[LEN + 1];
   int i;

   util_Assert (k > 0, "uvaria_CreateTindo:   must have k > 0");
   util_Assert (k <= 32, "uvaria_CreateTindo:   must have k <= 32");
   util_Assert (l > 0, "uvaria_CreateTindo:   must have l > 0");
   util_Assert (l <= 32, "uvaria_CreateTindo:   must have l <= 32");
   util_Assert (b < 32768, "uvaria_CreateTindo:   must have b < 2^15");
   util_Assert (Delta > 0, "uvaria_CreateTindo:   must have Delta > 0");
   util_Assert (Delta < b - 1,
      "uvaria_CreateTindo:   must have Delta < b - 1");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Tindo_param));
   state = util_Malloc (sizeof (Tindo_state));

   strcpy (name, "uvaria_CreateTindo:");
   addstr_Long (name, "   b = ", b);
   addstr_Long (name, ",   Delta = ", Delta);
   addstr_Int (name, ",   s = ", l);
   addstr_Int (name, ",   k = ", k);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->Current = 0;
   state->L = l;
   for (i = 1; i <= l; i++)
      state->C[i] = i % b;
   param->B = b;
   param->K = k;
   param->a1 = Delta + 1;
   param->a0 = b - Delta;
   param->BK[1] = 1.0 / b;
   for (i = 2; i <= k; i++)
      param->BK[i] = param->BK[i - 1] * param->BK[1];

   gen->GetBits = &Tindo_Bits;
   gen->GetU01  = &Tindo_U01;
   gen->Write   = &WrTindo;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

static double CSD_U01 (void *junk, void *vsta)
{
   CSD_state *state = vsta;
   double u;
   unsigned long d1, d2, temp, c1, c4;
   long v;

   v = state->S / 127773;
   state->S = 16807 * (state->S - v * 127773) - v * 2836;
   if (state->S < 0)
      state->S += 2147483647;
   u = state->S * INV31;

   d1 = 10.0 * u;
   d2 = ((unsigned long) (100.0 * u)) - 10 * d1;

   state->V = (state->V + d1) % 10000;
   temp = ((state->V * state->V) % 10000) * state->V;
   c1 = temp % 10;
   temp /= 1000;
   c4 = temp % 10;
   state->V = c1 * 1000 + d1 * 100 + c4 * 10 + d2;
   return state->V * 0.0001;
}

/*-----------------------------------------------------------------------*/

static unsigned long CSD_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (CSD_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrCSD (void *vsta)
{
   CSD_state *state = vsta;
   printf (" V = %1lu,     S = %1ld\n", state->V, state->S);
}

/*-----------------------------------------------------------------------*/

unif01_Gen * uvaria_CreateCSD (long v, long s)
{
   unif01_Gen *gen;
   CSD_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert (v >= 0, "uvaria_CreateCSD:   must have v >= 0");
   util_Assert (v <= 9999, "uvaria_CreateCSD:   must have v <= 9999");
   util_Assert (s > 0, "uvaria_CreateCSD:   must have s > 0");
   util_Assert (s < 2147483647,
      "uvaria_CreateCSD:   must have s < 2^31 - 1");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (CSD_state));

   strcpy (name, "uvaria_CreateCSD:");
   addstr_Long (name, "   v = ", v);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->V = v;
   state->S = s;

   gen->GetBits = &CSD_Bits;
   gen->GetU01  = &CSD_U01;
   gen->Write   = &WrCSD;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

/***********************************************************************/
/*
   This code is from Agner Fog:  http://www.agner.org
*/
/************************* RANROTB.C ******************** AgF 1999-03-03 *
*  Random Number generator 'RANROT' type B                               *
*                                                                        *
*  This is a lagged-Fibonacci type of random number generator with       *
*  rotation of bits.  The algorithm is:                                  *
*  X[n] = ((X[n-j] rotl r1) + (X[n-k] rotl r2)) modulo 2^b               *
*                                                                        *
*  The last k values of X are stored in a circular buffer named          *
*  randbuffer.                                                           *
*                                                                        *
*  This version works with any integer size: 16, 32, 64 bits etc.        *
*  The integers must be unsigned. The resolution depends on the integer  *
*  size.                                                                 *
*                                                                        *
*  Note that the function RanrotBInit must be called before the first    *
*  call to RanrotB or iRanrotB                                           *
*                                                                        *
*  The theory of the RANROT type of generators is described at           *
*  www.agner.org/random/ranrot.htm                                       *
*                                                                        *
*************************************************************************/

/* names of header files may differ
#include <stdio.h>
#include <conio.h>
#include <mem.h>
#include <math.h>
#include <time.h>
*/

/* define desired integer type */
typedef  unsigned int my_uint;

/* If your system doesn't have a rotate function for the chosen integer type
   then define it thus:                                                 */
static my_uint rotl (my_uint x, my_uint r) {
  return (x << r) | (x >> (sizeof(x)*8-r));}

/* define parameters (R1 and R2 must be smaller than the integer size): */
#define KK  17
#define JJ  10
#define R1   5
#define R2   3

/* global variables */
static my_uint randbuffer[KK];  /* history buffer */
static int r_p1, r_p2;          /* indexes into history buffer */
static float scale;             /* 2^(- integer size) */


/* returns a random number between 0 and 1 */
static double RanrotB(void) {
  my_uint x;
  /* generate next random number */
  x = randbuffer[r_p1] = rotl(randbuffer[r_p2], R1) + rotl(randbuffer[r_p1], R2);
  /* rotate list pointers */
  if (--r_p1 < 0) r_p1 = KK - 1;
  if (--r_p2 < 0) r_p2 = KK - 1;
  /* conversion to float */
  return x * scale;
}

/* get integer random number in interval from min to max */
static int iRanrotB(int min, int max) {
  int i, r;
  i = max - min + 1;
  r = i * RanrotB();
  if (r >= i) r = i-1;
  return min + r;
}

/* this function initializes the random number generator.      */
/* Must be called before the first call to RanrotB or iRanrotB */
static void RanrotBInit (my_uint seed) {
  int i;

  /* put semi-random numbers into the buffer */
  for (i=0; i<KK; i++) {
    randbuffer[i] = seed;
    seed = rotl(seed,5) + 97;}

  /* initialize pointers to circular buffer */
  r_p1 = 0;  r_p2 = JJ;

  /* randomize */
  for (i = 0;  i < 300;  i++) RanrotB();

  /* compute 2^(- integer size) */
  scale = ldexp(1., -8* (int) sizeof(my_uint));
}


/*********************** End of code from Agner Fog **********************/


/*------------------------ Our code restarts here -----------------------*/

static void WrRanrotB (void *junk)
{  
}

/*-----------------------------------------------------------------------*/

static double RanrotB_U01 (void *junk1, void *junk2)
{
   return RanrotB ();
}

/*-----------------------------------------------------------------------*/

static unsigned long RanrotB_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (RanrotB_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

unif01_Gen * uvaria_CreateRanrotB (unsigned int s)
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];

   util_Assert (co == 0,
      "uvaria_CreateRanrotB:   only 1 generator at a time can be in use");
   co++;

   RanrotBInit (s);

   gen = util_Malloc (sizeof (unif01_Gen));
   strcpy (name, "uvaria_CreateRanrotB:");
   addstr_Uint (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &RanrotB_Bits;
   gen->GetU01 = &RanrotB_U01;
   gen->Write = &WrRanrotB;
   gen->param = NULL;
   gen->state = NULL;
   return gen;
}

/*-----------------------------------------------------------------------*/

void uvaria_DeleteRanrotB (unif01_Gen *gen)
{
   unif01_DeleteGen (gen);
   co--;
}


/**************************************************************************/

static double Rey97_U01 (void *vpar, void *vsta)
{
   Rey97_state *state = vsta;
   Rey97_param *param = vpar;
   double U, Z;

   Z = param->A1 * sin (param->B1 * state->n);
   Z = modf (Z, &U);        /* U is unused  */
   if (Z < 0.0)
      Z += 1.0;
   ++state->n;

   U = (param->A2 + Z) * sin (param->B2 * Z);
   U = modf (U, &Z);        /* Z is unused  */
   if (U < 0.0)
      return U + 1.0;
   else
      return U;
}

/*-----------------------------------------------------------------------*/

static unsigned long Rey97_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (Rey97_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrRey97 (void *vsta)
{
   Rey97_state *state = vsta;
   printf (" n = %1lu\n", state->n);
}

/*-----------------------------------------------------------------------*/

unif01_Gen * uvaria_CreateRey97 (double a1, double a2, double b2, long n0)
{
   unif01_Gen *gen;
   Rey97_state *state;
   Rey97_param *param;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (Rey97_state));
   param = util_Malloc (sizeof (Rey97_param));

   strcpy (name, "uvaria_CreateRey97:");
   addstr_Double (name, "   a1 = ", a1);
   addstr_Double (name, ",   a2 = ", a2);
   addstr_Double (name, ",   b2 = ", b2);
   addstr_Long (name, ",   n0 = ", n0);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->n = n0;
   param->A1 = a1;
   param->A2 = a2;
   param->B2 = b2;
   param->B1 = (sqrt (5.0) - 1) * num_Pi / 2;

   gen->GetBits = &Rey97_Bits;
   gen->GetU01  = &Rey97_U01;
   gen->Write   = &WrRey97;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

void uvaria_DeleteGen (unif01_Gen *gen)
{
   unif01_DeleteGen (gen);
}
