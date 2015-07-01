/*************************************************************************\
 *
 * Package:        TestU01
 * File:           umarsa.c
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

#include "umarsa.h"
#include "unif01.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>





/*============================== Constants ==============================*/

#define LEN  200                  /* Max length of strings */

#define MASK8  0xff               /* Mask of 8 bits */
#define MASK16 0xffff             /* Mask of 16 bits */
#define MASK31 0x7fffffffU        /* Mask of 31 bits */
#define MASK32 0xffffffffUL       /* Mask of 32 bits */
#define MA     4294967291UL       /* 2^32 - 5 */

#define INV32L  2.32830643653869628906250E-10L   /* 1 / 2^32 */
#define INV32   2.3283064365386963E-10           /* 1 / 2^32 */
#define INV32m1 2.3283064370807973E-10           /* 1 / (2^32 - 1) */ 



/*================================ Types ================================*/

typedef struct {
   unsigned int S1, S2, S3, S4, carry;
} KISS_state;


typedef struct {
   unsigned int r, s;            /* Lags */
   unsigned int C;              /* Carry */
   unsigned int X[43];
   unsigned int W;              /* State of the Weyl generator */
} Marsa90a_state;


typedef struct {
   double SeqD, SeqM;
} RANMAR_param;


typedef struct {
   double X[98];
   unsigned int r, s;
   double Seq;
} RANMAR_state;


#ifdef USE_LONGLONG
typedef struct {
   ulonglong S1, S2, S3, S4;    /* State */
   ulonglong C;                 /* Carry */
} Mother0_state;
#endif

typedef struct {
   unsigned long X1, X2, Y1;
} Combo_state;


typedef struct {
   double X1, X2, X3;
} ECG_state;

typedef struct {
   double Norm;
} ECG_param;


typedef struct{
   unsigned int I1, I2;
} MWC97R_state;


typedef struct{
   unsigned int X[256];
   unsigned int r;
   unsigned int b;
} LFIB4_99_state;

#define ULTRA_R 99

typedef struct {
   unsigned long X[1 + ULTRA_R];
   int r, s;
   unsigned long I;
} ULTRA_state;


#ifdef USE_LONGLONG
typedef struct{
   ulonglong x, y;
} SupDup64_state;

typedef struct{
   ulonglong a, c;
   unsigned int s1, s2, s3;
} SupDup64_param;
#endif


typedef MWC97R_state FIB_99_state;
typedef MWC97R_state SupDup73_state;
typedef MWC97R_state SupDup96_state;
typedef LFIB4_99_state SWB_99_state;





/*============================== Functions ==============================*/



static unsigned long KISS93_Bits (void *junk, void *vsta)
{
   KISS_state *state = vsta;
   unsigned int b;

   state->S1 = 69069 * state->S1 + 23606797;
   b = state->S2 ^ (state->S2 << 17);
   state->S2 = (b >> 15) ^ b;
   b = ((state->S3 << 18) ^ state->S3) & MASK31;
   state->S3 = (b >> 13) ^ b;
   return state->S1 + state->S2 + state->S3;
}

/*-----------------------------------------------------------------------*/

static double KISS93_U01 (void *vpar, void *vsta)
{
   return KISS93_Bits (vpar, vsta) * INV32;
}

/*-----------------------------------------------------------------------*/

static void WrKISS93 (void *vsta)
{
   KISS_state *state = vsta;
   printf (" x0 = %1u,    y0 = %1u,    z0 = %1u\n\n",
      state->S1, state->S2, state->S3);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateKISS93 (unsigned int s1, unsigned int s2,
   unsigned int s3)
{
   unif01_Gen *gen;
   KISS_state *state;
   size_t leng;
   char name[LEN + 1];

   if (s3 > 2147483647)
      util_Error ("umarsa_CreateKISS93:   s3 >= 2^31");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (KISS_state));

   strcpy (name, "umarsa_CreateKISS93:");
   addstr_Uint (name, "   x0 = ", s1);
   addstr_Uint (name, ",   y0 = ", s2);
   addstr_Uint (name, ",   z0 = ", s3);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S1 = s1;
   state->S2 = s2;
   state->S3 = s3;

   gen->GetBits = &KISS93_Bits;
   gen->GetU01  = &KISS93_U01;
   gen->Write   = &WrKISS93;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

static unsigned long KISS96_Bits (void *junk, void *vsta)
{
   KISS_state *state = vsta;
   unsigned int y, k;

   state->S1 = 69069 * state->S1 + 1;
   y = state->S2 ^ (state->S2 << 13);
   y ^= y >> 17;
   state->S2 = y ^ (y << 5);
   k = (state->S3 >> 2) + (state->S4 >> 3) + (state->carry >> 2);
   y = state->S4 + state->S4 + state->S3 + state->carry;
   state->S3 = state->S4;
   state->S4 = y;
   state->carry = k >> 30;
   return state->S1 + state->S2 + state->S4;
}

/*-----------------------------------------------------------------------*/

static double KISS96_U01 (void *vpar, void *vsta)
{
   return KISS96_Bits (vpar, vsta) * INV32;
}

/*-----------------------------------------------------------------------*/

static void WrKISS96 (void *vsta)
{
   KISS_state *state = vsta;
   printf (" x = %u,    y = %u,    z1 = %u,    z2 = %u\n\n",
      state->S1, state->S2, state->S3, state->S4);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateKISS96 (unsigned int x, unsigned int y,
   unsigned int z1, unsigned int z2)
{
   unif01_Gen *gen;
   KISS_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (KISS_state));

   strcpy (name, "umarsa_CreateKISS96:");
   addstr_Uint (name, "   x = ", x);
   addstr_Uint (name, ",   y = ", y);
   addstr_Uint (name, ",   z1 = ", z1);
   addstr_Uint (name, ",   z2 = ", z2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S1 = x;
   state->S2 = y;
   state->S3 = z1;
   state->S4 = z2;
   state->carry = 0;

   gen->GetBits = &KISS96_Bits;
   gen->GetU01  = &KISS96_U01;
   gen->Write   = &WrKISS96;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

static unsigned long KISS99_Bits (void *junk, void *vsta)
{
   KISS_state *state = vsta;
   unsigned int b;

   state->S1 = 69069 * state->S1 + 1234567;
   b = state->S2 ^ (state->S2 << 17);
   b ^= (b >> 13);
   state->S2 = b ^ (b << 5);
   state->S3 = 36969 * (state->S3 & MASK16) + (state->S3 >> 16);
   state->S4 = 18000 * (state->S4 & MASK16) + (state->S4 >> 16);
   b = (state->S3 << 16) + state->S4;
   return state->S2 + (state->S1 ^ b);
}

/*-----------------------------------------------------------------------*/

static double KISS99_U01 (void *vpar, void *vsta)
{
   return KISS99_Bits (vpar, vsta) * INV32;
}

/*-----------------------------------------------------------------------*/

static void WrKISS99 (void *vsta)
{
   KISS_state *state = vsta;
   printf (" x = %1u,    y = %1u,    I1 = %1u,    I2 = %1u\n\n",
      state->S1, state->S2, state->S3, state->S4);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateKISS99 (unsigned int x, unsigned int y,
   unsigned int I1, unsigned int I2)
{
   unif01_Gen *gen;
   KISS_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (KISS_state));

   strcpy (name, "umarsa_CreateKISS99:");
   addstr_Uint (name, "   x0 = ", x);
   addstr_Uint (name, ",   y0 = ", y);
   addstr_Uint (name, ",   I1 = ", I1);
   addstr_Uint (name, ",   I2 = ", I2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S1 = x;
   state->S2 = y;
   state->S3 = I1;
   state->S4 = I2;

   gen->GetBits = &KISS99_Bits;
   gen->GetU01  = &KISS99_U01;
   gen->Write   = &WrKISS99;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

static unsigned long ULTRA_Bits (void *junk, void *vsta)
{
   ULTRA_state *state = vsta;
   unsigned long t;

   state->X[state->r] = state->X[state->r] * state->X[state->s];
   t = state->X[state->r] &= MASK32;
   if (--state->r < 0)
      state->r = 96;
   if (--state->s < 0)
      state->s = 96;

   state->I = (state->I & MASK16) * 30903 + (state->I >> 16);
   return (t + state->I) & MASK32;
}

/*-----------------------------------------------------------------------*/

static double ULTRA_U01 (void *vpar, void *vsta)
{
   return ULTRA_Bits (vpar, vsta) * INV32;
}

/*-----------------------------------------------------------------------*/

static void WrULTRA (void *junk)
{
   unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/
/*
 * I believe there are errors in that gen, but I have chosen to make it
 * output the same numbers as those from the ULTRA in Diehard. The 
 * Lag Fib is supposed to be of order 99, but the programmed version in
 * Diehard uses only 97 elements of the array.
 */
unif01_Gen *umarsa_CreateULTRA (unsigned int x, unsigned int y,
   unsigned int z, unsigned int w)
{
   unif01_Gen *gen;
   ULTRA_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned int js;
   int i;

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (ULTRA_state));

   strcpy (name, "umarsa_CreateULTRA:");
   addstr_Uint (name, "   s1 = ", x);
   addstr_Uint (name, ",   s2 = ", y);
   addstr_Uint (name, ",   s3 = ", z);
   addstr_Uint (name, ",   s4 = ", w);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->I = x + y + z + w;
   for (i = 0; i < 99; ++i) {
      x = (x & 65535) * 18273 + (x >> 16);
      y = (y & 65535) * 23163 + (y >> 16);
      z = (z & 65535) * 24984 + (z >> 16);
      w = (w & 65535) * 28854 + (w >> 16);
      js = (x << 16) + (y & 65535) + (z << 16) + (w & 65535);
      state->X[i] = js | 1;
   }
   state->r = 98;
   state->s = 32;
   gen->GetBits = &ULTRA_Bits;
   gen->GetU01  = &ULTRA_U01;
   gen->Write   = &WrULTRA;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

static unsigned long SupDup73_Bits (void *junk, void *vsta)
{
   SupDup73_state *state = vsta;
   unsigned int b;

   state->I1 = 69069 * state->I1;
   b = state->I2 ^ (state->I2 >> 15);
   state->I2 = b ^ (b << 17);
   return  state->I2 ^ state->I1;
}

/*-----------------------------------------------------------------------*/

static double SupDup73_U01 (void *vpar, void *vsta)
{
   return SupDup73_Bits (vpar, vsta) * INV32m1;
}

/*-----------------------------------------------------------------------*/

static void WrSupDup73 (void *vsta)
{
   SupDup73_state *state = vsta;
   printf (" x = %1u,    y = %1u\n\n", state->I1, state->I2);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateSupDup73 (unsigned int x, unsigned int y)
{
   unif01_Gen *gen;
   SupDup73_state *state;
   size_t leng;
   char name[LEN + 1];


   util_Warning (!(x & 1), "umarsa_CreateSupDup73:   x reset to odd");
   x |= 1;
   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (SupDup73_state));
   state->I1 = x;
   state->I2 = y;

   strcpy (name, "umarsa_CreateSupDup73:");
   addstr_Uint (name, "   x0 = ", x);
   addstr_Uint (name, ",   y0 = ", y);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &SupDup73_Bits;
   gen->GetU01  = &SupDup73_U01;
   gen->Write   = &WrSupDup73;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

static unsigned long SupDup96ADD_Bits (void *vpar, void *vsta)
{
   SupDup96_state *state = vsta;
   unsigned int *pc = vpar;
   unsigned int b;

   state->I1 = 69069U * state->I1 + *pc;
   b = state->I2 ^ (state->I2 << 13);
   b ^= (b >> 17);
   state->I2 = b ^ (b << 5);
   return state->I2 + state->I1;
}

/*-----------------------------------------------------------------------*/

static double SupDup96ADD_U01 (void *vpar, void *vsta)
{
   return SupDup96ADD_Bits (vpar, vsta) * INV32;
}

/*-----------------------------------------------------------------------*/

static unsigned long SupDup96XOR_Bits (void *vpar, void *vsta)
{
   SupDup96_state *state = vsta;
   unsigned int *pc = vpar;
   unsigned int b;

   state->I1 = 69069U * state->I1 + *pc;
   b = state->I2 ^ (state->I2 << 13);
   b ^= (b >> 17);
   state->I2 = b ^ (b << 5);
   return  state->I2 ^ state->I1;
}

/*-----------------------------------------------------------------------*/

static double SupDup96XOR_U01 (void *vpar, void *vsta)
{
   return SupDup96XOR_Bits (vpar, vsta) * INV32;
}

/*-----------------------------------------------------------------------*/

static void WrSupDup96 (void *vsta)
{
   SupDup96_state *state = vsta;
   printf (" x = %1u,    y = %1u\n\n", state->I1, state->I2);
}

/*-----------------------------------------------------------------------*/

static unif01_Gen *CreateSupDup96 (unsigned int x, unsigned int y,
   unsigned int c, char op)
{
   unif01_Gen *gen;
   SupDup96_state *state;
   unsigned int *p;   
   size_t leng;
   char name[LEN + 1];

   util_Assert (op == '+' || op == 'x',
      "umarsa_CreateSupDup96:   op must be '+' or 'x'");
   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (SupDup96_state));
   p = util_Malloc (sizeof (unsigned int));
   state->I1 = x;
   state->I2 = y;
   *p = c | 1;

   if (op == '+')
      strcpy (name, "umarsa_CreateSupDup96Add:");
   else
      strcpy (name, "umarsa_CreateSupDup96Xor:");
   addstr_Uint (name, "   x0 = ", x);
   addstr_Uint (name, ",   y0 = ", y);
   addstr_Uint (name, ",   c = ", *p);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if (op == 'x') {
      gen->GetBits = &SupDup96XOR_Bits;
      gen->GetU01 = &SupDup96XOR_U01;
   } else {
      gen->GetBits = &SupDup96ADD_Bits;
      gen->GetU01 = &SupDup96ADD_U01;
   }
   gen->Write   = &WrSupDup96;
   gen->param   = p;
   gen->state   = state;
   return gen;
}

unif01_Gen * umarsa_CreateSupDup96Add (unsigned int x0, unsigned int y0,
                                       unsigned int c)
{
   return CreateSupDup96 (x0, y0, c, '+');
}

unif01_Gen * umarsa_CreateSupDup96Xor (unsigned int x0, unsigned int y0,
                                       unsigned int c)
{
   return CreateSupDup96 (x0, y0, c, 'x');
}

/**************************************************************************/
#ifdef USE_LONGLONG

static unsigned long SupDup64ADD_Bits (void *vpar, void *vsta)
{
   SupDup64_state *state = vsta;
   SupDup64_param *param = vpar;
   ulonglong b;

   state->x = param->a * state->x + param->c;
   b = state->y ^ (state->y << param->s1);
   b ^= (b >> param->s2);
   state->y = b ^ (b << param->s3);
   return (state->y + state->x) >> 32;
}

/*-----------------------------------------------------------------------*/

static double SupDup64ADD_U01 (void *vpar, void *vsta)
{
   return SupDup64ADD_Bits (vpar, vsta) * INV32;
}

/*-----------------------------------------------------------------------*/

static unsigned long SupDup64XOR_Bits (void *vpar, void *vsta)
{
   SupDup64_state *state = vsta;
   SupDup64_param *param = vpar;
   ulonglong b;

   state->x = param->a * state->x + param->c;
   b = state->y ^ (state->y << param->s1);
   b ^= b >> param->s2;
   state->y = b ^ (b << param->s3);
   return (state->y ^ state->x) >> 32;
}

/*-----------------------------------------------------------------------*/

static double SupDup64XOR_U01 (void *vpar, void *vsta)
{
   return SupDup64XOR_Bits (vpar, vsta) * INV32;
}

/*-----------------------------------------------------------------------*/

static void WrSupDup64 (void *vsta)
{
   SupDup64_state *state = vsta;
   printf (" x = %1" PRIuLEAST64 ",    y = %1" PRIuLEAST64 "\n\n",
           state->x, state->y);
}

/*-----------------------------------------------------------------------*/

static unif01_Gen *CreateSupDup64 (ulonglong x, ulonglong y, ulonglong a,
    ulonglong c, int s1, int s2, int s3, char op)
{
   unif01_Gen *gen;
   SupDup64_state *state;
   SupDup64_param *param;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (SupDup64_state));
   param = util_Malloc (sizeof (SupDup64_param));
   state->x = x;
   state->y = y;
   param->a = a;
   param->c = c;
   param->s1 = s1;
   param->s2 = s2;
   param->s3 = s3;

   util_Assert ((a % 8 == 3) || (a % 8 == 5),
      "umarsa_CreateSupDup64:   a must be 3 mod 8  or  5 mod 8");
   if (op == '+')
      strcpy (name, "umarsa_CreateSupDup64Add:");
   else
      strcpy (name, "umarsa_CreateSupDup64Xor:");
   addstr_ULONG (name, "   x0 = ", x);
   addstr_ULONG (name, ",   y0 = ", y);
   addstr_ULONG (name, ",   a = ", a);
   addstr_ULONG (name, ",   c = ", c);
   addstr_Uint (name, ",   s1 = ", (unsigned int) s1);
   addstr_Uint (name, ",   s2 = ", (unsigned int) s2);
   addstr_Uint (name, ",   s3 = ", (unsigned int) s3);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if (op == 'x') {
      gen->GetBits = &SupDup64XOR_Bits;
      gen->GetU01 = &SupDup64XOR_U01;
   } else {
      gen->GetBits = &SupDup64ADD_Bits;
      gen->GetU01 = &SupDup64ADD_U01;
   }
   gen->Write   = &WrSupDup64;
   gen->param   = param;
   gen->state   = state;
   return gen;
}

unif01_Gen * umarsa_CreateSupDup64Add (ulonglong x, ulonglong y, ulonglong a,
    ulonglong c, int s1, int s2, int s3)
{
   return CreateSupDup64 (x, y, a, c, s1, s2, s3, '+');
}

unif01_Gen * umarsa_CreateSupDup64Xor (ulonglong x, ulonglong y, ulonglong a,
    ulonglong c, int s1, int s2, int s3)
{
   return CreateSupDup64 (x, y, a, c, s1, s2, s3, 'x');
}

#endif
/**************************************************************************/

static unsigned long Marsa90a_Bits (void *junk, void *vsta)
{
   Marsa90a_state *state = vsta;
   unsigned int x, y;

   x = state->X[state->s];
   y = state->X[state->r] + state->C;
   if (x < y) {
      x += (MA - y);
      state->C = 1;
   } else {
      x -= y;
      state->C = 0;
   }
   state->X[state->r] = x;
   state->r = (state->r + 1) % 43;
   state->s = (state->s + 1) % 43;

   state->W -= 362436069;
   return x - state->W;
}

/*-----------------------------------------------------------------------*/

static double Marsa90a_U01 (void *vpar, void *vsta)
{
   return Marsa90a_Bits (vpar, vsta) * INV32;
}

/*-----------------------------------------------------------------------*/

static void WrMarsa90a (void *vsta)
{
   Marsa90a_state *state = vsta;
   int j;

   if (unif01_WrLongStateFlag) {
      printf (" X = {\n");
      for (j = 0; j < 43; j++) {
         printf ("   %10u\n", state->X[j]);
      }
      printf ("   }\n\n Weyl:   W = %10u\n", state->W);
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateMarsa90a (int s1, int s2, int s3, int s4,
   unsigned int s5)
{
   unif01_Gen *gen;
   Marsa90a_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned int s;
   unsigned int t, M, j;

   /* For the generation of the 43 first values of the SWB generator, */
   /* ACARRYPC uses a combination of a 3-Lag Fibonacci generator */
   /* with a LCG generator LCG. */

   if ((s1 >= 179) || (s2 >= 179) || (s3 >= 179) ||
      (s1 <= 0) || (s2 <= 0) || (s3 <= 0) || (s4 >= 169))
      util_Error ("umarsa_CreateMarsa90a:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (Marsa90a_state));

   strcpy (name, "umarsa_CreateMarsa90a:");
   addstr_Uint (name, "   y1 = ", (unsigned) s1);
   addstr_Uint (name, ",   y2 = ", (unsigned) s2);
   addstr_Uint (name, ",   y3 = ", (unsigned) s3);
   addstr_Uint (name, ",   z0 = ", (unsigned) s4);
   addstr_Uint (name, ",   Y0 = ", s5);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->r = 0;
   state->s = 21;
   for (j = 0; j < 43; j++) {
      s = 0;
      for (t = 0; t < 32; t++) {
         M = (((s1 * s2) % 179) * s3) % 179;
         s1 = s2;
         s2 = s3;
         s3 = M;
         s4 = ((53 * s4) + 1) % 169;
         if (((s4 * M) % 64) >= 32)
            s |= 1 << t;
      }
      state->X[j] = s % MA;
   }
   M = (((s1 * s2) % 179) * s3) % 179;
   s1 = s2;
   s2 = s3;
   s3 = M;
   s4 = ((53 * s4) + 1) % 169;
   if (((s4 * M) % 64) >= 32)
      state->C = 1;
   else
      state->C = 0;
   state->W = s5;

   gen->GetBits = &Marsa90a_Bits;
   gen->GetU01  = &Marsa90a_U01;
   gen->Write   = &WrMarsa90a;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

static double RANMAR_U01 (void *vpar, void *vsta)
/*
 * Generator proposed by G. Marsaglia. Combines a lag-Fibonacci generator
 * F(97,33,-) with the following arithmetic sequence modulo (2^24 - 3): I
 * - k, I - 2k, I - 3k,... (for an initial integer I).
 */
{
   RANMAR_state *state = vsta;
   RANMAR_param *param = vpar;
   double temp;

   temp = state->X[state->r] - state->X[state->s];
   if (temp < 0.0)
      temp += 1.0;
   state->X[state->r] = temp;

   state->r--;
   if (state->r == 0)
      state->r = 97;

   state->s--;
   if (state->s == 0)
      state->s = 97;

   state->Seq -= param->SeqD;
   if (state->Seq < 0.0)
      state->Seq += param->SeqM;
   temp -= state->Seq;
   if (temp < 0.0)
      return temp + 1.0;
   return temp;
}

/*-----------------------------------------------------------------------*/

static unsigned long RANMAR_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (RANMAR_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrRANMAR (void *vsta)
{
   RANMAR_state *state = vsta;
   int j;

   if (unif01_WrLongStateFlag) {
      printf (" X = {\n");
      for (j = 1; j <= 97; j++) {
         printf ("  %12.9f\n", state->X[j]);
      }
      printf ("   }\n\n");
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateRANMAR (int s1, int s2, int s3, int s4)
/*
 * For the generation of its 97 first values, Marsa90 uses
 * a combination of a 3-Lag Fib. gen. with a LCG gen. 
 */
{
   unif01_Gen *gen;
   RANMAR_state *state;
   RANMAR_param *param;
   size_t leng;
   char name[LEN + 1];
   double s, t;
   unsigned int M, k, j;


   if ((s1 >= 179) || (s2 >= 179) || (s3 >= 179) ||
      (s1 <= 0) || (s2 <= 0) || (s3 <= 0) || (s4 >= 169)) {
      util_Error ("umarsa_CreateRANMAR:   Invalid parameter");
   }

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (RANMAR_state));
   param = util_Malloc (sizeof (RANMAR_param));

   strcpy (name, "umarsa_CreateRANMAR:");
   addstr_Uint (name, "   y1 = ", (unsigned) s1);
   addstr_Uint (name, ",   y2 = ", (unsigned) s2);
   addstr_Uint (name, ",   y3 = ", (unsigned) s3);
   addstr_Uint (name, ",   z0 = ", (unsigned) s4);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->r = 97;
   state->s = 33;
   state->Seq = 362436.0 / 16777216.0;
   param->SeqD = 7654321.0 / 16777216.0;
   param->SeqM = 16777213.0 / 16777216.0;

   for (j = 1; j <= 97; j++) {
      s = 0.0;
      t = 0.5;
      for (k = 1; k <= 24; k++) {
         M = (((s1 * s2) % 179) * s3) % 179;
         s1 = s2;
         s2 = s3;
         s3 = M;
         s4 = ((53 * s4) + 1) % 169;
         if (((s4 * M) % 64) >= 32)
            s += t;
         t *= 0.5;
      }
      state->X[j] = s;
   }

   gen->GetBits = &RANMAR_Bits;
   gen->GetU01  = &RANMAR_U01;
   gen->Write   = &WrRANMAR;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/**************************************************************************/
#ifdef USE_LONGLONG

static unsigned long Mother0_Bits (void *junk, void *vsta)
{
   Mother0_state *state = vsta;
   ulonglong V;

   V = 2111111111ULL * state->S1 + 1492ULL * state->S2 + 1776ULL * state->S3 +
       5115ULL * state->S4 + state->C;
   state->S1 = state->S2;
   state->S2 = state->S3;
   state->S3 = state->S4;
   state->S4 = V & MASK32;
   state->C = V >> 32;
   return state->S4;
}

/*-----------------------------------------------------------------------*/

static double Mother0_U01 (void *vpar, void *vsta)
{
   return Mother0_Bits (vpar, vsta) * INV32;
}

/*-----------------------------------------------------------------------*/

static void WrMother0 (void *vsta)
{
   Mother0_state *state = vsta;
   printf
      (" (S1, S2, S3, S4, C) = ( %1" PRIuLEAST64 ",  %1" PRIuLEAST64 ",  %1"
       PRIuLEAST64 ",\n" "                         %1" PRIuLEAST64
       ",  %1" PRIuLEAST64 " )\n",
      state->S1, state->S2, state->S3, state->S4, state->C);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateMother0 (unsigned long s1, unsigned long s2,
   unsigned long s3, unsigned long s4, unsigned long c)
{
   unif01_Gen *gen;
   Mother0_state *state;
   size_t leng;
   char name[LEN + 1];

   if (c > 2111119494) {
      util_Error ("umarsa_CreateMother0:   Invalid parameter");
   }

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (Mother0_state));

   strcpy (name, "umarsa_CreateMother0:");
   addstr_Ulong (name, "   x1 = ", s1);
   addstr_Ulong (name, ",   x2 = ", s2);
   addstr_Ulong (name, ",   x3 = ", s3);
   addstr_Ulong (name, ",   x4 = ", s4);
   addstr_Ulong (name, ",   c = ", c);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S1 = s1;
   state->S2 = s2;
   state->S3 = s3;
   state->S4 = s4;
   state->C = c;

   gen->GetBits = &Mother0_Bits;
   gen->GetU01  = &Mother0_U01;
   gen->Write   = &WrMother0;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}


#endif
/*****************************************************************************/

static unsigned long Combo_Bits (void *junk, void *vsta)
{
   Combo_state *state = vsta;
   unsigned long V, Z;

   V = (state->X1 * state->X2) & MASK32;
   state->X1 = state->X2;
   state->X2 = V;
   Z = state->Y1;
   state->Y1 = 30903 * (MASK16 & Z);
   state->Y1 = (state->Y1 + (Z >> 16)) & MASK32;
   return (state->X2 + state->Y1) & MASK32;
}

/*-----------------------------------------------------------------------*/

static double Combo_U01 (void *vpar, void *vsta)
{
   return Combo_Bits (vpar, vsta) * INV32;
}

/*-----------------------------------------------------------------------*/

static void WrCombo (void *vsta)
{
   Combo_state *state = vsta;
   printf (" (x1, x2, y1, c) = ( %1lu, %1lu,  %1lu, %1lu )\n\n",
      state->X1, state->X2, state->Y1 % 65536, state->Y1 / 65536);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateCombo (unsigned int x1, unsigned int x2,
   unsigned int y1, unsigned int c)
{
   unif01_Gen *gen;
   Combo_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((y1 >= 65536) || (c > 30903)) {
      util_Error ("umarsa_CreateCombo:   Invalid parameter");
   }
   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (Combo_state));

   strcpy (name, "umarsa_CreateCombo:");
   addstr_Uint (name, "   x1 = ", x1);
   addstr_Uint (name, ",   x2 = ", x2);
   addstr_Uint (name, ",   y1 = ", y1);
   addstr_Uint (name, ",   c = ", c);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->X1 = 2 * x1 + 1;
   state->X1 = 3 * state->X1 * state->X1;
   state->X2 = 2 * x2 + 1;
   state->Y1 = y1 + c;

   gen->GetBits = &Combo_Bits;
   gen->GetU01  = &Combo_U01;
   gen->Write   = &WrCombo;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}


/***************************************************************************/

static double ECG1_U01 (void *vpar, void *vsta)
{
   ECG_state *state = vsta;
   ECG_param *param = vpar;
   double y;

   y = 65065.0 * state->X1 + 67067.0 * state->X2 + 69069.0 * state->X3;
   y = fmod (y, param->Norm);
   state->X1 = state->X2;
   state->X2 = state->X3;
   state->X3 = y;
   return (y / param->Norm);
}

/*-----------------------------------------------------------------------*/

static unsigned long ECG1_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (ECG1_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static double ECG2_U01 (void *vpar, void *vsta)
{
   ECG_param *param = vpar;
   ECG_state *state = vsta;
   double y;

   y = 1024.0 * (state->X1 + state->X2 + state->X3);
   y = fmod (y, param->Norm);
   state->X1 = state->X2;
   state->X2 = state->X3;
   state->X3 = y;
   return (y / param->Norm);
}

/*-----------------------------------------------------------------------*/

static unsigned long ECG2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (ECG2_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static double ECG3_U01 (void *vpar, void *vsta)
{
   ECG_param *param = vpar;
   ECG_state *state = vsta;
   double y;

   y = 2000.0 * state->X1 + 1950.0 * state->X2 + 1900.0 * state->X3;
   y = fmod (y, param->Norm);
   state->X1 = state->X2;
   state->X2 = state->X3;
   state->X3 = y;
   return (y / param->Norm);
}

/*-----------------------------------------------------------------------*/

static unsigned long ECG3_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (ECG3_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static double ECG4_U01 (void *vpar, void *vsta)
{
   ECG_param *param = vpar;
   ECG_state *state = vsta;
   double y;

   y = 1048576.0 * (state->X1 + state->X2 + state->X3);
   y = fmod (y, param->Norm);
   state->X1 = state->X2;
   state->X2 = state->X3;
   state->X3 = y;
   return (y / param->Norm);
}

/*-----------------------------------------------------------------------*/

static unsigned long ECG4_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (ECG4_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrECG (void *vsta)
{
   ECG_state *state = vsta;
   printf (" X1 = %12.0f,    X2 = %12.0f,    X3 = %12.0f\n\n",
      state->X1, state->X2, state->X3);
}

/*-----------------------------------------------------------------------*/

static unif01_Gen *CreateECG (char *nom, unsigned int x1, unsigned int x2,
   unsigned int x3)
{
   unif01_Gen *gen;
   ECG_param *param;
   ECG_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (ECG_state));
   param = util_Malloc (sizeof (ECG_param));

   strcpy (name, nom);
   addstr_Uint (name, "   x1 = ", x1);
   addstr_Uint (name, ",   x2 = ", x2);
   addstr_Uint (name, ",   x3 = ", x3);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->X1 = x1;
   state->X2 = x2;
   state->X3 = x3;

   gen->Write = &WrECG;
   gen->param = param;
   gen->state = state;
   return gen;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateECG1 (unsigned int x1, unsigned int x2,
   unsigned int x3)
{
   unif01_Gen *gen;
   ECG_param *param;

   gen = CreateECG ("umarsa_CreateECG1:", x1, x2, x3);
   param = gen->param;
   param->Norm = 4294967291.0;
   gen->GetBits = &ECG1_Bits;
   gen->GetU01  = &ECG1_U01;
   return gen;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateECG2 (unsigned int x1, unsigned int x2,
   unsigned int x3)
{
   unif01_Gen *gen;
   ECG_param *param;

   gen = CreateECG ("umarsa_CreateECG2:", x1, x2, x3);
   param = gen->param;
   param->Norm = 4294967291.0;
   gen->GetBits = &ECG2_Bits;
   gen->GetU01  = &ECG2_U01;
   return gen;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateECG3 (unsigned int x1, unsigned int x2,
   unsigned int x3)
{
   unif01_Gen *gen;
   ECG_param *param;

   gen = CreateECG ("umarsa_CreateECG3:", x1, x2, x3);
   param = gen->param;
   param->Norm = 4294967087.0;
   gen->GetBits = &ECG3_Bits;
   gen->GetU01  = &ECG3_U01;
   return gen;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateECG4 (unsigned int x1, unsigned int x2,
   unsigned int x3)
{
   unif01_Gen *gen;
   ECG_param *param;

   gen = CreateECG ("umarsa_CreateECG4:", x1, x2, x3);
   param = gen->param;
   param->Norm = 4294967087.0;
   gen->GetBits = &ECG4_Bits;
   gen->GetU01  = &ECG4_U01;
   return gen;
}


/*****************************************************************************/

static unsigned long MWC97R_Bits (void *junk, void *sta)
{
   MWC97R_state *state = sta;
   state->I1 = 36969 * (state->I1 & MASK16) + (state->I1 >> 16);
   state->I2 = 18000 * (state->I2 & MASK16) + (state->I2 >> 16);
   return (state->I1 << 16) + (state->I2 & MASK16);
}

/*-----------------------------------------------------------------------*/

static double MWC97R_U01 (void *par, void *sta)
{
   return MWC97R_Bits (par, sta) * 2.328306437080797e-10;
}

/*-----------------------------------------------------------------------*/

static void WrMWC97R (void *sta)
{
   MWC97R_state *state = sta;
   printf (" I1 = %u,   I2 = %u\n", state->I1, state->I2);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateMWC97R (unsigned int I1, unsigned int I2)
{
   unif01_Gen *gen;
   MWC97R_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   gen->state = state = util_Malloc (sizeof (MWC97R_state));
   state->I1 = I1;
   state->I2 = I2;
   gen->param = NULL;
   gen->Write = WrMWC97R;
   gen->GetU01 = MWC97R_U01;
   gen->GetBits = MWC97R_Bits;

   strcpy (name, "umarsa_CreateMWC97R:");
   addstr_Uint (name, "   x0 = ", I1);
   addstr_Uint (name, ",   y0 = ", I2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);
   return gen;
}


/***************************************************************************/

static unsigned long LFIB4_99_Bits (void *junk, void *vsta)
{
   LFIB4_99_state *state = vsta;

   ++state->r;
   state->r &= MASK8;
   state->X[state->r] += state->X[(state->r + 58) & MASK8]
    + state->X[(state->r + 119) & MASK8] + state->X[(state->r + 178) & MASK8];
   return state->X[state->r];
}

/*-------------------------------------------------------------------------*/

static double LFIB4_99_U01 (void *vpar, void *vsta)
{
   return LFIB4_99_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void WrLFIB4_99 (void *vsta)
{
   LFIB4_99_state *state = vsta;
   unsigned int j;

   if (unif01_WrLongStateFlag) {
      printf ("T = {\n");
      for (j = 0; j < 256; j++) {
         printf (" %12u", state->X[(state->r + j) & MASK8]);
         if (j < 255)
            printf (",");
         if (((j + 1) % 5) == 0)
            printf ("\n");
      };
      printf ("\n};\n");
   } else
      unif01_WrLongStateDef ();
}


/*-------------------------------------------------------------------------*/

unif01_Gen *umarsa_Create4LFIB99 (unsigned int T[256])
{
   unif01_Gen *gen;
   LFIB4_99_state *state;
   size_t leng;
   char name[LEN + 1];
   int j;

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (LFIB4_99_state));

   strcpy (name, "umarsa_Create4LFIB99:");
   addstr_ArrayUint (name, "   T = ", 256, T);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->r = 0;
   gen->GetBits = &LFIB4_99_Bits;
   gen->GetU01 = &LFIB4_99_U01;

   for (j = 0; j < 256; j++) {
      state->X[j] = T[j];
   }
   gen->param = NULL;
   gen->state = state;
   gen->Write = &WrLFIB4_99;

   return gen;
}


/***************************************************************************/

static unsigned long SHR3_99_Bits (void *junk, void *sta)
{
   unsigned long *I1 = sta;
   unsigned long t = *I1;

   t ^= t << 17;
#ifndef IS_ULONG32
   t &= MASK32;
#endif
   t ^= t >> 13;
   t ^= t << 5;
   *I1 = t & MASK32;
   return *I1;
}

/*-----------------------------------------------------------------------*/

static double SHR3_99_U01 (void *par, void *sta)
{
   return SHR3_99_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

static void WrSHR3_99 (void *sta)
{
   unsigned long *I1 = sta;
   printf (" x = %lu\n", *I1);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *umarsa_Create3SHR99 (unsigned int I1)
{
   unif01_Gen *gen;
   unsigned long *p1;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   p1 = gen->state = util_Malloc (sizeof (unsigned long));
   *p1 = I1;
   gen->param = NULL;
   gen->Write = WrSHR3_99;
   gen->GetU01 = SHR3_99_U01;
   gen->GetBits = SHR3_99_Bits;

   strcpy (name, "umarsa_Create3SHR99:");
   addstr_Uint (name, "   x0 = ", I1);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);
   return gen;
}


/***************************************************************************/

static unsigned long SWB_99_Bits (void *junk, void *vsta)
{
   SWB_99_state *state = vsta;

   ++state->r;
   state->r &= MASK8;
   state->b = (state->X[(state->r + 33) & MASK8] <
      state->X[(state->r + 18) & MASK8] + state->b) ? 1 : 0;
   state->X[state->r] = state->X[(state->r + 34) & MASK8]
      - state->X[(state->r + 19) & MASK8] - state->b;
   return state->X[state->r];
}

/*-------------------------------------------------------------------------*/

static double SWB_99_U01 (void *vpar, void *vsta)
{
   return SWB_99_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void WrSWB_99 (void *vsta)
{
   SWB_99_state *state = vsta;
   unsigned int j;

   if (unif01_WrLongStateFlag) {
      printf ("b = %u\n", state->b);
      printf ("T = {\n");
      for (j = 0; j < 256; j++) {
         printf (" %12u", state->X[(state->r + j) & MASK8]);
         if (j < 255)
            printf (",");
         if (((j + 1) % 5) == 0)
            printf ("\n");
      };
      printf ("\n};\n");
   } else
      unif01_WrLongStateDef ();
}


/*-------------------------------------------------------------------------*/

unif01_Gen *umarsa_CreateSWB99 (unsigned int T[256], int b)
{
   unif01_Gen *gen;
   SWB_99_state *state;
   size_t leng;
   char name[LEN + 1];
   int j;

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (SWB_99_state));

   strcpy (name, "umarsa_CreateSWB99:");
   addstr_Uint (name, "   b = ", (unsigned) b);
   addstr_ArrayUint (name, ",   T = ", 256, T);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->r = 0;
   state->b = b;
   gen->GetBits = &SWB_99_Bits;
   gen->GetU01 = &SWB_99_U01;

   for (j = 0; j < 256; j++) {
      state->X[j] = T[j];
   }
   gen->param = NULL;
   gen->state = state;
   gen->Write = &WrSWB_99;

   return gen;
}


void umarsa_DeleteGen (unif01_Gen *gen)
{
   unif01_DeleteGen (gen);
}
