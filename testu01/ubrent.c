/*************************************************************************\
 *
 * Package:        TestU01
 * File:           ubrent.c
 * Environment:    ANSI C
 *
 * Comment: I have made mostly cosmetic changes on Brent's code to adapt
 *    his generators for testing by TestU01. The random numbers
 *    generated should be the same. The code was taken from Brent's file
 *    xorgens201.c available on his Web site for version 2004 of the RNG,
 *    and from xorgens304 for the version 2006.
 *
 * Richard P. Brent kindly gave us the explicit permission to include
 *    his code in TestU01 and distribute it. (Richard Simard)     
 *
\*************************************************************************/


#include "util.h"
#include "addstr.h"
#include "ubrent.h"

#include <stdio.h>
#include <string.h>

#define LEN  200                  /* Max length of strings */


typedef struct {
   unsigned int r, s, a, b, c, d;
   unsigned int Mask;
   unsigned int weil;
   lebool hasWeyl;
} Xorgen32_param;


typedef struct {
   unsigned int *x;
   unsigned int w, r2;
   int i;
} Xorgen32_state;


#ifdef USE_LONGLONG
typedef struct {
   unsigned int r, s, a, b, c, d;
   unsigned int Mask;
   ulonglong weil;
   lebool hasWeyl;
} Xorgen64_param;


typedef struct {
   ulonglong *x;
   ulonglong w;
   unsigned int r2;
   int i;
} Xorgen64_state;
#endif


static int co1 = 0, co2 = 0, co3 = 0, co4 = 0, co5 = 0;      /* Counters */


/* xorgens.c - Some long-period random number generators
               generalising Marsaglia's Xorshift RNGs
/*
==========================================================================
|                                                                        |
|  Copyright (C) 2004, 2006 R. P. Brent.                                 |
|                                                                        |
|  This program is free software; you can redistribute it and/or         |
|  modify it under the terms of the GNU General Public License,          |
|  version 2, June 1991, as published by the Free Software Foundation.   |
|  For details see http://www.gnu.org/copyleft/gpl.html .                |
|                                                                        |
|  If you would like to use this software but the GNU GPL creates legal  |
|  problems, then please contact the author to negotiate a special       |
|  agreement.                                                            |
|                                                                        |
==========================================================================
*/
/*
 
Author:         Richard P. Brent (random@rpbrent.co.uk)

Version:        2.01

Contents:       xor4096s        "unsigned long" 32-bit integer RNG,
                                period (2^4096 - 1)*2^32.
                                
                xor4096f        "float" 32-bit real RNG (based on xor4096s).

                xor4096l        "unsigned long long" 64-bit integer RNG,
                                period (2^4096 - 1)*2^64.

                xor4096d        "double" 64-bit real RNG (based on xor4096l).


Comments:       Some fast RNGs with very long periods, using minimal
                storage.  They use generalisations of (but not exactly the
                same as) the xorshift RNGs described by Marsaglia[4].  The
                output is combined with a "Weil generator" although
                this extra precaution is probably unnecessary.
                
                The generators are fast, although quality has not been
                sacrificed for speed.  To speed them up further, the code
                could easily be modified to return an array of random numbers
                instead of a single number, thus reducing function call
                overheads. For the sake of simplicity I have not done this.
                
                The generators are believed to be of high quality, and have
                passed Marsaglia's DIEHARD tests. Of course, users should
                apply their own tests.  Please inform the author if you
                discover any genuine problems.
*/


/*========================================================================*/
   /*
unsigned long xor4096s (unsigned long seed)

   32-bit random number generator with period at least 2^4096-1.

      Is is assumed that "unsigned long" is a 32-bit integer.

      The method is a generalisation of Marsaglia's xorgen generators, see G.
      Marsaglia, "Xorshift RNGs", JSS 8, 14 (2003), Sec. 3.1 (available from
      http://www.jstatsoft.org/ ).

      The primary recurrence is x_k = x_{k-r}A + x_{k-s}B, where

      A = (I + L^a)(I + R^b), B = (I + L^c)(I + R^d)

      in the notation of Marsaglia [JSS 8,14 (2003)].

      We choose r > 1, 0 < s < r, n = 32*r, and (a, b, c, d) such that the n
      times n matrix T defining the recurrence has minimal polynomial P which
      is of degree n and primitive over GF(2).  Thus the period is 2^n-1 and
      all possible n-bit sequences (except all zeros) occur in the output over 
      a full period.  Storage is (r + constant) 32-bit words.

      Our generalisation is:

      (1) The block matrix has r times r blocks, where r is a parameter. To
      test primitivity we need to be able to factor 2^n-1.  Marsaglia
      considered some cases with small r.  In our implementation r is a power
      of 2, but this restriction is easily removed. We use a circular array,
      which makes the number of memory references independent of r.

      (2) Marsaglia's block (I + R^c) is replaced by (I + L^c)(I + R^d).

      (3) The block mentioned in (2) is in position (r+1-s,r) instead of
      position (r,r) of the block matrix. (Marsaglia considers the case s =
      1.)

      By introducing the additional parameters d and s we increase the set of
      candidate recurrences to compensate for the fact that primitive
      recurrences get harder to find as n increases.

      The output of the primary recurrence is combined with a Weil generator
      to avoid problems with the binary rank test and correlations related to
      the sparsity of T. This increases the period by a factor of 2^32.

      Should be called once with nonzero seed, thereafter with zero seed.

      The generator implemented here has parameters (wlen 32, r 128, s 95, a
      17, b 12, c 13, d 15) with Weight(P) 251. Other parameters can be used
      by changing the definitions below.

      R. P. Brent, 20040720, last revised 20040802.
    */

#define wlen 32
#define r 128
#define s 95
#define a 17
#define b 12
#define c 13
#define d 15

static unsigned long xor4096s_Bits (void *junk, void *vsta)
{
   static unsigned int w, x[r], weil = 0x61c88647;
   unsigned int t, v;
   static int i = -1;             /* i < 0 indicates first call */
   int k;
   if (i < 0) {  /* Initialisation necessary */
      int seed = *((unsigned int *) vsta);
      v = (seed != 0) ? seed : ~seed; /* v must be nonzero */

      for (k = wlen; k > 0; k--)  /* Avoid correlations for close seeds */
         v ^= (v ^= (v ^= v << 13) >> 17) << 5; /* This recurrence has period
                                                   2^32-1 */

      for (w = v, k = 0; k < r; k++)
         /* Initialise circular array */
         x[k] = (v ^= (v ^= (v ^= v << 13) >> 17) << 5) + (w += weil);

      for (i = r - 1, k = 4 * r; k > 0; k--) { /* Discard first 4*r results
                                                  (Gimeno) */
         t = x[i = (i + 1) & (r - 1)];
         v = x[(i + (r - s)) & (r - 1)];
         t ^= (t ^= t << a) >> b;
         v ^= v << c;
         x[i] = (v ^= t ^ (v >> d));
      }
   }

   /* Apart from initialisation (above), this is the generator */
   t = x[i = (i + 1) & (r - 1)];  /* Increment i mod r (r is a power of 2) */
   v = x[(i + (r - s)) & (r - 1)]; /* Index is (i - s) mod r */
   t ^= (t ^= t << a) >> b;       /* (I + L^a)(I + R^b) */
   v ^= v << c;                   /* I + L^c */
   x[i] = (v ^= t ^ (v >> d));    /* Update circular array */
   return (v + (w += weil));      /* Return combination with Weil generator */
}

#undef wlen
#undef r
#undef s
#undef a
#undef b
#undef c
#undef d


/*-----------------------------------------------------------------------*/

static double xor4096s_U01 (void *vpar, void *vsta)
{
   return xor4096s_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

static void Wr_xor4096 (void *junk)
{
}

/*-----------------------------------------------------------------------*/

unif01_Gen * ubrent_CreateXor4096s (unsigned int seed)
{
   unif01_Gen *gen;
   unsigned int *pseed;
   size_t leng;
   char name[LEN + 1];

   util_Assert (co1 == 0,
   "ubrent_CreateXor4096s:\n   only 1 such generator can be used at a time");
   co1++;

   gen = util_Malloc (sizeof (unif01_Gen));
   pseed = util_Malloc (sizeof (unsigned int));
   *pseed = seed;
   strcpy (name, "ubrent_CreateXor4096s:");
   addstr_Uint (name, "   seed = ", seed);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &xor4096s_Bits;
   gen->GetU01 = &xor4096s_U01;
   gen->Write = &Wr_xor4096;
   gen->param = NULL;
   gen->state = pseed;
   /* Discard the first value
   xor4096s_Bits (NULL, pseed); */
   return gen;
}


/*-----------------------------------------------------------------------*/

void ubrent_DeleteXor4096s (unif01_Gen *gen)
{
   if (NULL == gen)
      return;
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
   co1--;
}


/*========================================================================*/

static unsigned long Xorgen32_Bits (void *vpar, void *vsta)
{
   Xorgen32_param *param = vpar;
   Xorgen32_state *state = vsta;
   unsigned int t, v;

   /* Increment i mod ... */
   t = state->x[state->i = (state->i + 1) & param->Mask];

   /* Index is (i - s) mod ... */
   v = state->x[(state->i + (param->r - param->s)) & param->Mask];

   t ^= (t ^= t << param->a) >> param->b;       /* (I + L^a)(I + R^b) */
   v ^= v << param->c;                          /* I + L^c */

   /* Update circular array */
   state->x[state->i] = (v ^= t ^ (v >> param->d));

   if (param->hasWeyl)
      /* combination with Weil generator */
      return (v + (state->w += param->weil));
   else
      return v;
}


/*-----------------------------------------------------------------------*/

static double Xorgen32_U01 (void *vpar, void *vsta)
{
   return Xorgen32_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

static void WrXorgen32 (void *vsta)
{
   Xorgen32_state *state = vsta;
   unsigned int j;
   unsigned int s = state->i;

   if (unif01_WrLongStateFlag) {
      printf (" i = %u,   w = %u\n", state->i, state->w);
      printf (" x = {\n ");
      for (j = 0; j < state->r2; j++) {
         if (++s >= state->r2)
            s = 0;
         printf (" %12u", state->x[s]);
         if (j < state->r2 - 1)
            printf (",");
         if ((j % 5) == 4)
            printf ("\n ");
      };
      printf ("    }\n");
   } else
      unif01_WrLongStateDef ();

}

/*-----------------------------------------------------------------------*/

unif01_Gen * ubrent_CreateXorgen32 (int r, int s, int a, int b, int c, int d,
                                    lebool hasWeyl, unsigned int seed)
{
   unif01_Gen *gen;
   Xorgen32_param *param;
   Xorgen32_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned int t, v;
   const unsigned int wlen = 32;
   int i, k;

   util_Assert (r > 1, "ubrent_CreateXorgen32:   r <= 1");
   util_Assert (s > 0, "ubrent_CreateXorgen32:   s <= 0");
   util_Assert (r > s, "ubrent_CreateXorgen32:   r <= s");
   util_Assert (a < 32 && b < 32 && c < 32 && d < 32,
      "ubrent_CreateXorgen32:   one of a, b, c, d >= 32");
   util_Assert (a > 0 && b > 0 && c > 0 && d > 0,
      "ubrent_CreateXorgen32:   one of a, b, c, d <= 0");

   i = 1;
   while (i < r)
      i <<= 1;
   util_Assert (r == i, "ubrent_CreateXorgen32:   r is not a power of 2");

   gen = util_Malloc (sizeof (unif01_Gen));

   strcpy (name, "ubrent_CreateXorgen32:");
   addstr_Int (name, "   r = ", r);
   addstr_Int (name, ",  s = ", s);
   addstr_Int (name, ",  a = ", a);
   addstr_Int (name, ",  b = ", b);
   addstr_Int (name, ",  c = ", c);
   addstr_Int (name, ",  d = ", d);
   if (hasWeyl)
      strncat (name, ",  hasWeyl = TRUE", 20);
   else
      strncat (name, ",  hasWeyl = FALSE", 20);
   addstr_Uint (name, ",  seed = ", seed);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param = util_Malloc (sizeof (Xorgen32_param));
   state = util_Malloc (sizeof (Xorgen32_state));

   gen->GetBits = &Xorgen32_Bits;
   gen->GetU01 = &Xorgen32_U01;
   gen->Write = &WrXorgen32;
   gen->param = param;
   gen->state = state;
   
   param->hasWeyl = hasWeyl;
   param->Mask = r - 1;
   state->r2 = r;
   state->x = util_Calloc ((size_t) r, sizeof (unsigned int));
   
   v = (seed != 0) ? seed : ~seed; /* v must be nonzero */

   for (k = wlen; k > 0; k--)  /* Avoid correlations for close seeds */
      v ^= (v ^= (v ^= v << 13) >> 17) << 5; /* This recurrence has period
                                                2^32-1 */
   if (hasWeyl) {
      param->weil = 0x61c88647;
      for (state->w = v, k = 0; k < r; k++)
         /* Initialise circular array */
         state->x[k] = (v ^= (v ^= (v ^= v << 13) >> 17) << 5) +
            (state->w += param->weil);
   } else {
      param->weil = 0;
      for (k = 0; k < r; k++)
         /* Initialise circular array */
         state->x[k] = (v ^= (v ^= (v ^= v << 13) >> 17) << 5);
   }

   for (i = r - 1, k = 4 * r; k > 0; k--) { /* Discard first 4*r results
                                                      (Gimeno) */
      t = state->x[i = (i + 1) & param->Mask];
      v = state->x[(i + (r - s)) & param->Mask];
      t ^= (t ^= t << a) >> b;
      v ^= v << c;
      state->x[i] = (v ^= t ^ (v >> d));
   }

   state->i = i;
   param->a = a;
   param->b = b;
   param->c = c;
   param->d = d;
   param->r = r;
   param->s = s;

   return gen;
}


/*-----------------------------------------------------------------------*/

void ubrent_DeleteXorgen32 (unif01_Gen *gen)
{
   Xorgen32_state *state;
   if (NULL == gen)
      return;
   state = gen->state;
   util_Free (state->x);
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/*========================================================================*/
#ifdef USE_LONGLONG

/* unsigned long long xor4096l (unsigned long long seed)

      64-bit random number generator with period at least 2^4096-1.

      Is is assumed that "unsigned long long" is a 64-bit integer.

      The method is a generalisation of Marsaglia's xorgen generators, see G.
      Marsaglia, "Xorshift RNGs", JSS 8, 14 (2003), Sec. 3.1 (available from
      http://www.jstatsoft.org/ ).

      The primary recurrence is x_k = x_{k-r}A + x_{k-s}B, where

      A = (I + L^a)(I + R^b), B = (I + L^c)(I + R^d)

      in the notation of Marsaglia [JSS 8,14 (2003)].

      We choose r > 1, 0 < s < r, n = 64*r, and (a, b, c, d) such that the n
      times n matrix T defining the recurrence has minimal polynomial P which
      is of degree n and primitive over GF(2).  Thus the period is 2^n-1 and
      all possible n-bit sequences (except all zeros) occur in the output over 
      a full period.  Storage is (r + constant) 64-bit words.

      Our generalisation is:

      (1) The block matrix has r times r blocks, where r is a parameter. To
      test primitivity we need to be able to factor 2^n-1.  Marsaglia
      considered some cases with small r.  In our implementation r is a power
      of 2, but this restriction is easily removed. We use a circular array,
      which makes the number of memory references independent of r.

      (2) Marsaglia's block (I + R^c) is replaced by (I + L^c)(I + R^d).

      (3) The block mentioned in (2) is in position (r+1-s,r) instead of
      position (r,r) of the block matrix. (Marsaglia considers the case s =
      1.)

      By introducing the additional parameters d and s we increase the set of
      candidate recurrences to compensate for the fact that primitive
      recurrences get harder to find as n increases.

      The output of the primary recurrence is combined with a Weil generator
      to avoid problems with the binary rank test and correlations related to
      the sparsity of T. This increases the period by a factor of 2^64.

      Should be called once with nonzero seed, thereafter with zero seed.

      The generator implemented here has parameters (wlen 64, r 64, s 53, a
      33, b 26, c 27, d 29) with Weight(P) 961. Other parameters can be used
      by changing the definitions below.

      R. P. Brent, 20040721, last revised 20040802.
    */

#define wlen 64
#define r 64
#define s 53
#define a 33
#define b 26
#define c 27
#define d 29

static unsigned long xor4096l_Bits (void *junk, void *vsta)
{
   static ulonglong w, x[r],
      weil = ((longlong) 0x61c88646 << 32) + (longlong) 0x80b583eb;
   ulonglong t, v;
   static int i = -1;             /* i < 0 indicates first call */
   int k;
   if (i < 0) {  /* Initialisation necessary */
      int seed = *((ulonglong *) vsta);
      v = (seed != 0) ? seed : ~seed; /* v must be nonzero */

      for (k = wlen; k > 0; k--)  /* Avoid correlations for close seeds */
         v ^= (v ^= v << 7) >> 9; /* This recurrence has period 2^64-1 */

      for (w = v, k = 0; k < r; k++)
         /* Initialise circular array */
         x[k] = (v ^= (v ^= v << 7) >> 9) + (w += weil);

      for (i = r - 1, k = 4 * r; k > 0; k--) { /* Discard first 4*r results
                                                  (Gimeno) */
         t = x[i = (i + 1) & (r - 1)];
         v = x[(i + (r - s)) & (r - 1)];
         t ^= (t ^= t << a) >> b;
         v ^= v << c;
         x[i] = (v ^= t ^ (v >> d));
      }
   }

   /* Apart from initialisation (above), this is the generator */

   t = x[i = (i + 1) & (r - 1)];  /* Increment i mod r (r is a power of 2) */
   v = x[(i + (r - s)) & (r - 1)]; /* Index is (i - s) mod r */
   t ^= (t ^= t << a) >> b;       /* (I + L^a)(I + R^b) */
   v ^= v << c;                   /* I + L^c */
   x[i] = (v ^= t ^ (v >> d));    /* Update circular array */
   t = (v + (w += weil));      /* Return combination with Weil generator */
   return (unsigned long) (t >> 32);
}

#undef wlen
#undef r
#undef s
#undef a
#undef b
#undef c
#undef d


/*-----------------------------------------------------------------------*/

static double xor4096l_U01 (void *vpar, void *vsta)
{
   return xor4096l_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

unif01_Gen * ubrent_CreateXor4096l (ulonglong seed)
{
   unif01_Gen *gen;
   ulonglong *pseed;
   size_t leng;
   char name[LEN + 1];

   util_Assert (co2 == 0,
   "ubrent_CreateXor4096l:\n   only 1 such generator can be used at a time");
   co2++;

   gen = util_Malloc (sizeof (unif01_Gen));
   pseed = util_Malloc (sizeof (ulonglong));
   *pseed = seed;
   strcpy (name, "ubrent_CreateXor4096l:");
   addstr_ULONG (name, "   seed = ", seed);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &xor4096l_Bits;
   gen->GetU01 = &xor4096l_U01;
   gen->Write = &Wr_xor4096;
   gen->param = NULL;
   gen->state = pseed;
   /* Discard the first value
   xor4096l_Bits (NULL, pseed); */
   return gen;
}


/*-----------------------------------------------------------------------*/

void ubrent_DeleteXor4096l (unif01_Gen *gen)
{
   if (NULL == gen)
      return;
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
   co2--;
}


/*========================================================================*/

static unsigned long Xorgen64_Bits (void *vpar, void *vsta)
{
   Xorgen64_param *param = vpar;
   Xorgen64_state *state = vsta;
   ulonglong t, v;

   /* Increment i mod r */
   t = state->x[state->i = (state->i + 1) & param->Mask];

   /* Index is (i - s) mod r */
   v = state->x[(state->i + (param->r - param->s)) & param->Mask];

   t ^= (t ^= t << param->a) >> param->b;       /* (I + L^a)(I + R^b) */
   v ^= v << param->c;                          /* I + L^c */

   /* Update circular array */
   state->x[state->i] = (v ^= t ^ (v >> param->d));

   if (param->hasWeyl)
      /* combination with Weil generator */
      return (v + (state->w += param->weil)) >> 32;
   else
      return v >> 32;
}


/*-----------------------------------------------------------------------*/

static double Xorgen64_U01 (void *vpar, void *vsta)
{
   return Xorgen64_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

static void WrXorgen64 (void *vsta)
{
   Xorgen64_state *state = vsta;
   unsigned int j;
   unsigned int s = state->i;

   if (unif01_WrLongStateFlag) {
      printf (" i = %d,   w = %llu\n", state->i, state->w);
      printf (" x = {\n ");
      for (j = 0; j < state->r2; j++) {
         if (++s >= state->r2)
            s = 0;
         printf ("  %20llu", state->x[s]);
         if (j < state->r2 - 1)
            printf (",");
         if ((j % 3) == 2)
            printf ("\n ");
      };
      printf ("    }\n");
   } else
      unif01_WrLongStateDef ();

}

/*-----------------------------------------------------------------------*/

unif01_Gen * ubrent_CreateXorgen64 (int r, int s, int a, int b, int c, int d,
                                    lebool hasWeyl, ulonglong seed)
{
   unif01_Gen *gen;
   Xorgen64_param *param;
   Xorgen64_state *state;
   size_t leng;
   char name[LEN + 1];
   ulonglong t, v;
   const unsigned int wlen = 64;
   int i, k;

   util_Assert (r > 1, "ubrent_CreateXorgen64:   r <= 1");
   util_Assert (s > 0, "ubrent_CreateXorgen64:   s <= 0");
   util_Assert (r > s, "ubrent_CreateXorgen64:   r <= s");
   util_Assert (a < 64 && b < 64 && c < 64 && d < 64,
      "ubrent_CreateXorgen64:   one of a, b, c, d >= 64");
   util_Assert (a > 0 && b > 0 && c > 0 && d > 0,
      "ubrent_CreateXorgen64:   one of a, b, c, d <= 0");
   i = 1;
   while (i < r)
      i <<= 1;
   util_Assert (r == i, "ubrent_CreateXorgen64:   r is not a power of 2");

   gen = util_Malloc (sizeof (unif01_Gen));

   strcpy (name, "ubrent_CreateXorgen64:");
   addstr_Int (name, "   r = ", r);
   addstr_Int (name, ",  s = ", s);
   addstr_Int (name, ",  a = ", a);
   addstr_Int (name, ",  b = ", b);
   addstr_Int (name, ",  c = ", c);
   addstr_Int (name, ",  d = ", d);
   strncat (name, ",  hasWeyl = ", 20);
   if (hasWeyl)
      strncat (name, "TRUE", 5);
   else
      strncat (name, "FALSE", 5);
   addstr_ULONG (name, ",  seed = ", seed);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param = util_Malloc (sizeof (Xorgen64_param));
   state = util_Malloc (sizeof (Xorgen64_state));

   gen->GetBits = &Xorgen64_Bits;
   gen->GetU01 = &Xorgen64_U01;
   gen->Write = &WrXorgen64;
   gen->param = param;
   gen->state = state;
   
   param->hasWeyl = hasWeyl;
   state->r2 = r;
   param->Mask = state->r2 - 1;
   state->x = util_Calloc ((size_t) state->r2, sizeof (ulonglong));

   v = (seed != 0) ? seed : ~seed; /* v must be nonzero */ 
   for (k = wlen; k > 0; k--)  /* Avoid correlations for close seeds */
      v ^= (v ^= v << 7) >> 9; /* This recurrence has period 2^64-1 */

   if (hasWeyl) {
      param->weil = ((longlong) 0x61c88646 << 32) + (longlong) 0x80b583eb;
      for (state->w = v, k = 0; k < r; k++)
         /* Initialise circular array */
         state->x[k] = (v ^= (v ^= v << 7) >> 9) + (state->w += param->weil);
   } else {
      param->weil = 0;
      for (k = 0; k < r; k++)
         /* Initialise circular array */
         state->x[k] = (v ^= (v ^= v << 7) >> 9);
   }
   for (k = r; k < (int) state->r2; k++)
       state->x[k] = state->x[k - r];

   for (i = r - 1, k = 4 * r; k > 0; k--) { /* Discard first 4*r results
                                               (Gimeno) */
      t = state->x[i = (i + 1) & param->Mask];
      v = state->x[(i + (r - s)) & param->Mask];
      t ^= (t ^= t << a) >> b;
      v ^= v << c;
      state->x[i] = (v ^= t ^ (v >> d));
   }

   state->i = i;
   param->a = a;
   param->b = b;
   param->c = c;
   param->d = d;
   param->r = r;
   param->s = s;

   return gen;
}


/*-----------------------------------------------------------------------*/

void ubrent_DeleteXorgen64 (unif01_Gen *gen)
{
   Xorgen64_state *state;
   if (NULL == gen)
      return;
   state = gen->state;
   util_Free (state->x);
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/*========================================================================*/
/* 
double xor4096d (unsigned long long seed)

    64-bit real random number generator with period at least 2^4096-1.

    Is is assumed that "unsigned long long" is a 64-bit integer and "double" 
      is an IEEE standard 64-bit floating-point number with 52 explicit + 1
      implicit bits in the fraction.

      The method used is as for the 64-bit integer RNG xor4096l, then the high 
      53 bits (if nonzero) are scaled to (0.0, 1.0).

      Should be called once with nonzero seed, thereafter with zero seed.

      The result should uniformly distributed in (0.0, 1.0) to the resolution
      of the floating-point system.  The result is never exactly 0.0 or 1.0.

      R. P. Brent, 20040802.

*/

#define wlen 64
#define r 64
#define s 53
#define a 33
#define b 26
#define c 27
#define d 29

static double xor4096d_U01 (void *junk, void *vsta)
{
   static ulonglong w, x[r],
      weil = ((longlong) 0x61c88646 << 32) + (longlong) 0x80b583eb;
   ulonglong t, v;
   static double t53 = (double) 1 / (double)((longlong) 1 << 53); /* 0.5**53*/
   static int i = -1;             /* i < 0 indicates first call */
   int k;
   if (i < 0) {  /* Initialisation necessary */
      int seed = *((ulonglong *) vsta);
      v = (seed != 0) ? seed : ~seed; /* v must be nonzero */

      for (k = wlen; k > 0; k--)  /* Avoid correlations for close seeds */
         v ^= (v ^= v << 7) >> 9; /* This recurrence has period 2^64-1 */

      for (w = v, k = 0; k < r; k++)
         /* Initialise circular array */
         x[k] = (v ^= (v ^= v << 7) >> 9) + (w += weil);

      for (i = r - 1, k = 4 * r; k > 0; k--) { /* Discard first 4*r results
                                                  (Gimeno) */
         t = x[i = (i + 1) & (r - 1)];
         v = x[(i + (r - s)) & (r - 1)];
         t ^= (t ^= t << a) >> b;
         v ^= v << c;
         x[i] = (v ^= t ^ (v >> d));
      }
   }

   /* Apart from initialisation (above), this is the generator */

   v = 0;                         /* Usually execute while loop once */
   while (v == (longlong) 0) {    /* Loop until result nonzero */
      t = x[i = (i + 1) & (r - 1)]; /* Increment i mod r (r = power of 2) */
      v = x[(i + (r - s)) & (r - 1)]; /* Index is (i - s) mod r */
      t ^= (t ^= t << a) >> b;    /* (I + L^a)(I + R^b) */
      v ^= v << c;                /* I + L^c */
      x[i] = (v ^= t ^ (v >> d)); /* Update circular array */
      v += (w += weil);           /* 64-bit unsigned integer */
      v >>= 11;                   /* 53-bit integer (possibly zero) */
   }
   return t53 * (double) v;       /* Scale to (0.0, 1.0) */
   }

#undef wlen
#undef r
#undef s
#undef a
#undef b
#undef c
#undef d

/*-----------------------------------------------------------------------*/

static unsigned long xor4096d_Bits (void *vpar, void *vsta)
{
   return xor4096d_U01 (vpar, vsta) * unif01_NORM32;
}

/*-----------------------------------------------------------------------*/

unif01_Gen * ubrent_CreateXor4096d (ulonglong seed)
{
   unif01_Gen *gen;
   ulonglong *pseed;
   size_t leng;
   char name[LEN + 1];

   util_Assert (co3 == 0,
   "ubrent_CreateXor4096d:\n   only 1 such generator can be used at a time");
   co3++;

   gen = util_Malloc (sizeof (unif01_Gen));
   pseed = util_Malloc (sizeof (ulonglong));
   *pseed = seed;
   strcpy (name, "ubrent_CreateXor4096d:");
   addstr_ULONG (name, "   seed = ", seed);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &xor4096d_Bits;
   gen->GetU01 = &xor4096d_U01;
   gen->Write = &Wr_xor4096;
   gen->param = NULL;
   gen->state = pseed;
   /* Discard the first value
   xor4096d_Bits (NULL, pseed); */
   return gen;
}


/*-----------------------------------------------------------------------*/

void ubrent_DeleteXor4096d (unif01_Gen *gen)
{
   if (NULL == gen)
      return;
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
   co3--;
}


/*========================================================================*/
#endif


/**************************************************************************/
/********** HERE BEGIN BRENT'S RNGs FROM THE 2006 VERSION *****************/


typedef unsigned long UINT;       /* Type for random 32 or 64-bit integer,
                                     e.g. unsigned long, unsigned long long,
                                     uint64_t, unsigned int or uint32_t */

typedef double UREAL;             /* Type for random 32 or 64-bit real, e.g.
                                     double or float */

/* xorgens.c version 3.04, R. P. Brent, 20060628. */

/* For type definitions see xorgens.h */

static UINT xor4096i (UINT seed)
{
   /* 32-bit or 64-bit integer random number generator with period at least
      2**4096-1.

      It is assumed that "UINT" is a 32-bit or 64-bit integer (see typedef
      statements in xorgens.h).

      xor4096i should be called exactly once with nonzero seed, and thereafter 
      with zero seed.

      One random number uniformly distributed in [0..2**wlen) is returned,
      where wlen = 8*sizeof(UINT) = 32 or 64.

      R. P. Brent, 20060628. */

   /* UINT64 is TRUE if 64-bit UINT, UINT32 is TRUE otherwise (assumed to be
      32-bit UINT). */

#define UINT64 (sizeof(UINT)>>3)
#define UINT32 (1 - UINT64)

#define wlen (64*UINT64 +  32*UINT32)
#define r    (64*UINT64 + 128*UINT32)
#define s    (53*UINT64 +  95*UINT32)
#define a    (33*UINT64 +  17*UINT32)
#define b    (26*UINT64 +  12*UINT32)
#define c    (27*UINT64 +  13*UINT32)
#define d    (29*UINT64 +  15*UINT32)
#define ws   (27*UINT64 +  16*UINT32)

   static UINT w, weyl, zero = 0, x[r];
   UINT t, v;
   static int i = -1;             /* i < 0 indicates first call */
   int k;
   if ((i < 0) || (seed != zero)) { /* Initialisation necessary */

      /* weyl = odd approximation to 2**wlen*(sqrt(5)-1)/2. */

      if (UINT32)
         weyl = 0x61c88647;
      else
         weyl = ((((UINT) 0x61c88646) << 16) << 16) + (UINT) 0x80b583eb;

      v = (seed != zero) ? seed : ~seed; /* v must be nonzero */

      for (k = wlen; k > 0; k--) { /* Avoid correlations for close seeds */
         v ^= v << 10;
         v ^= v >> 15;            /* Recurrence has period 2**wlen-1 */
         v ^= v << 4;
         v ^= v >> 13;            /* for wlen = 32 or 64 */
      }
      for (w = v, k = 0; k < r; k++) { /* Initialise circular array */
         v ^= v << 10;
         v ^= v >> 15;
         v ^= v << 4;
         v ^= v >> 13;
         x[k] = v + (w += weyl);
      }
      for (i = r - 1, k = 4 * r; k > 0; k--) { /* Discard first 4*r results */
         t = x[i = (i + 1) & (r - 1)];
         t ^= t << a;
         t ^= t >> b;
         v = x[(i + (r - s)) & (r - 1)];
         v ^= v << c;
         v ^= v >> d;
         x[i] = t ^ v;
      }
   }

   /* Apart from initialisation (above), this is the generator */

   t = x[i = (i + 1) & (r - 1)];  /* Assumes that r is a power of two */
   v = x[(i + (r - s)) & (r - 1)]; /* Index is (i-s) mod r */
   t ^= t << a;
   t ^= t >> b;                   /* (I + L^a)(I + R^b) */
   v ^= v << c;
   v ^= v >> d;                   /* (I + L^c)(I + R^d) */
   x[i] = (v ^= t);               /* Update circular array */
   w += weyl;                     /* Update Weyl generator */
   return (v + (w ^ (w >> ws)));  /* Return combination */

#undef UINT64
#undef UINT32
#undef wlen
#undef r
#undef s
#undef a
#undef b
#undef c
#undef d
#undef ws
}


/*-----------------------------------------------------------------------*/

static UINT xor4096i_Bits (void *junk, void *vsta)
{
   UINT seed = *((UINT *) vsta);
   UINT Z = xor4096i(seed);
   *((UINT *) vsta) = 0;
#if ULONG_MAX <= 4294967295UL
   return Z;
#else
   return Z >> 32;
#endif
}


static double xor4096i_U01 (void *vpar, void *vsta)
{
   return xor4096i_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

unif01_Gen * ubrent_CreateXor4096i (unsigned long seed)
{
   unif01_Gen *gen;
   unsigned long *pseed;
   size_t leng;
   char name[LEN + 1];

   util_Assert (co4 == 0,
      "ubrent_CreateXor4096i:\n   only 1 such generator can be used at a time");
   co4++;

   gen = util_Malloc (sizeof (unif01_Gen));
   pseed = util_Malloc (sizeof (unsigned long));
   *pseed = seed;
   strcpy (name, "ubrent_CreateXor4096i:");
   addstr_Ulong (name, "   seed = ", seed);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &xor4096i_Bits;
   gen->GetU01 = &xor4096i_U01;
   gen->Write = &Wr_xor4096;
   gen->param = NULL;
   gen->state = pseed;
   return gen;
}


/*-----------------------------------------------------------------------*/

void ubrent_DeleteXor4096i (unif01_Gen *gen)
{
   if (NULL == gen)
      return;
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
   co4--;
}


/*========================================================================*/

static UREAL xor4096r (UINT seed)
{
   /* 64-bit or 32-bit real random number generator with period at least
      2**4096-1.

      It is assumed that "UINT" is a 32-bit or 64-bit integer and "UREAL" is
      "double" or "float". If "double" this is an IEEE standard floating-point 
      number with 53 bits in the fraction; if "single" it has 24 bits in the
      fraction (including 1 implicit bit in each case).

      In the 64-bit integer case, the method used is to call xor4096i to get
      64 random bits, then the high 53 (for double) or 24 bits (for float) are 
      scaled to the open interval (0.0, 1.0), except that they are discarded
      if all zero.

      In the 32-bit integer case, one or two calls to xor4096i are made to get 
      32 or 64 random bits, some are discarded, and the remaining bits (if
      nonzero) are scaled to the open interval (0.0, 1.0).

      xor4096r should be called exactly once with nonzero seed, and thereafter 
      with zero seed.

      One random number of type UREAL is returned per call.

      The results be should uniformly distributed in (0.0, 1.0) to the
      resolution of the floating-point system (0.5**53 or 0.5**24).

      The results are never exactly 0.0 or 1.0.

      R. P. Brent, 20060628. */

#define UINT64 (sizeof(UINT)>>3)
#define UINT32 (1 - UINT64)
#define UREAL64 (sizeof(UREAL)>>3)
#define UREAL32 (1 - UREAL64)

/* sr = number of bits discarded = 11 for double, 40 or 8 for float */

#define sr (11*UREAL64 +(40*UINT64 + 8*UINT32)*UREAL32)

/* ss (used for scaling) is 53 or 21 for double, 24 for float */

#define ss ((53*UINT64 + 21*UINT32)*UREAL64 + 24*UREAL32)

/* SCALE is 0.5**ss, SC32 is 0.5**32 */

#define SCALE ((UREAL)1/(UREAL)((UINT)1<<ss))
#define SC32  ((UREAL)1/((UREAL)65536*(UREAL)65536))

   UREAL res;

   res = (UREAL) 0;
   while (res == (UREAL) 0) {     /* Loop until nonzero result.  *//* Usually only one iteration . */
      res = (UREAL) (xor4096i (seed) >> sr); /* Discard sr random bits.  */
      seed = (UINT) 0;            /* Zero seed for next time. */
      if (UINT32 && UREAL64)      /* Need another call to xor4096i. */
         res += SC32 * (UREAL) xor4096i (seed); /* Add low-order 32 bits. */
   }
   return (SCALE * res);          /* Return result in (0.0, 1.0). */

#undef UINT64
#undef UINT32
#undef UREAL64
#undef UREAL32
#undef SCALE
#undef SC32
#undef sr
#undef ss
}


/*-----------------------------------------------------------------------*/

static double xor4096r_U01 (void *junk, void *vsta)
{
   UINT seed = *((UINT *) vsta);
   double u = xor4096r(seed);
   *((UINT *) vsta) = 0;
   return u;
}


static unsigned long xor4096r_Bits (void *vpar, void *vsta)
{
   return xor4096r_U01 (vpar, vsta) * unif01_NORM32;
}

/*-----------------------------------------------------------------------*/

unif01_Gen * ubrent_CreateXor4096r (unsigned long seed)
{
   unif01_Gen *gen;
   unsigned long *pseed;
   size_t leng;
   char name[LEN + 1];

   util_Assert (co5 == 0,
   "ubrent_CreateXor4096r:\n   only 1 such generator can be used at a time");
   co5++;

   gen = util_Malloc (sizeof (unif01_Gen));
   pseed = util_Malloc (sizeof (unsigned long));
   *pseed = seed;
   strcpy (name, "ubrent_CreateXor4096r:");
   addstr_ULONG (name, "   seed = ", seed);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &xor4096r_Bits;
   gen->GetU01 = &xor4096r_U01;
   gen->Write = &Wr_xor4096;
   gen->param = NULL;
   gen->state = pseed;
   return gen;
}


/*-----------------------------------------------------------------------*/

void ubrent_DeleteXor4096r (unif01_Gen *gen)
{
   if (NULL == gen)
      return;
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
   co5--;
}
