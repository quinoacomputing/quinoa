/*************************************************************************\
 *
 * Package:        TestU01
 * File:           utouzin.c
 * Environment:    ANSI C
 * Programmer:     R. Touzin
 * Adapted to TestU01 by R. Simard
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
#include "utouzin.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>


typedef struct {
   unsigned long x1, x2, x3, x4, x5, x6, x7, x8;
} MRG00a_state;

typedef MRG00a_state MRG00b_state;
typedef MRG00a_state MRG00c_state;
typedef MRG00a_state MRG00d_state;


typedef struct {
   unsigned long x10, x11, x12, x20, x21, x22, x13, x23;
} MRG00e_state;

typedef MRG00e_state MRG00f_state;
typedef MRG00e_state MRG00h_state;

typedef struct {
   unsigned long x10, x11, x12, x20, x21, x22, x30, x31, x32;
} MRG00g_state;

#define  LEN 200




/*=========================================================================*/

#define MASK1a  127               /* 2^(31 - 24) - 1 */
#define MASK3a  8191              /* 2^(31 - 18) - 1 */
#define MASK4a  134217727         /* 2^(31 - 4) - 1 */
#define MASK5a  1048575           /* 2^(31 - 11) - 1 */
#define m1      2147483647        /* 2^31 - 1 */
#define norm    4.656612875245797e-10 /* 1/m1 */

static double MRG00a_U01 (void *junk, void *vsta)
{
   MRG00a_state *state = vsta;
   unsigned long Resultat;

   /* Coefficient 1 */
   Resultat = m1 - (((state->x1 & MASK1a) << 24) + (state->x1 >> 7));
   Resultat += state->x1;
   if (Resultat >= m1)
      Resultat -= m1;
   Resultat += state->x1;
   if (Resultat >= m1)
      Resultat -= m1;

   /* Coefficient 3 */
   Resultat += (m1 - (((state->x3 & MASK3a) << 18) + (state->x3 >> 13)));
   if (Resultat >= m1)
      Resultat -= m1;

   /* Coefficient 4 */
   Resultat += (m1 - (((state->x4 & MASK4a) << 4) + (state->x4 >> 27)));
   if (Resultat >= m1)
      Resultat -= m1;

   /* Coefficient 5 */
   Resultat += (((state->x5 & MASK5a) << 11) + (state->x5 >> 20));
   if (Resultat >= m1)
      Resultat -= m1;
   Resultat = Resultat + (m1 - state->x5);
   if (Resultat >= m1)
      Resultat -= m1;

   state->x5 = state->x4;
   state->x4 = state->x3;
   state->x3 = state->x2;
   state->x2 = state->x1;

   state->x1 = Resultat;
   return (norm * state->x1);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG00a_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG00a_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG00a (void *vsta)
{
   MRG00a_state *state = vsta;
   printf
      ("  x1 = %10lu,   x2 = %10lu,   x3 = %10lu,\n  x4 = %10lu,   x5 = %10lu\n\n",
      state->x1, state->x2, state->x3, state->x4, state->x5);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * utouzin_CreateMRG00a (long s1, long s2, long s3, long s4,
   long s5)
{
   unif01_Gen *gen;
   MRG00a_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG00a_state));

   strncpy (name, "utouzin_CreateMRG00a:", LEN);
   addstr_Long (name, "   s1 = ", s1);
   addstr_Long (name, ",  s2 = ", s2);
   addstr_Long (name, ",  s3 = ", s3);
   addstr_Long (name, ",  s4 = ", s4);
   addstr_Long (name, ",  s5 = ", s5);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->x1 = s1 % m1;
   state->x2 = s2 % m1;
   state->x3 = s3 % m1;
   state->x4 = s4 % m1;
   state->x5 = s5 % m1;

   gen->GetBits = &MRG00a_Bits;
   gen->GetU01 = &MRG00a_U01;
   gen->Write = &WrMRG00a;
   gen->param = NULL;
   gen->state = state;
   return gen;
}

#undef MASK1a
#undef MASK3a
#undef MASK4a
#undef MASK5a
#undef m1
#undef norm


/*=========================================================================*/

#define MASK1a  1023              /* 2^(31 - 21) -1 */
#define MASK2a  524287            /* 2^(31 - 12) -1 */
#define MASK3a  32767             /* 2^(31 - 16) -1 */
#define MASK5a  16777215          /* 2^(31 - 7) -1 */
#define MASK6a  15                /* 2^(31 - 27) -1 */
#define m1      2147483647        /* 2^31 -1 */
#define norm    4.656612875245797e-10 /* 1/m1 */

static double MRG00b_U01 (void *junk, void *vsta)
{
   MRG00b_state *state = vsta;
   unsigned long Resultat;

   /* Coefficient 1 */
   Resultat = m1 - (((state->x1 & MASK1a) << 21) + (state->x1 >> 10));
   Resultat += (m1 - state->x1);
   if (Resultat >= m1)
      Resultat -= m1;

   /* Coefficient 2 */
   Resultat += (m1 - (((state->x2 & MASK2a) << 12) + (state->x2 >> 19)));
   if (Resultat >= m1)
      Resultat -= m1;

   /* Coefficient 3 */
   Resultat += (((state->x3 & MASK3a) << 16) + (state->x3 >> 15));
   if (Resultat >= m1)
      Resultat -= m1;

   /* Coefficient 5 */
   Resultat += (((state->x5 & MASK5a) << 7) + (state->x5 >> 24));
   if (Resultat >= m1)
      Resultat -= m1;

   /* Coefficient 6 */
   Resultat += (m1 - (((state->x6 & MASK6a) << 27) + (state->x6 >> 4)));
   if (Resultat >= m1)
      Resultat -= m1;
   Resultat += state->x6;
   if (Resultat >= m1)
      Resultat -= m1;

   state->x6 = state->x5;
   state->x5 = state->x4;
   state->x4 = state->x3;
   state->x3 = state->x2;
   state->x2 = state->x1;

   state->x1 = Resultat;
   return (norm * state->x1);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG00b_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG00b_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG00b (void *vsta)
{
   MRG00b_state *state = vsta;
   printf
      ("  x1 = %10lu,   x2 = %10lu,   x3 = %10lu,\n  x4 = %10lu,   x5 = %10lu,   x6 = %10lu\n\n",
      state->x1, state->x2, state->x3, state->x4, state->x5, state->x6);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * utouzin_CreateMRG00b (long s1, long s2, long s3,
   long s4, long s5, long s6)
{
   unif01_Gen *gen;
   MRG00b_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG00b_state));

   strncpy (name, "utouzin_CreateMRG00b:", LEN);
   addstr_Long (name, "   s1 = ", s1);
   addstr_Long (name, ",  s2 = ", s2);
   addstr_Long (name, ",  s3 = ", s3);
   addstr_Long (name, ",  s4 = ", s4);
   addstr_Long (name, ",  s5 = ", s5);
   addstr_Long (name, ",  s6 = ", s6);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->x1 = s1 % m1;
   state->x2 = s2 % m1;
   state->x3 = s3 % m1;
   state->x4 = s4 % m1;
   state->x5 = s5 % m1;
   state->x6 = s6 % m1;

   gen->GetBits = &MRG00b_Bits;
   gen->GetU01 = &MRG00b_U01;
   gen->Write = &WrMRG00b;
   gen->param = NULL;
   gen->state = state;
   return gen;
}

#undef MASK1a
#undef MASK2a
#undef MASK3a
#undef MASK5a
#undef MASK6a
#undef m1
#undef norm


/*=========================================================================*/

#define MASK1a  524287            /* 2^(31 - 12) -1 */
#define MASK2a  2047              /* 2^(31 - 20) -1 */
#define MASK3a  131071            /* 2^(31 - 14) -1 */
#define MASK5a  63                /* 2^(31 - 25) -1 */
#define MASK6a  33554431          /* 2^(31 - 6) -1 */
#define MASK7a  134217727         /* 2^(31-4)-1 */
#define m1      2147483629        /* 2^31 -19 */
#define norm    4.656612914277075e-10  /* 1/m1 */

static double MRG00c_U01 (void *junk, void *vsta)
{
   MRG00c_state *state = vsta;
   unsigned long Resultat1, Resultat2;

   /* Coefficient 1 */
   Resultat1 = (((state->x1 & MASK1a) << 12) + 19 * (state->x1 >> 19));
   if (Resultat1 >= m1)
      Resultat1 = m1 - (Resultat1 - m1);
   else
      Resultat1 = m1 - Resultat1;

   /* Coefficient 2 */
   Resultat2 = (((state->x2 & MASK2a) << 20) + 19 * (state->x2 >> 11));
   if (Resultat2 >= m1)
      Resultat2 = m1 - (Resultat2 - m1);
   else
      Resultat2 = m1 - Resultat2;

   Resultat2 += Resultat1;
   if (Resultat2 >= m1)
      Resultat2 -= m1;

   /* Coefficient 3 */
   Resultat1 = (((state->x3 & MASK3a) << 14) + 19 * (state->x3 >> 17));
   if (Resultat1 >= m1)
      Resultat1 -= m1;

   Resultat2 += Resultat1;
   if (Resultat2 >= m1)
      Resultat2 -= m1;

   /* Coefficient 5 */
   Resultat1 = (((state->x5 & MASK5a) << 25) + 19 * (state->x5 >> 6));
   if (Resultat1 >= m1)
      Resultat1 -= m1;

   Resultat2 += Resultat1;
   if (Resultat2 >= m1)
      Resultat2 -= m1;

   /* Coefficient 6 */
   Resultat1 = (((state->x6 & MASK6a) << 6) + 19 * (state->x6 >> 25));
   if (Resultat1 >= m1)
      Resultat1 = m1 - (Resultat1 - m1);
   else
      Resultat1 = m1 - Resultat1;

   Resultat2 += Resultat1;
   if (Resultat2 >= m1)
      Resultat2 -= m1;

   /* Coefficient 7 */
   Resultat1 = (((state->x7 & MASK7a) << 4) + 19 * (state->x7 >> 27));
   if (Resultat1 >= m1)
      Resultat1 -= m1;
   Resultat1 += state->x7;
   if (Resultat1 >= m1)
      Resultat1 -= m1;

   Resultat2 += Resultat1;
   if (Resultat2 >= m1)
      Resultat2 -= m1;

   state->x7 = state->x6;
   state->x6 = state->x5;
   state->x5 = state->x4;
   state->x4 = state->x3;
   state->x3 = state->x2;
   state->x2 = state->x1;

   state->x1 = Resultat2;
   return (norm * state->x1);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG00c_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG00c_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG00c (void *vsta)
{
   MRG00c_state *state = vsta;
   printf
      ("  x1 = %10lu,   x2 = %10lu,   x3 = %10lu,   x4 = %10lu,\n  x5 = %10lu,   x6 = %10lu,   x7 = %10lu\n\n",
      state->x1, state->x2, state->x3, state->x4, state->x5,
      state->x6, state->x7);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * utouzin_CreateMRG00c (long s1, long s2, long s3,
   long s4, long s5, long s6, long s7)
{
   unif01_Gen *gen;
   MRG00c_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG00c_state));

   strncpy (name, "utouzin_CreateMRG00c:", LEN);
   addstr_Long (name, "   s1 = ", s1);
   addstr_Long (name, ",  s2 = ", s2);
   addstr_Long (name, ",  s3 = ", s3);
   addstr_Long (name, ",  s4 = ", s4);
   addstr_Long (name, ",  s5 = ", s5);
   addstr_Long (name, ",  s6 = ", s6);
   addstr_Long (name, ",  s7 = ", s7);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->x1 = s1 % m1;
   state->x2 = s2 % m1;
   state->x3 = s3 % m1;
   state->x4 = s4 % m1;
   state->x5 = s5 % m1;
   state->x6 = s6 % m1;
   state->x7 = s7 % m1;

   gen->GetBits = &MRG00c_Bits;
   gen->GetU01 = &MRG00c_U01;
   gen->Write = &WrMRG00c;
   gen->param = NULL;
   gen->state = state;
   return gen;
}

#undef MASK1a
#undef MASK2a
#undef MASK3a
#undef MASK5a
#undef MASK6a
#undef MASK7a
#undef m1
#undef norm


/*=========================================================================*/

#define MASK1a  134217727         /* 2^(31 - 4) -1 */
#define MASK3a  65535             /* 2^(31 - 15) -1 */
#define MASK4a  524287            /* 2^(31 - 12) -1 */
#define MASK5a  511               /* 2^(31 - 22) -1 */
#define MASK6a  4194303           /* 2^(31 - 9) -1 */
#define MASK7a  15                /* 2^(31 - 27) -1 */
#define MASK8a  8191              /* 2^(31 - 18) -1 */
#define MASK8b  1073741823        /* 2^(31 - 1) -1 */
#define m1      2147483647        /* 2^31 -1 */
#define norm    4.656612875245797e-10 /* 1/m1 */

static double MRG00d_U01 (void *junk, void *vsta)
{
   MRG00d_state *state = vsta;
   unsigned long Resultat;

   /* Coefficient 1 */
   Resultat = m1 - (((state->x1 & MASK1a) << 4) + (state->x1 >> 27));

   /* Coefficient 3 */
   Resultat += (((state->x3 & MASK3a) << 15) + (state->x3 >> 16));
   if (Resultat >= m1)
      Resultat -= m1;

   /* Coefficient 4 */
   Resultat += (m1 - (((state->x4 & MASK4a) << 12) + (state->x4 >> 19)));
   if (Resultat >= m1)
      Resultat -= m1;

   /* Coefficient 5 */
   Resultat += (((state->x5 & MASK5a) << 22) + (state->x5 >> 9));
   if (Resultat >= m1)
      Resultat -= m1;

   /* Coefficient 6 */
   Resultat += (((state->x6 & MASK6a) << 9) + (state->x6 >> 22));
   if (Resultat >= m1)
      Resultat -= m1;

   /* Coefficient 7 */
   Resultat += (((state->x7 & MASK7a) << 27) + (state->x7 >> 4));
   if (Resultat >= m1)
      Resultat -= m1;

   /* Coefficient 8 */
   Resultat += (((state->x8 & MASK8a) << 18) + (state->x8 >> 13));
   if (Resultat >= m1)
      Resultat -= m1;

   Resultat += (m1 - state->x8);
   if (Resultat >= m1)
      Resultat -= m1;
   Resultat += (m1 - state->x8);
   if (Resultat >= m1)
      Resultat -= m1;

   state->x8 = state->x7;
   state->x7 = state->x6;
   state->x6 = state->x5;
   state->x5 = state->x4;
   state->x4 = state->x3;
   state->x3 = state->x2;
   state->x2 = state->x1;
   state->x1 = Resultat;
   return (norm * state->x1);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG00d_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG00d_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG00d (void *vsta)
{
   MRG00d_state *state = vsta;
   printf
      ("  x1 = %10lu,   x2 = %10lu,   x3 = %10lu,   x4 = %10lu,\n  x5 = %10lu,   x6 = %10lu,   x7 = %10lu,   x8 = %10lu\n\n",
      state->x1, state->x2, state->x3, state->x4, state->x5,
      state->x6, state->x7, state->x8);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * utouzin_CreateMRG00d (long s1, long s2, long s3,
   long s4, long s5, long s6, long s7, long s8)
{
   unif01_Gen *gen;
   MRG00d_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG00d_state));

   strncpy (name, "utouzin_CreateMRG00d:", LEN);
   addstr_Long (name, "   s1 = ", s1);
   addstr_Long (name, ",  s2 = ", s2);
   addstr_Long (name, ",  s3 = ", s3);
   addstr_Long (name, ",  s4 = ", s4);
   addstr_Long (name, ",  s5 = ", s5);
   addstr_Long (name, ",  s6 = ", s6);
   addstr_Long (name, ",  s7 = ", s7);
   addstr_Long (name, ",  s8 = ", s8);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->x1 = s1 % m1;
   state->x2 = s2 % m1;
   state->x3 = s3 % m1;
   state->x4 = s4 % m1;
   state->x5 = s5 % m1;
   state->x6 = s6 % m1;
   state->x7 = s7 % m1;
   state->x8 = s8 % m1;

   gen->GetBits = &MRG00d_Bits;
   gen->GetU01 = &MRG00d_U01;
   gen->Write = &WrMRG00d;
   gen->param = NULL;
   gen->state = state;
   return gen;
}

#undef MASK1a
#undef MASK3a
#undef MASK4a
#undef MASK5a
#undef MASK6a
#undef MASK7a
#undef MASK8a
#undef MASK8b
#undef m1
#undef norm


/*=========================================================================*/

#define MASK11  511               /* 2^(31-22) - 1 */
#define MASK12  16777215          /* 2^(31-7) - 1 */
#define MASK20  65535             /* 2^(31-15) - 1 */
#define MASK22  65535             /* 2^(31-15) - 1 */
#define m1      2147483647        /* 2^31 -1 */
#define m2      2147462579        /* 2^31 - 21069 */
#define norm    4.656612873077393e-10  /* 1/(m1 + 1) */

static double MRG00e_U01 (void *junk, void *vsta)
{
   MRG00e_state *state = vsta;
   unsigned long Resultat1;
   unsigned long Resultat2;

   /* Component 1 x10 = ((2^22)* x11 + (2^7 + 1)*x12) mod 2^31 -1 */
   Resultat1 = (((state->x11 & MASK11) << 22) + (state->x11 >> 9))
      + (((state->x12 & MASK12) << 7) + (state->x12 >> 24));
   if (Resultat1 >= m1)
      Resultat1 -= m1;
   Resultat1 += state->x12;
   if (Resultat1 >= m1)
      Resultat1 -= m1;
   state->x12 = state->x11;
   state->x11 = state->x10;
   state->x10 = Resultat1;

   /* Component 2 x20 = ((2^15)* x20 + (2^15+1) * x22) mod (2^31 - 21069) */
   Resultat1 = ((state->x20 & MASK20) << 15) + 21069 * (state->x20 >> 16);
   if (Resultat1 >= m2)
      Resultat1 -= m2;
   Resultat2 = ((state->x22 & MASK22) << 15) + 21069 * (state->x22 >> 16);
   if (Resultat2 >= m2)
      Resultat2 -= m2;
   Resultat1 += Resultat2;
   if (Resultat1 >= m2)
      Resultat1 -= m2;
   Resultat1 += state->x22;
   if (Resultat1 >= m2)
      Resultat1 -= m2;
   state->x22 = state->x21;
   state->x21 = state->x20;
   state->x20 = Resultat1;

   /* Combinaison */
   if (state->x10 <= state->x20)
      return ((state->x10 - state->x20 + m1) * norm);
   else
      return ((state->x10 - state->x20) * norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG00e_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG00e_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG00e (void *vsta)
{
   MRG00e_state *state = vsta;
   printf
      ("  x10 = %10lu,   x11 = %10lu,   x12 = %10lu,\n  x20 = %10lu,   x21 = %10lu,   x22 = %10lu \n\n",
      state->x10, state->x11, state->x12, state->x20, state->x21, state->x22);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * utouzin_CreateMRG00e (long s10, long s11, long s12,
   long s20, long s21, long s22)
{
   unif01_Gen *gen;
   MRG00e_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG00e_state));

   strncpy (name, "utouzin_CreateMRG00e:", LEN);
   addstr_Long (name, "  s10 = ", s10);
   addstr_Long (name, ",  s11 = ", s11);
   addstr_Long (name, ",  s12 = ", s12);
   addstr_Long (name, ",  s20 = ", s20);
   addstr_Long (name, ",  s21 = ", s21);
   addstr_Long (name, ",  s22 = ", s22);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->x10 = s10;
   state->x11 = s11;
   state->x12 = s12;
   state->x20 = s20;
   state->x21 = s21;
   state->x22 = s22;

   gen->GetBits = &MRG00e_Bits;
   gen->GetU01 = &MRG00e_U01;
   gen->Write = &WrMRG00e;
   gen->param = NULL;
   gen->state = state;
   return gen;
}

#undef MASK11
#undef MASK12
#undef MASK20
#undef MASK22
#undef m1
#undef m2
#undef norm


/*=========================================================================*/

#define MASK11  131071            /* 2^(31-14) - 1 */
#define MASK12  31                /* 2^(31-26) - 1 */
#define MASK20  16383             /* 2^(31-17) - 1 */
#define MASK22  1048575           /* 2^(31-11) - 1 */
#define m1      2147483647        /* 2^31 -1 */
#define m2      2147483629        /* 2^31 - 21069 */
#define norm    4.656612873077393e-10  /* 1/(m1 + 1) */

static double MRG00f_U01 (void *junk, void *vsta)
{
   MRG00f_state *state = vsta;
   unsigned long Resultat1;
   unsigned long Resultat2;

   /* Component 1 x10 = ((2^14)* x11 + (-2^26 + 1)*x12) mod 2^31 -1 */
   Resultat1 = (((state->x11 & MASK11) << 14) + (state->x11 >> 17))
      + (m1 - (((state->x12 & MASK12) << 26) + (state->x12 >> 5)));
   if (Resultat1 >= m1)
      Resultat1 -= m1;
   Resultat1 += state->x12;
   if (Resultat1 >= m1)
      Resultat1 -= m1;
   state->x12 = state->x11;
   state->x11 = state->x10;
   state->x10 = Resultat1;

   /* Component 2 x20 = ((2^17)* x20 + (2^11) * x22) mod (2^31 - 19) */
   Resultat1 = ((state->x20 & MASK20) << 17) + (19 * (state->x20 >> 14));
   if (Resultat1 >= m2)
      Resultat1 -= m2;
   Resultat2 = ((state->x22 & MASK22) << 11) + (19 * (state->x22 >> 20));
   if (Resultat2 >= m2)
      Resultat2 -= m2;
   Resultat1 += Resultat2;
   if (Resultat1 >= m2)
      Resultat1 -= m2;

   state->x22 = state->x21;
   state->x21 = state->x20;
   state->x20 = Resultat1;

   /* Combinaison */
   if (state->x10 <= state->x20)
      return ((state->x10 - state->x20 + m1) * norm);
   else
      return ((state->x10 - state->x20) * norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG00f_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG00f_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG00f (void *vsta)
{
   MRG00f_state *state = vsta;
   printf
      ("  x10 = %10lu,   x11 = %10lu,   x12 = %10lu,\n  x20 = %10lu,   x21 = %10lu,   x22 = %10lu \n\n",
      state->x10, state->x11, state->x12, state->x20, state->x21, state->x22);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * utouzin_CreateMRG00f (long s10, long s11, long s12,
   long s20, long s21, long s22)
{
   unif01_Gen *gen;
   MRG00f_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG00f_state));

   strncpy (name, "utouzin_CreateMRG00f:", LEN);
   addstr_Long (name, "  s10 = ", s10);
   addstr_Long (name, ",  s11 = ", s11);
   addstr_Long (name, ",  s12 = ", s12);
   addstr_Long (name, ",  s20 = ", s20);
   addstr_Long (name, ",  s21 = ", s21);
   addstr_Long (name, ",  s22 = ", s22);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->x10 = s10;
   state->x11 = s11;
   state->x12 = s12;
   state->x20 = s20;
   state->x21 = s21;
   state->x22 = s22;

   gen->GetBits = &MRG00f_Bits;
   gen->GetU01 = &MRG00f_U01;
   gen->Write = &WrMRG00f;
   gen->param = NULL;
   gen->state = state;
   return gen;
}

#undef MASK11
#undef MASK12
#undef MASK20
#undef MASK22
#undef m1
#undef m2
#undef norm


/*=========================================================================*/

#define MASK10  1                 /* 2^(31-30) - 1 */
#define MASK12  4095              /* 2^(31-19) - 1 */
#define MASK21  255               /* 2^(31-23) - 1 */
#define MASK22  4095              /* 2^(31-19) - 1 */
#define MASK30  1048575           /* 2^(31-11) - 1 */
#define MASK31  4194303           /* 2^(31-9 ) - 1 */
#define m1      2147483647        /* 2^31 -1 */
#define m2      2147483629        /* 2^31 -19 */
#define m3      2147483587        /* 2^31 -61 */
#define norm    4.656612873077393e-10  /* 1/(m1 + 1) */

static double MRG00g_U01 (void *junk, void *vsta)
{
   MRG00g_state *state = vsta;
   unsigned long Resultat1;
   unsigned long Resultat2;

   /* Component 1 x10 = ((2^14)* x11 + (-2^26 + 1)*x12) mod 2^31 -1 */
   Resultat1 = (((state->x10 & MASK10) << 30) + (state->x10 >> 1))
      + (((state->x12 & MASK12) << 19) + (state->x12 >> 12));
   if (Resultat1 >= m1)
      Resultat1 -= m1;
   Resultat1 += (m1 - state->x12);
   if (Resultat1 >= m1)
      Resultat1 -= m1;
   state->x12 = state->x11;
   state->x11 = state->x10;
   state->x10 = Resultat1;

   /* Component 2 x20 = ((2^17)* x20 + (2^11) * x22) mod (2^31 - 19) */
   Resultat1 = ((state->x21 & MASK21) << 23) + (19 * (state->x21 >> 8));
   if (Resultat1 >= m2)
      Resultat1 -= m2;
   Resultat2 = ((state->x22 & MASK22) << 19) + (19 * (state->x22 >> 12));
   if (Resultat2 >= m2)
      Resultat2 -= m2;
   Resultat1 += Resultat2;
   if (Resultat1 >= m2)
      Resultat1 -= m2;
   state->x22 = state->x21;
   state->x21 = state->x20;
   state->x20 = Resultat1;

   /* Component 3 x20 = ((2^17)* x20 + (2^11) * x22) mod (2^31 - 19) */
   Resultat1 = ((state->x30 & MASK30) << 11) + (61 * (state->x30 >> 20));
   if (Resultat1 >= m3)
      Resultat1 -= m3;
   Resultat2 = ((state->x31 & MASK31) << 9) + (61 * (state->x31 >> 22));
   if (Resultat2 >= m3)
      Resultat2 -= m3;
   Resultat1 += Resultat2;
   if (Resultat1 >= m3)
      Resultat1 -= m3;
   Resultat1 += state->x32;
   if (Resultat1 >= m3)
      Resultat1 -= m3;
   Resultat1 += state->x32;
   if (Resultat1 >= m3)
      Resultat1 -= m3;
   state->x32 = state->x31;
   state->x31 = state->x30;
   state->x30 = Resultat1;

   /* Combinaison */
   if (state->x10 + state->x30 <= state->x20)
      return ((state->x10 - state->x20 + state->x30 + m1) * norm);
   else if (state->x10 - state->x20 + state->x30 > m1)
      return ((state->x10 - state->x20 + state->x30 - m1) * norm);
   else
      return ((state->x10 - state->x20 + state->x30) * norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG00g_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG00g_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG00g (void *vsta)
{
   MRG00g_state *state = vsta;
   printf
      ("  x10 = %10lu,   x11 = %10lu,   x12 = %10lu,\n  x20 = %10lu,   x21 = %10lu,   x22 = %10lu,\n  x30 = %10lu,   x31 = %10lu,   x32 = %10lu\n\n",
      state->x10, state->x11, state->x12, state->x20, state->x21, state->x22,
      state->x30, state->x31, state->x32);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * utouzin_CreateMRG00g (long s10, long s11, long s12,
   long s20, long s21, long s22,
   long s30, long s31, long s32)
{
   unif01_Gen *gen;
   MRG00g_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG00g_state));

   strncpy (name, "utouzin_CreateMRG00g:", LEN);
   addstr_Long (name, "  s10 = ", s10);
   addstr_Long (name, ",  s11 = ", s11);
   addstr_Long (name, ",  s12 = ", s12);
   addstr_Long (name, ",  s20 = ", s20);
   addstr_Long (name, ",  s21 = ", s21);
   addstr_Long (name, ",  s22 = ", s22);
   addstr_Long (name, ",  s30 = ", s30);
   addstr_Long (name, ",  s31 = ", s31);
   addstr_Long (name, ",  s32 = ", s32);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->x10 = s10;
   state->x11 = s11;
   state->x12 = s12;
   state->x20 = s20;
   state->x21 = s21;
   state->x22 = s22;
   state->x30 = s30;
   state->x31 = s31;
   state->x32 = s32;

   gen->GetBits = &MRG00g_Bits;
   gen->GetU01 = &MRG00g_U01;
   gen->Write = &WrMRG00g;
   gen->param = NULL;
   gen->state = state;
   return gen;
}

#undef MASK10
#undef MASK12
#undef MASK21
#undef MASK22
#undef MASK30
#undef MASK31
#undef m1
#undef m2
#undef m3
#undef norm


/*=========================================================================*/

#define MASK11  262143            /* 2^(31-13) - 1 */
#define MASK13  255               /* 2^(31-23) - 1 */
#define MASK20  2097151           /* 2^(31-10) - 1 */
#define MASK22  2047              /* 2^(31-20) - 1 */
#define MASK23  16777215          /* 2^(31 -7) - 1 */
#define m1      2147483647        /* 2^31 -1 */
#define m2      2147483629        /* 2^31 - 21069 */
#define norm    4.656612873077393e-10  /* 1/(m1 + 1) */

static double MRG00h_U01 (void *junk, void *vsta)
{
   MRG00h_state *state = vsta;
   unsigned long Resultat1;
   unsigned long Resultat2;

   /* Component 1 x10 = ((2^14)* x11 + (-2^26 + 1)*x12) mod 2^31 -1 */
   Resultat1 = (m1 - state->x10) + (m1 - (((state->x11 & MASK11) << 13) +
               (state->x11 >> 18)));
   if (Resultat1 >= m1)
      Resultat1 -= m1;
   Resultat1 += ((state->x13 & MASK13) << 23) + (state->x13 >> 8);
   if (Resultat1 >= m1)
      Resultat1 -= m1;
   Resultat1 += state->x13;
   if (Resultat1 >= m1)
      Resultat1 -= m1;

   state->x13 = state->x12;
   state->x12 = state->x11;
   state->x11 = state->x10;
   state->x10 = Resultat1;

   /* Component 2 x20 = ((2^17)* x20 + (2^11) * x22) mod (2^31 - 19) */
   Resultat1 = ((state->x20 & MASK20) << 10) + (19 * (state->x20 >> 21));
   if (Resultat1 >= m2)
      Resultat1 -= m2;

   Resultat2 = (((state->x22 & MASK22) << 20) + (19 * (state->x22 >> 11)));
   if (Resultat2 >= m2)
      Resultat2 = m2 - (Resultat2 - m2);
   else
      Resultat2 = m2 - Resultat2;
   if (Resultat2 >= m2)
      Resultat2 -= m2;

   Resultat1 += Resultat2;
   if (Resultat1 >= m2)
      Resultat1 -= m2;
   Resultat2 = ((state->x23 & MASK23) << 7) + (19 * (state->x23 >> 24));
   if (Resultat2 >= m2)
      Resultat2 -= m2;
   Resultat1 += Resultat2;
   if (Resultat1 >= m2)
      Resultat1 -= m2;

   state->x23 = state->x22;
   state->x22 = state->x21;
   state->x21 = state->x20;
   state->x20 = Resultat1;

   /* Combinaison */
   if (state->x10 <= state->x20)
      return ((state->x10 - state->x20 + m1) * norm);
   else
      return ((state->x10 - state->x20) * norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG00h_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG00h_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG00h (void *vsta)
{
   MRG00h_state *state = vsta;
   printf
      (" x10 = %10lu,  x11 = %10lu,  x12 = %10lu,  x13 = %10lu,\n x20 = %10lu,  x21 = %10lu,  x22 = %10lu,  x23 = %10lu\n\n",
      state->x10, state->x11, state->x12, state->x13,
      state->x20, state->x21, state->x22, state->x23);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * utouzin_CreateMRG00h (long s10, long s11, long s12,
   long s13, long s20, long s21, long s22, long s23)
{
   unif01_Gen *gen;
   MRG00h_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG00h_state));

   strncpy (name, "utouzin_CreateMRG00h:", LEN);
   addstr_Long (name, "  s10 = ", s10);
   addstr_Long (name, ",  s11 = ", s11);
   addstr_Long (name, ",  s12 = ", s12);
   addstr_Long (name, ",  s13 = ", s13);
   addstr_Long (name, ",  s20 = ", s20);
   addstr_Long (name, ",  s21 = ", s21);
   addstr_Long (name, ",  s22 = ", s22);
   addstr_Long (name, ",  s23 = ", s23);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->x10 = s10;
   state->x11 = s11;
   state->x12 = s12;
   state->x13 = s13;
   state->x20 = s20;
   state->x21 = s21;
   state->x22 = s22;
   state->x23 = s23;
   gen->GetBits = &MRG00h_Bits;
   gen->GetU01 = &MRG00h_U01;
   gen->Write = &WrMRG00h;
   gen->param = NULL;
   gen->state = state;
   return gen;
}


/*=========================================================================*/

void utouzin_DeleteGen (unif01_Gen * gen)
{
   if (NULL == gen)
      return;
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}
