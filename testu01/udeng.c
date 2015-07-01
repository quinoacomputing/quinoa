/*************************************************************************\
 *
 * Package:        TestU01
 * File:           udeng.c
 * Environment:    ANSI C
 * Programmer:     Richard Simard.
 *
\*************************************************************************/

#include "gdef.h"
#include "util.h"
#include "addstr.h"
#include "num.h"

#include "udeng.h"
#include "unif01.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#define LEN  200                  /* Max length of strings */

static int co1 = 0;               /* Counters */


typedef struct {
   double b, m;
} DX02_param;

typedef struct {
   double *X;
   unsigned int s, K; 
} DX02_state;







/*=========================================================================*/
#define MASK1 127
#define LAC

static unif01_Gen *CreateDenga (unsigned long m, unsigned long b, int k,
				unsigned long S[], int id);

#ifndef LAC
static double DX02a_U01 (void *vpar, void *vsta)
{
   DX02_state *state = vsta;
   DX02_param *param = vpar;
   double Z;
   ++state->s;

   Z = param->b * (state->X[(state->s - 1) & MASK1] +
                   state->X[(state->s - state->K) & MASK1]);
   Z = state->X[state->s & MASK1] = fmod (Z, param->m);
   return  Z / param->m;
   /*
   Z = state->X[state->s & MASK1] = num_MultModD (param->b, 
       state->X[(state->s - 1) & MASK1] + state->X[(state->s - state->K) &
                MASK1], 0.0, param->m);
   return  Z / param->m;
   */
}

#else
static double DX02a_U01 (void *vpar, void *vsta)
{
  /* Generate 3 numbers; discard the next 100 numbers; and so on ... */
   DX02_state *state = vsta;
   DX02_param *param = vpar;
   double Z;
   int i;
   static int co = 2;
   ++state->s;

   if (co % 3 == 0) {
      for (i = 0; i < 100; i++) {
         state->X[state->s & MASK1] = num_MultModD (param->b, 
         state->X[(state->s - 1) & MASK1] + state->X[(state->s - state->K) &
            MASK1], 0.0, param->m);
         ++state->s;
      }
      co = 0;
   }
   co++;
   Z = state->X[state->s & MASK1] = num_MultModD (param->b, 
    state->X[(state->s - 1) & MASK1] + state->X[(state->s - state->K) &
     MASK1], 0.0, param->m);
   return  Z / param->m;
}
#endif

/*-----------------------------------------------------------------------*/

static unsigned long DX02a_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * DX02a_U01 (vpar, vsta));
}


/*-----------------------------------------------------------------------*/

unif01_Gen *udeng_CreateDX02a (unsigned long m, unsigned long b, int k,
   unsigned long S[])
{
   return CreateDenga (m, b, k, S, 2);
}


/*=========================================================================*/

static double DL00a_U01 (void *vpar, void *vsta)
{
   DX02_state *state = vsta;
   DX02_param *param = vpar;
   double Z;
   ++state->s;

   Z = param->b * state->X[(state->s - state->K) & MASK1] -
          state->X[(state->s - 1) & MASK1];
   Z = state->X[state->s & MASK1] = fmod (Z + param->m, param->m);

   /*   Z = state->X[state->s & MASK1] = num_MultModD (param->b, 
       state->X[(state->s - state->K) & MASK1],
       -state->X[(state->s - 1) & MASK1], param->m); */
   return  Z / param->m;
}


/*-----------------------------------------------------------------------*/

static unsigned long DL00a_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * DL00a_U01 (vpar, vsta));
}


/*-----------------------------------------------------------------------*/

unif01_Gen *udeng_CreateDL00a (unsigned long m, unsigned long b, int k,
   unsigned long S[])
{
   return CreateDenga (m, b, k, S, 0);
}


/*=========================================================================*/

static void WrDX02a (void *vsta)
{
   DX02_state *state = vsta;
   unsigned int j;
   int s = state->s & MASK1;

   if (unif01_WrLongStateFlag || (state->K < 8)) {
      printf (" S = {\n ");
      for (j = 0; j < state->K; j++) {
         printf (" %12.0f", state->X[s]);
         if (--s < 0)
            s = MASK1;
         if (j < state->K - 1)
            printf (",");
         if ((j % 5) == 4)
            printf ("\n ");
      };
      printf ("   }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

static unif01_Gen *CreateDenga (unsigned long m, unsigned long b, int k,
   unsigned long S[], int id)
{
   unif01_Gen *gen;
   DX02_state *state;
   DX02_param *param;
   size_t leng;
   char name[LEN + 1];
   int j;

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (DX02_state));
   param = util_Malloc (sizeof (DX02_param));
   if (0 == id) {
      util_Assert (k < 129, "udeng_CreateDL00a:   k > 128");
   } else {
      util_Assert (k < 129, "udeng_CreateDX02a:   k > 128");
   }

   state->X = util_Calloc (MASK1 + 1, sizeof (double));
   for (j = 0; j < k; j++)
      state->X[k - j - 1] = S[j] % m;
   state->s = k - 1;
   state->K = k;
   param->b = b;
   param->m = m;

   if (0 == id)
      strcpy (name, "udeng_CreateDL00a:");
   else
#ifdef LAC
      strcpy (name, "udeng_CreateDX02a, Lac = {0, 101, 102}:");
#else
      strcpy (name, "udeng_CreateDX02a:");
#endif
   addstr_Ulong (name, "   m = ", m);
   addstr_Ulong (name, ",   b = ", b);
   addstr_Uint (name, ",   k = ", (unsigned) k);
   addstr_ArrayUlong (name, ",   S = ", k, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->param   = param;
   gen->state   = state;
   if (0 == id) {
      gen->GetBits = &DL00a_Bits;
      gen->GetU01  = &DL00a_U01;
   } else {
      gen->GetBits = &DX02a_Bits;
      gen->GetU01  = &DX02a_U01;
   }
   gen->Write   = &WrDX02a;
   return gen;
}


/*=========================================================================*/

void udeng_DeleteGen (unif01_Gen * gen)
{
   DX02_state *state;
   if (NULL == gen)
      return;
   state = gen->state;
   util_Free (state->X);
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/*=========================================================================*/
