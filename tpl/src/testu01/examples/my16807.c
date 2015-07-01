#include "my16807.h"
#include "unif01.h"
#include "util.h"
#include "addstr.h"
#include <string.h>

typedef struct { double S; } My16807_state;

static double My16807_U01 (void *par, void *sta)
{
   My16807_state *state = sta;
   long k;
   state->S *= 16807.0;
   k = state->S / 2147483647.0;
   state->S -= k * 2147483647.0;
   return (state->S * 4.656612875245797E-10);
}

static unsigned long My16807_Bits (void *par, void *sta)
{
   return (unsigned long) (My16807_U01 (par, sta) * 4294967296.0);
}

static void WrMy16807 (void *sta)
{
   My16807_state *state = sta;
   printf (" S = %.0f\n", state->S);
}

unif01_Gen *CreateMy16807 (int s)
{
   unif01_Gen *gen;
   My16807_state *state;
   size_t leng;
   char name[60];

   gen = util_Malloc (sizeof (unif01_Gen));
   gen->state = state = util_Malloc (sizeof (My16807_state));
   state->S = s;
   gen->param = NULL;
   gen->Write = WrMy16807;
   gen->GetU01 = My16807_U01;
   gen->GetBits = My16807_Bits;

   strcpy (name, "My LCG implementation for a = 16807:");
   addstr_Int (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);
   return gen;
}

void DeleteMy16807 (unif01_Gen * gen)
{
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}
