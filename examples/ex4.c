
#include "unif01.h"
#include "utaus.h"
#include <stdio.h>

int main (void)
{
   unif01_Gen *gen1, *gen2, *gen3;
   long I[3] = { 3, 7, 9 };
   int i, n = 20;
   double x;

   gen1 = utaus_CreateTaus (31, 3, 12, 12345);
   for (i = 0; i < n; i++)
      printf ("%f\n", unif01_StripD (gen1, 0));
   utaus_DeleteGen (gen1);
   printf ("\n");

   gen1 = utaus_CreateTaus (31, 3, 12, 12345);
   gen2 = unif01_CreateLacGen (gen1, 3, I);
   for (i = 0; i < n; i++)
      printf ("%f\n", unif01_StripD (gen2, 0));

   gen3 = unif01_CreateDoubleGen (gen2, 24);
   for (i = 0; i < n; i++)
      x = unif01_StripD (gen3, 0);
   unif01_DeleteDoubleGen (gen3);
   unif01_DeleteLacGen (gen2);

   gen2 = utaus_CreateTaus (28, 7, 14, 12345);
   gen3 = unif01_CreateCombXor2 (gen1, gen2, "A Combined Tausworthe Gener.");
   for (i = 0; i < n; i++)
      x = unif01_StripD (gen3, 0);
   unif01_DeleteCombGen (gen3);
   utaus_DeleteGen (gen2);
   utaus_DeleteGen (gen1);
   return 0;
}
