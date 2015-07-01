
#include "unif01.h"
#include "ulcg.h"
#include "ulec.h"
#include "my16807.h"
#include <stdio.h>

int main (void) 
{
   unif01_Gen *gen;
   double x = 0.0;
   int i;

   gen = ulcg_CreateLCGFloat (2147483647, 16807, 0, 12345);
   unif01_TimerSumGenWr (gen, 10000000, TRUE);
   ulcg_DeleteGen (gen);

   gen = CreateMy16807 (12345);
   unif01_TimerSumGenWr (gen, 10000000, TRUE);
   DeleteMy16807 (gen);

   gen = ulec_CreateMRG32k3a (123., 123., 123., 123., 123., 123.);
   unif01_TimerSumGenWr (gen, 10000000, TRUE);
   ulec_DeleteGen (gen);

   gen = ulec_Createlfsr113 (12345, 12345, 12345, 12345);
   unif01_TimerSumGenWr (gen, 10000000, TRUE);
   for (i = 0; i < 100; i++)
      x += unif01_StripD (gen, 0);
   printf ("Sum = %14.10f\n", x);
   ulec_DeleteGen (gen);

   return 0;
}
