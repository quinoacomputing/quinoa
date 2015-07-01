#include "usoft.h"
#include "unif01.h"
#include "bbattery.h"

int main (void)
{
   unif01_Gen *gen;
   int n = 1024*1024;

   /* Test the first n bits of binary file vax.bin */
   bbattery_AlphabitFile ("vax.bin", n);

   /* Test the Java random number generator */
   gen = usoft_CreateJava48 (1234567, 1);
   bbattery_Alphabit (gen, n, 0, 32);

   return 0;
}
