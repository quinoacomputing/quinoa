
#include "unif01.h"
#include "bbattery.h"

unsigned int xorshift (void);
double MRG32k3a (void);

int main (void) 
{
   unif01_Gen *gen;

   gen = unif01_CreateExternGen01 ("MRG32k3a", MRG32k3a);
   bbattery_SmallCrush (gen);
   unif01_DeleteExternGen01 (gen);

   gen = unif01_CreateExternGenBits ("xorshift", xorshift);
   bbattery_SmallCrush (gen);
   unif01_DeleteExternGenBits (gen);

   return 0;
}
