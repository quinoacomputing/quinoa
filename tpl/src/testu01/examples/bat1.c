#include "ulcg.h"
#include "unif01.h"
#include "bbattery.h"

int main (void)
{
   unif01_Gen *gen;
   gen = ulcg_CreateLCG (2147483647, 16807, 0, 12345);
   bbattery_SmallCrush (gen);
   ulcg_DeleteGen (gen);
   return 0;
}
