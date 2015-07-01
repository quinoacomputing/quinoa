
#include<stdio.h>
#include<lfsr113.h>

#define NN       100000000UL

int main(void){ 
   long i; unsigned int sum=0;
   lfsr113_state state;
   lfsr113_init_(&state);
   lfsr113_print_state_(&state);
   for(i=0;i<NN;i++) sum+=lfsr113_generate_(&state);
   printf("%ld LFSR113 pseudorandom numbers generated.\n",NN);
   printf("Fractional part of the total sum of generated numbers: %f\n",sum/4294967296.);
   printf("Next output value: %f\n",lfsr113_generate_(&state)/4294967296.);
   return 0;
}
