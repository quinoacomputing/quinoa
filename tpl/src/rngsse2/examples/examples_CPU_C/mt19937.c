
#include<stdio.h>
#include<mt19937.h>

#define NN       100000000UL

int main(void){ 
   long i; unsigned int sum=0;
   mt19937_state state;
   mt19937_init_(&state);
   mt19937_print_state_(&state);
   for(i=0;i<NN;i++) sum+=mt19937_generate_(&state);
   printf("%ld MT19937 pseudorandom numbers generated.\n",NN);
   printf("Fractional part of the total sum of generated numbers: %f\n",sum/4294967296.);
   printf("Next output value: %f\n",mt19937_generate_(&state)/4294967296.);
   return 0;
}
