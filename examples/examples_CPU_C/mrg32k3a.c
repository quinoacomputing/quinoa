
#include<stdio.h>
#include<mrg32k3a.h>

#define NN       100000000UL

int main(void){ 
   long i; unsigned int sum=0;
   mrg32k3a_state state;
   mrg32k3a_init_(&state);
   mrg32k3a_print_state_(&state);
   for(i=0;i<NN;i++) sum+=mrg32k3a_generate_(&state);
   printf("%ld MRG32K3A pseudorandom numbers generated.\n",NN);
   printf("Fractional part of the total sum of generated numbers: %f\n",sum/4294967296.);
   printf("Next output value: %f\n",mrg32k3a_generate_(&state)/4294967296.);
   return 0;
}
