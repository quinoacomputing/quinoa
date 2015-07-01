
#include<stdio.h>
#include<gm31.h>

#define NN       100000000UL

int main(void){ 
   long i; unsigned int sum=0;
   gm31_state state;
   gm31_init_(&state);
   gm31_print_state_(&state);
   for(i=0;i<NN;i++) sum+=gm31_generate_(&state);
   printf("%ld GM31 pseudorandom numbers generated.\n",NN);
   printf("Fractional part of the total sum of generated numbers: %f\n",sum/4294967296.);
   printf("Next output value: %f\n",gm31_generate_(&state)/4294967296.);
   return 0;
}
