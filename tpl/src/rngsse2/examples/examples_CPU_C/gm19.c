
#include<stdio.h>
#include<gm19.h>

#define NN       100000000UL

int main(void){ 
   long i; unsigned int sum=0;
   gm19_state state;
   gm19_init_(&state);
   gm19_print_state_(&state);
   for(i=0;i<NN;i++) sum+=gm19_generate_(&state);
   printf("%ld GM19 pseudorandom numbers generated.\n",NN);
   printf("Fractional part of the total sum of generated numbers: %f\n",sum/4294967296.);
   printf("Next output value: %f\n",gm19_generate_(&state)/4294967296.);
   return 0;
}
