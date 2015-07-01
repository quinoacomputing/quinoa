
#include<stdio.h>
#include<mt19937.h>

#define NN       100000000UL

int main(void){ 
   long i; unsigned int sum=0;
   mt19937_state state; mt19937_sse_state sse_state;
   mt19937_init_(&state); mt19937_get_sse_state_(&state,&sse_state);
   mt19937_print_state_(&state);
   for(i=0;i<NN;i++) sum+=mt19937_sse_generate_(&sse_state);
   printf("%ld MT19937_SSE pseudorandom numbers generated.\n",NN);
   printf("Fractional part of the total sum of generated numbers: %f\n",sum/4294967296.);
   printf("Next output value: %f\n",mt19937_sse_generate_(&sse_state)/4294967296.);
   return 0;
}
