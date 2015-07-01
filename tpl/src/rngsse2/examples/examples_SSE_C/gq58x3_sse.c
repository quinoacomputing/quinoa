
#include<stdio.h>
#include<gq58x3.h>

#define NN       100000000UL

int main(void){ 
   long i; unsigned int sum=0;
   gq58x3_state state; gq58x3_sse_state sse_state;
   gq58x3_init_(&state); gq58x3_get_sse_state_(&state,&sse_state);
   gq58x3_print_state_(&state);
   for(i=0;i<NN;i++) sum+=gq58x3_sse_generate_(&sse_state);
   printf("%ld GQ58X3_SSE pseudorandom numbers generated.\n",NN);
   printf("Fractional part of the total sum of generated numbers: %f\n",sum/4294967296.);
   printf("Next output value: %f\n",gq58x3_sse_generate_(&sse_state)/4294967296.);
   return 0;
}
