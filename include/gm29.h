// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#define gm29_g 536870909U

typedef struct{
  unsigned xN[32] __attribute__ ((aligned(16))), 
           xP[32] __attribute__ ((aligned(16)));
} gm29_state;

typedef gm29_state gm29_sse_state;

void gm29_skipahead_(gm29_state* state, unsigned long long offset);

void gm29_init_(gm29_state* state);

void gm29_init_short_sequence_(gm29_state* state,unsigned SequenceNumber);

void gm29_init_medium_sequence_(gm29_state* state,unsigned SequenceNumber);

void gm29_init_long_sequence_(gm29_state* state,unsigned SequenceNumber);

unsigned int gm29_generate_(gm29_state* state);

float gm29_generate_uniform_float_(gm29_state* state);

unsigned int gm29_sse_generate_(gm29_sse_state* state);

void gm29_get_sse_state_(gm29_state* state,gm29_sse_state* sse_state);

void gm29_print_state_(gm29_state* state);

void gm29_print_sse_state_(gm29_sse_state* state);
