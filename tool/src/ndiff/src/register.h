#ifndef REGISTER_H
#define REGISTER_H

/*
 o---------------------------------------------------------------------o
 |
 | Ndiff
 |
 | Copyright (c) 2012+ laurent.deniau@cern.ch
 | Gnu General Public License
 |
 o---------------------------------------------------------------------o
  
   Purpose:
     manage access to registers
     handle arithmetic operations (+, -, *, /) and range (~)
 
 o---------------------------------------------------------------------o
*/

#include <stdlib.h>
#include "types.h"
#include "error.h"

// ----- constants

#define REG_MAX (1 << 13)
#define REG_UNARY_OP  "-/\\^|[]"
#define REG_BINARY_OP "+-*/%^<>~"

// ----- attributes

#ifdef __GNUC__
static inline short reg_encode(short rn, char  op) __attribute__((always_inline));
static inline short reg_decode(short rn, char *op) __attribute__((always_inline));
#endif

// ----- interface

static inline bool
reg_isvalid(short rn)
{
  return rn > 0 && rn < REG_MAX;
}

static inline short
reg_encode(short rn, char op)
{
  ensure(reg_isvalid(rn), "invalid register number %d", rn);

  switch(op) {
  case  0  : return 0 * REG_MAX + rn;
  case '-' : return 1 * REG_MAX + rn;
  case '/' : return 2 * REG_MAX + rn;
  case '\\': return 3 * REG_MAX + rn;
  case '^' : return 4 * REG_MAX + rn;
  case '|' : return 5 * REG_MAX + rn;
  case '[' : return 6 * REG_MAX + rn;
  case ']' : return 7 * REG_MAX + rn;
  default: error("invalid register unary operation '%c'", op); exit(EXIT_FAILURE);
  }
}

static inline short
reg_decode(short rn, char *op)
{
  switch(rn / REG_MAX) {
  case 0: if (op) *op =  0  ; break;
  case 1: if (op) *op = '-' ; break;
  case 2: if (op) *op = '/' ; break;
  case 3: if (op) *op = '\\'; break;
  case 4: if (op) *op = '^' ; break;
  case 5: if (op) *op = '|' ; break;
  case 6: if (op) *op = '[' ; break;
  case 7: if (op) *op = ']' ; break;
  default: error("invalid register unary operation '%c'", rn/REG_MAX);
  }
  rn &= REG_MAX-1;
  ensure(reg_isvalid(rn), "invalid register number %d", rn);
  return rn;
}

static inline void
reg_setval(double *reg, short reg_n, short rn, double val)
{
  ensure(rn > 0 && rn <= reg_n, "invalid register number %d", rn);
  reg[rn-1] = val;
}

double reg_getval(const double *reg, short reg_n, short rn);
void   reg_eval  (double       *reg, short reg_n, short dst, short src, short src2, char op);

#endif

