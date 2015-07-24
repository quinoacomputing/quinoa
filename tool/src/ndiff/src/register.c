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

#include <math.h>
#include <stdio.h>

#include "args.h"
#include "register.h"

// ----- interface

double
reg_getval(const double *reg, short reg_n, short rn)
{
  short r = rn & (REG_MAX-1);
  ensure(r > 0 && r <= reg_n, "invalid register R%d", r);

  switch(rn / REG_MAX) {
  case 0: return reg[r-1];
  case 1: return -reg[r-1];
  case 2: return 1/reg[r-1];
  case 3: return -1/reg[r-1];
  case 4: return exp(reg[r-1]);
  case 5: return fabs(reg[r-1]);
  case 6: return floor(reg[r-1]);
  case 7: return ceil(reg[r-1]);
  default: error("invalid register unary operation '%c'", rn/REG_MAX); exit(EXIT_FAILURE);
  }
}

#ifdef __GNUC__
__attribute__((always_inline))
#endif
static inline void
reg_eval_assign(double *reg, short reg_n, short dst, short src, short src2, char op)
{
  ensure(dst > 0 && dst <= reg_n, "invalid register R%d", dst);

  if (!op) { // nop = getval
    reg[dst-1] = reg_getval(reg, reg_n, src);
  }
  else {
    ensure(src  > 0 && src  <= reg_n, "invalid register R%d", src);
    ensure(src2 > 0 && src2 <= reg_n, "invalid register R%d", src2);

    switch(op) {
    case '+': reg[dst-1] = reg[src-1] + reg[src2-1]; break;
    case '-': reg[dst-1] = reg[src-1] - reg[src2-1]; break;
    case '*': reg[dst-1] = reg[src-1] * reg[src2-1]; break;
    case '/': reg[dst-1] = reg[src-1] / reg[src2-1]; break;
    case '%': reg[dst-1] = fmod(reg[src-1], reg[src2-1]); break;
    case '^': reg[dst-1] = pow(reg[src-1], reg[src2-1]); break;
    case '<': reg[dst-1] = reg[src-1] < reg[src2-1] ? reg[src-1] : reg[src2-1]; break;
    case '>': reg[dst-1] = reg[src-1] > reg[src2-1] ? reg[src-1] : reg[src2-1]; break;
    case '~': {
      short end = dst+src2-src;
      ensure(end > 0 && end <= reg_n, "invalid range of registers R%d~R%d", dst, end);
      ensure(src < src2, "invalid range of registers R%d~R%d", src, src2);
      if (dst <= src) // take care of overlapping registers
        for (short i=0; i <= src2-src; i++)    
          reg[dst+i-1] = reg[src+i-1];
      else
        for (short i=src2-src; i >= 0; i--)    
          reg[dst+i-1] = reg[src+i-1];
    } break;
    default:
      error("invalid register operation R%d'%c'R%d", src, op, src2);
    }
  }
}

#ifdef __GNUC__
__attribute__((always_inline))
#endif
static inline void
reg_eval_print(double *reg, short reg_n, short src, short src2, char op)
{
  if (!op) { // nop = getval
    printf(option.rfmt, reg_getval(reg, reg_n, src));
  }
  else {
    ensure(src  > 0 && src  <= reg_n, "invalid register R%d", src);
    ensure(src2 > 0 && src2 <= reg_n, "invalid register R%d", src2);

    switch(op) {
    case '+': printf(option.rfmt, reg[src-1] + reg[src2-1]); break;
    case '-': printf(option.rfmt, reg[src-1] - reg[src2-1]); break;
    case '*': printf(option.rfmt, reg[src-1] * reg[src2-1]); break;
    case '/': printf(option.rfmt, reg[src-1] / reg[src2-1]); break;
    case '%': printf(option.rfmt, fmod(reg[src-1], reg[src2-1])); break;
    case '^': printf(option.rfmt, pow(reg[src-1], reg[src2-1]));  break;
    case '<': printf(option.rfmt, reg[src-1] < reg[src2-1] ? reg[src-1] : reg[src2-1]); break;
    case '>': printf(option.rfmt, reg[src-1] > reg[src2-1] ? reg[src-1] : reg[src2-1]); break;
    case '~':
      ensure(src < src2, "invalid range of registers R%d~R%d", src, src2);
      for (short i=0; i <= src2-src; i++)    
        printf(option.rfmt, reg[src+i-1]);
      break;
    default:
      error("invalid register operation R%d'%c'R%d", src, op, src2);
    }
  }
  putchar('\n'); 
}

void
reg_eval(double *reg, short reg_n, short dst, short src, short src2, char op)
{
  ensure(dst >= 0 && dst <= reg_n, "invalid register %d", dst);

  if (dst)
    reg_eval_assign(reg, reg_n, dst, src, src2, op);
  else
    reg_eval_print (reg, reg_n,      src, src2, op);
}

