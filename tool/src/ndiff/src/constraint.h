#ifndef CONSTRAINT_H
#define CONSTRAINT_H

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
     create constraints content
     print, scan constraints from file
 
 o---------------------------------------------------------------------o
*/

#include <stdio.h>
#include <string.h>
#include "slice.h"
#include "error.h"

// ----- constants

#ifndef MAXREGOP
#define MAXREGOP 10

#elif   MAXREGOP < 5
#undef  MAXREGOP
#define MAXREGOP 5

#elif   MAXREGOP > 32
#undef  MAXREGOP
#define MAXREGOP 32
#endif

#ifndef MAXTAGLEN
#define MAXTAGLEN 64

#elif   MAXTAGLEN < 48
#undef  MAXTAGLEN
#define MAXTAGLEN 48

#elif   MAXTAGLEN > 256
#undef  MAXTAGLEN
#define MAXTAGLEN 256
#endif

enum eps_cmd {
  eps_invalid = 0u,       // invalid command

// must be firsts (constraints)
  eps_abs    = 1u <<  0,  // absolute eps
  eps_rel    = 1u <<  1,  // relative eps
  eps_dig    = 1u <<  2,  // relative input eps
  eps_equ    = 1u <<  3,  // equal string
  eps_ign    = 1u <<  4,  // ignore value
  eps_istr   = 1u <<  5,  // ignore strings

// intermediate (commands & qualifiers)
  eps_lhs    = 1u <<  6,  // load x command
  eps_rhs    = 1u <<  7,  // load y command

  eps_any    = 1u <<  8,  // any qualifier
  eps_alt    = 1u <<  9,  // alt qualifier
  eps_eval   = 1u << 10,  // eval qualifier
  eps_omit   = 1u << 11,  // omit qualifier
  eps_nofail = 1u << 12,  // no failure qualifier
  eps_trace  = 1u << 13,  // trace qualifier
  eps_traceR = 1u << 14,  // trace Rn qualifier

  eps_swap   = 1u << 15,  // swap lhs and rhs (reserved)
  eps_onfail = 1u << 16,  // check alt on failure (reserved)

// must be lasts (actions)
  eps_skip   = 1u << 17,  // skip line, must be the first action!!
  eps_goto   = 1u << 18,  // go to tag
  eps_gonum  = 1u << 19,  // go to number

// marker & mask
  eps_last   = 1u << 20,  // the end
  eps_mask   = eps_last - 1,

// temporary commands & qualifiers
  eps_large  = eps_last << 0,  // large tolerance

// unions
  eps_dra    = eps_abs  | eps_rel  | eps_dig,
  eps_chk    = eps_dra  | eps_equ  | eps_ign,
  eps_sgg    = eps_skip | eps_goto | eps_gonum,
};

// ----- types

struct eps {
  enum eps_cmd cmd;

  // values
  double  lhs,  rhs;
  double  scl,  off;
  double  abs,  rel,  dig;
  double _abs, _rel, _dig;
  double  num;

  // loads
  short   lhs_reg,  rhs_reg;
  short   scl_reg,  off_reg;
  short   abs_reg,  rel_reg,  dig_reg;
  short  _abs_reg, _rel_reg, _dig_reg;
  short   gto_reg;

  // operations
  short   dst[MAXREGOP], src[MAXREGOP], src2[MAXREGOP], op_n;
  char    op [MAXREGOP];

  // tags
  char    tag[MAXTAGLEN];
}; // ~ 250 bytes

struct constraint {
  struct slice row;
  struct slice col;
  struct eps   eps;
  int    idx, line;
}; // ~ 280 bytes

// ----- interface

#define T struct constraint

static inline struct eps
eps_initAllNum(enum eps_cmd cmd, double abs, double rel, double dig, double _abs, double _rel, double _dig, double scl, double off)
{
  ensure(cmd > eps_invalid && cmd < eps_last, "invalid eps command");
  return (struct eps) { .cmd=cmd, .abs= abs,  .rel= rel,  .dig= dig,
                                 ._abs=_abs, ._rel=_rel, ._dig=_dig, .scl=scl, .off=off };
}

static inline struct eps
eps_initNum(enum eps_cmd cmd, double abs, double rel, double dig, double scl, double off)
{
  return eps_initAllNum(cmd, abs, rel, dig, -abs, -rel, -dig, scl, off);
}

static inline struct eps
eps_init(enum eps_cmd cmd, double val)
{
  return eps_initNum(cmd, cmd&eps_abs?val:0, cmd&eps_rel?val:0, cmd&eps_dig?val:0, 1.0, 0.0);
}

static inline struct eps
eps_initStrTag(enum eps_cmd cmd, const char *tag)
{
  ensure((cmd & eps_goto) || (cmd & eps_omit), "invalid eps goto or omit command");
  struct eps eps = (struct eps) { .cmd=cmd };
  enum { sz = sizeof eps.tag };
  strncpy(eps.tag, tag, sz); eps.tag[sz-1] = 0;
  return eps;
}

static inline struct eps
eps_initNumTag(enum eps_cmd cmd, const char *tag)
{
  char *end;
  struct eps eps = (struct eps) { .cmd=cmd };
  eps.num = strtod(tag, &end);
  ensure((cmd & eps_gonum) && !*end, "invalid eps goto command");
  enum { sz = sizeof eps.tag };
  strncpy(eps.tag, tag, sz); eps.tag[sz-1] = 0;
  return eps;
}

static inline T
constraint_init(const struct slice row, const struct slice col, const struct eps eps, int idx, int line)
{
  bool allcol = eps.cmd & eps_skip || eps.cmd & eps_goto;
  return (T){ row, allcol ? slice_initAll() : col, eps, idx, line };
}

void constraint_print(const T* cst, FILE *out);
void constraint_scan (      T* cst, FILE *in, int *row);

#undef T

#endif

