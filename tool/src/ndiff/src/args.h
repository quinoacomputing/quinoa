#ifndef ARGS_H
#define ARGS_H

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
    manage arguments and options
    display the help
    run unit tests (immediate)
 
 o---------------------------------------------------------------------o
*/

#include <stdio.h>
#include <time.h>
#include "types.h"

struct option {
  int check, debug, nowarn, keep, lgopt;
  int serie, list, blank, utest, reset, trunc, nregs, recycle;
  const char *suite, *test;
  const char *fmt, *sfmt, *rfmt;
  const char *pchr, *cchr;
  const char *out_e, *ref_e, *cfg_e, *res_e;
  const char *unzip[3];
  char lhs_file[FILENAME_MAX];
  char rhs_file[FILENAME_MAX];
  char cfg_file[FILENAME_MAX];
  int  lhs_zip, rhs_zip, cfg_zip;
  int  lhs_res, rhs_res;
  int  argi;

  const char *accum;
  time_t dat_t0;
  double clk_t0, clk_t1;
};

extern struct option option;

void usage(void);
void invalid_file(const char*);
void invalid_option(const char*);
void parse_args(int argc, const char *argv[]);
void clear_args(void);

static inline bool
is_option(const char *arg)
{
  return arg[0] == '-' && (arg[1] == '-' || !arg[1] || !option.lgopt);
}

#endif
