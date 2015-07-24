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
    main file and loop
 
 o---------------------------------------------------------------------o
*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "main.h"
#include "args.h"
#include "utils.h"
#include "error.h"
#include "ndiff.h"
#include "context.h"
#include "constraint.h"

static int
diff_summary(const struct ndiff *dif)
{
  int n, c;
  ndiff_getInfo(dif, &n, 0, &c, 0);

  if (!ndiff_feof(dif, 1) && !option.trunc) {
    c += 1;
    warning("diff of '%s'|'%s' ended prematurely (truncated output?)", option.lhs_file, option.rhs_file);
  }

  if (c) {
//    if (option.test)
//    warning("(*) files '%s'|'%s' from test '%s' differ", option.lhs_file, option.rhs_file, option.test);
    warning("(=) % 6d lines have been diffed", n);
    warning("(=) % 6d diffs have been detected", c);
  } else {
    if (option.test)
    inform ("files '%s'|'%s' from test '%s' do not differ", option.lhs_file, option.rhs_file, option.test);
    inform ("% 6d lines have been diffed", n);
  }

  return c;
}

static void
test_summary(int total, int failed)
{
  double t = (option.clk_t1 - option.clk_t0) / CLOCKS_PER_SEC;
  printf(" + %-50s (%.2f s) - %2d/%2d : %s\n", option.test, t, total-failed, total,
          failed ? fail_str : pass_str);
}

static void
check_transition(const char* argv[], int *total, int *failed, long lines, long numbers)
{
  if (is_option(argv[option.argi]) && option.test && *total && (
      (!strcmp(argv[option.argi], "-t") && !option.lgopt) || !strcmp(argv[option.argi], "--test" ) ||
      (!strcmp(argv[option.argi], "-s") && !option.lgopt) || !strcmp(argv[option.argi], "--suite"))) {

    // stop timer
    option.clk_t1 = clock();

    // test stats
    test_summary(*total, *failed);

    // suites stats
    if (option.accum)
      accum_summary(*total, *failed, lines, numbers);

    // cleanup
    *total = *failed = 0;

    // start timer
    option.clk_t0 = clock();
  }
}

// counters
static int  total, failed;
static long lines, numbers;

void
quit(int exit_code)
{
  if (option.test)
    test_summary(total, failed);

  if (option.accum)
    accum_summary(total, failed, lines, numbers);

  exit(exit_code);
}

int
main(int argc_, char** argv_)
{
  // get const copy
  const int    argc = argc_;
  const char **argv = (void*)argv_;

  // start timers
  option.dat_t0 = time(0);
  option.clk_t0 = clock();

  // set logging level
  logmsg_config.level = inform_level;

  // argument list loop (too long, should refactored)
  while (option.argi < argc) {
    const char *lhs_s = 0, *rhs_s = 0, *cfg_s = 0;
    int n = 0;

    // check for test or test suite transition (chaining)
    check_transition(argv, &total, &failed, lines, numbers);

    // parse arguments [incremental]
    parse_args(argc, argv);

    // setup filenames [incremental]
    if (!option.list) {
      int i;
      for (i = option.argi; i < option.argi+3; i++)
        if (i >= argc || is_option(argv[i])) break;
      if (option.argi < i) lhs_s = argv[option.argi++];
      if (option.argi < i) rhs_s = argv[option.argi++];
      if (option.argi < i) cfg_s = argv[option.argi++];
    } else
      if (option.argi < argc) lhs_s = rhs_s = cfg_s = argv[option.argi++];

    trace("arguments: total=%d, left=%d, right=%d, curr=%s",
          argc, option.argi, argc-option.argi, option.argi < argc ? argv[option.argi] : "nil");

    // no more files
    if (!lhs_s || !rhs_s) exit(EXIT_SUCCESS);

    // suite title (first time only)
    if (option.suite) {
      fprintf(stdout, option.sfmt, option.suite);
      fputc('\n', stdout);
      option.suite = 0;
    }

    // serie loop
    while (option.serie || !n) {
      FILE *lhs_fp=0, *rhs_fp=0, *cfg_fp=0, *lhs_rfp=0, *rhs_rfp=0;
      int nn = n;

      // clean filenames
      *option.lhs_file = *option.rhs_file = *option.cfg_file = 0;
       option.lhs_zip  =  option.rhs_zip  =  option.cfg_zip  = 0;

      // open files
      lhs_fp = open_file(lhs_s, option.lhs_res ? &lhs_rfp : 0, &nn, option.out_e, 1, 0);
      if (!lhs_fp && n) break; // end of serie
      rhs_fp = open_file(rhs_s, option.rhs_res ? &rhs_rfp : 0, &nn, option.ref_e, !option.list, 1);
      cfg_fp = open_file(cfg_s,                             0, &nn, option.cfg_e, !option.list, 0);
      if (n != nn) { n = nn; --total; }

      if (!lhs_fp) {
        if (option.list) {
          warning("output file '%s[.out]' not found, skipping diff", lhs_s);
          close_file(rhs_fp, option.rhs_zip);
          close_file(cfg_fp, option.cfg_zip);
          close_file(rhs_rfp, 0);
          ++failed;
          break;
        } else
          invalid_file(lhs_s);
      }

      if (!rhs_fp) {
        if (option.list) {
          warning("reference file '%s.ref' not found, skipping diff", rhs_s);
          close_file(lhs_fp, option.lhs_zip);
          close_file(cfg_fp, option.cfg_zip);
          close_file(lhs_rfp, 0);
          ++failed;
          break;
        } else
          invalid_file(rhs_s);
      }

      inform("processing '%s'|'%s'", option.lhs_file, option.rhs_file);

      // create context of constraints (using default size)
      struct context *cxt = context_alloc(0);

      // load constraints
      if (cfg_fp) cxt = context_scan(cxt, cfg_fp);

      // show constraints
      if (option.debug) {
        debug("rules list:");
        context_print(cxt, stderr);
      }

      // ndiff loop
      struct ndiff *dif = ndiff_alloc(lhs_fp, rhs_fp, cxt, 0, option.nregs);
      ndiff_option(dif, &option.keep, &option.blank, &option.check, &option.recycle);
      ndiff_result(dif, lhs_rfp, rhs_rfp);
      ndiff_loop(dif);

      // print summary
      if (diff_summary(dif) > 0) ++failed;

      // collect stats
      { int row; long num;
        ndiff_getInfo(dif, &row, 0, 0, &num);
        lines += row-1; numbers += num;
      }

      // destroy components
      ndiff_free(dif);
      context_free(cxt);

      // close files
      close_file(lhs_fp, option.lhs_zip);
      close_file(rhs_fp, option.rhs_zip);
      close_file(cfg_fp, option.cfg_zip);
      close_file(lhs_rfp, 0);
      close_file(rhs_rfp, 0);

      // clear options
      clear_args();

      n += 1;

      // not a serie, stop this loop
      if (!option.serie) break;
    }

    total += n;
  }

  option.clk_t1 = clock();

  quit(EXIT_SUCCESS);
}

