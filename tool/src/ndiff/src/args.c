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

#include <stdlib.h>
#include <string.h>

#include "args.h"
#include "utils.h"
#include "utest.h"
#include "error.h"
#include "ndiff.h"
#include "context.h"
#include "register.h"

#ifndef VERSION
#define VERSION "2013.04.15"
#endif

#ifndef MAXREGS
#define MAXREGS 99

#elif   MAXREGS < 20
#undef  MAXREGS
#define MAXREGS 20

#elif   MAXREGS > REG_MAX
#undef  MAXREGS
#define MAXREGS REG_MAX
#endif

#ifndef MAXKEEP
#define MAXKEEP 25
#endif

#ifndef CMTCHRS
#define CMTCHRS ""
#endif

#ifndef PUNCTCHRS
#define PUNCTCHRS "._$"
#endif

#ifndef REG0FMT
#define REG0FMT "%g "
#endif

#ifndef SUITEFMT
#define SUITEFMT "[ %s ]"
#endif

#ifndef SERIEFMT
#define SERIEFMT "%d"
#endif

#ifndef OUTFILEEXT
#define OUTFILEEXT ".out"
#endif

#ifndef REFFILEEXT
#define REFFILEEXT ".ref"
#endif

#ifndef CFGFILEEXT
#define CFGFILEEXT ".cfg"
#endif

#ifndef RESFILEEXT
#define RESFILEEXT ".res"
#endif

#ifndef UNZIPCMD
#define UNZIPCMD "unzip -cq"
#endif

#ifndef GZIPCMD
#define GZIPCMD "gzip -cdq"
#endif

#ifndef BZIP2CMD
#define BZIP2CMD "bzip2 -cdq"
#endif

struct option option = {
  // index of processed option
  .argi = 1,

  // series numbering
  .fmt = SERIEFMT,

  // suite title
  .sfmt = SUITEFMT,

  // register print format
  .rfmt = REG0FMT,

  // comment characters
  .cchr = CMTCHRS,

  // punctuation characters part of identifiers
  .pchr = PUNCTCHRS,

  // number of diff displayed by default
  .keep = MAXKEEP,

  // number of registers allocated by default
  .nregs = MAXREGS,

  // file extensions
  .out_e = OUTFILEEXT, .ref_e = REFFILEEXT,
  .cfg_e = CFGFILEEXT, .res_e = RESFILEEXT,

  // unzip commands
  .unzip = { UNZIPCMD, GZIPCMD, BZIP2CMD}
};

static void
run_utest(void)
{
  struct utest *ut = utest_alloc(0);

  inform("Running unit tests (incomplete)");

  // list of unit tests: TODO: more utests
  context_utest(ut);
  ndiff_utest(ut);

  // stat
  utest_stat(ut);
  utest_free(ut);
}

void
invalid_option(const char *str)
{
  warning("invalid program options or arguments '%s'", str);
  usage();
}

void
invalid_file(const char *str)
{
  warning("invalid filename argument '%s'", str);
  usage();
}

void
usage(void)
{
  logmsg_config.level = inform_level;

  inform("usage:");
  inform("\tndiff [options] fileA[.out] fileB[.ref] [fileC[.cfg]]");
  inform("\tndiff [options] --list fileA fileB ...");
  inform("\tndiff [options] --test '1st' fileA fileB --test '2nd' fileC ...");

  inform("");
  inform("options:");
  inform("\t-a  --accum file    accumulate tests information in file");
  inform("\t-b  --blank         ignore blank spaces (space and tabs)");
  inform("\t    --cfgext ext    specify the config file extension, default is \"%s\"", option.cfg_e);
  inform("\t-c  --comment chrs  comment characters, default is \"%s\"", option.cchr);
  inform("\t-d  --debug         enable debug mode (include xcheck mode)");
  inform("\t-h  --help          display this help");
  inform("\t-i  --info          enable info mode (default)");
  inform("\t-k  --keep num      specify the number of diffs to display per file, default is %d", option.keep);
  inform("\t    --lhsrec        recycle next left file (exclusive with --rhsrec)");
  inform("\t    --lhsres        echo valid lines of next left file to its result file");
  inform("\t-l  --list          enable list mode (list of filenames)");
  inform("\t    --long          disable short options");
  inform("\t    --nocolor       disable color output for PASS/FAIL");
  inform("\t    --noloc         disable C file location during trace");
  inform("\t    --nowarn        disable warnings");
  inform("\t    --nregs num     specify the number of registers to allocate");
  inform("\t    --outext ext    specify the output file extension, default is \"%s\"", option.out_e);
  inform("\t    --punct chrs    punctuation characters part of identifiers, default is \"%s\"", option.pchr);
  inform("\t-q  --quiet         enable quiet mode (no output if no diff)");
  inform("\t    --refext ext    specify the reference file extension, default is \"%s\"", option.ref_e);
  inform("\t    --regfmt fmt    specify the (printf) format fmt for register 0, default is \"%s\"", option.rfmt);
  inform("\t-r  --reset         reset accumulated information");
  inform("\t    --resext ext    specify the result file extension, default is \"%s\"", option.res_e);
  inform("\t    --rhsrec        recycle next right file (exclusive with --lhsrec)");
  inform("\t    --rhsres        echo valid lines of next right file to its result file");
  inform("\t-n  --serie         enable serie mode (indexed filenames)");
  inform("\t    --seriefmt fmt  specify the (printf) format fmt for indexes, default is \"%s\"", option.fmt);
  inform("\t-s  --suite name    set test suite name for output message (title)");
  inform("\t    --suitefmt fmt  specify the (printf) format fmt for testsuite, default is \"%s\"", option.sfmt);
  inform("\t-t  --test name     set test name for output message (item)");
  inform("\t    --trace         enable trace mode (very verbose, include debug mode)");
  inform("\t    --trunc         allow premature ending of one of the input file");
  inform("\t    --utest         run the ndiff unit tests (still incomplete)");
  inform("\t-x  --xcheck        enable cross check mode (algorithms cross check)");

  inform("");
  inform("decompression:");
  inform("\t    --bzip2 cmd     command to uncompress .bz .bz2 .tbz .tbz2 files, default is \"%s\"", option.unzip[2]);
  inform("\t    --gzip  cmd     command to uncompress .gz .z .Z .tgz .taz .taZ files, default is \"%s\"", option.unzip[1]);
  inform("\t    --unzip cmd     command to uncompress .zip files, default is \"%s\"", option.unzip[0]);

  inform("");
  inform("rules (%s):", option.cfg_e);
  inform("\t#rows       cols          commands");
  inform("\t 1-5        *             skip                            # banner");
  inform("\t *          2-$           any abs=1e-15 rel=1e-12 dig=1.5 # global");
  inform("\t 41         *             goto='penalty function'         # jump");
  inform("\t 109:20/5   2-8/3         abs=1e-8                        # specific");

  inform("");
  inform("ranges:");
  inform("\tnum                 row or column number, num >= 0");
  inform("\trange               start - end  [/ stride]");
  inform("\tslice               start : size [/ stride]");
  inform("\t$, *                last row or column, alias for 0-$");

  inform("");
  inform("commands:");
  inform("\tabs=num or reg      absolute error (0 <= num <= 1)");
  inform("\t-abs=num or reg     negative absolute error (-1 <= num <= 0)");
  inform("\tall                 constraints are conjunctive (default, qualifier)");
  inform("\talt                 declare the rule as an alternate rule (qualifier)");
  inform("\tany                 constraints are disjunctive (qualifier)");
  inform("\tdig=num or reg      input-defined relative error (num >= 1)");
  inform("\t-dig=num or reg     input-defined negative relative error (num <= -1)");
  inform("\tequ                 strict numbers equality (same representation)");
  inform("\teval                perform operations even if rule fails");
  inform("\tgoto='tag'          skip lines until string 'tag' is found (action)");
  inform("\tgoto='num' or reg   skip lines until number 'num' is found (action)");
  inform("\tign                 ignore numbers, accept missing number if with istr");
  inform("\tistr                ignore strings while scanning for numbers");
  inform("\tlarge               allow num > 1 in  abs and  rel ");
  inform("\t                    and  num < -1 in -abs and -rel (qualifier)");
  inform("\tlhs=num or reg      set left hand side 'x'");
  inform("\tnofail              do not count nor display warning for failure");
  inform("\toff=num or reg      set error offset 'b'");
  inform("\tomit='tag'          ignore strings or numbers if preceded by 'tag'");
  inform("\trel=num or reg      relative error (0 <= num <= 1)");
  inform("\t-rel=num or reg     negative relative error (-1 <= num <= 0)");
  inform("\trhs=num or reg      set right hand side 'y'");
  inform("\tscl=num or reg      set error scaling factor 'a'");
  inform("\tskip                skip lines (action)");
  inform("\tsmall               forbid num > 1 in  abs and  rel and");
  inform("\t                    num < -1 in -abs and -rel (default, qualifier)");
  inform("\ttrace               trace rule when active (debug, qualifier)");
  inform("\ttraceR              trace rule and modified registers when active");

  inform("");
  inform("registers:");
  inform("\tR1..R9              contain lhs, rhs, dif, err, abs, rel, dig, min, prec");
  inform("\t=Rn                 load value from register n");
  inform("\t=-Rn                load negated value from register n");
  inform("\t=/Rn                load inverted value from register n");
  inform("\t=\\Rn                load negated and inverted value from register n");
  inform("\t=^Rn                load the exponential value from register n");
  inform("\t=|Rn                load the absolute value from register n");
  inform("\t=[Rn                load the value from register n rounded toward zero");
  inform("\t=]Rn                load the value from register n rounded toward infty");
  inform("\tR0=                 print the value(s) on the console");
  inform("\tRn=Rp+Rq            load the sum of registers p and q to register n");
  inform("\tRn=Rp-Rq            load the difference of registers p and q to register n");
  inform("\tRn=Rp*Rq            load the product of registers p and q to register n");
  inform("\tRn=Rp/Rq            load the ratio of registers p and q to register n");
  inform("\tRn=Rp%%Rq            load the reminder of registers p and q to register n");
  inform("\tRn=Rp^Rq            load the power of registers p and q to register n");
  inform("\tRn=Rp<Rq            load the min of registers p and q to register n");
  inform("\tRn=Rp>Rq            load the max of registers p and q to register n");
  inform("\tRn=Rp~Rq            move registers p..q to registers n..n+q-p");

  inform("");
  inform("info   :\thttp://cern.ch/mad/ndiff");
  inform("author :\tlaurent.deniau@cern.ch");
  inform("version:\t%s", VERSION);
  inform("license:\tGPLv3");

  exit(EXIT_FAILURE);
}

void
clear_args(void)
{
  if (option.lhs_res || option.rhs_res) {
    debug("results file(s) cleared");
    option.lhs_res = 0;
    option.rhs_res = 0;
  }

  if (option.recycle) {
    debug("recycling file(s) cleared");
    option.recycle = ndiff_norecycle;
  }
}

void
parse_args(int argc, const char *argv[])
{
  // parse command line arguments
  for (; option.argi < argc; option.argi++) {

    // not an option
    if (!is_option(argv[option.argi])) return;

    // set accumulation filename [setup]
    if (!strcmp(argv[option.argi], "--accum") || (!option.lgopt && !strcmp(argv[option.argi], "-a"))) {
      option.accum = argv[++option.argi]; 
      debug("accumulation filename set to '%s'", option.accum);
      continue;
    }

    // set blank mode [setup]
    if (!strcmp(argv[option.argi], "--blank") || (!option.lgopt && !strcmp(argv[option.argi], "-b"))) {
      debug("blank spaces ignored");
      option.blank = 1;
      continue;
    }

    // set config extension [setup]
    if (!strcmp(argv[option.argi], "--cfgext")) {
      option.cfg_e = argv[++option.argi]; 
      debug("config extension set to '%s'", option.cfg_e);
      continue;
    }

    // set comment characters [setup]
    if (!strcmp(argv[option.argi], "--comment") || (!option.lgopt && !strcmp(argv[option.argi], "-c"))) {
      option.cchr = argv[++option.argi]; 
      debug("comment characters set to '%s'", option.cchr);
      continue;
    }

    // set debug mode [setup]
    if (!strcmp(argv[option.argi], "--debug") || (!option.lgopt && !strcmp(argv[option.argi], "-d"))) {
      logmsg_config.level = debug_level;
      logmsg_config.locate = 1;
      debug("debug mode on");
      option.debug = 1;
      option.check = 1;
      continue;
    }

    // display help [action]
    if (!strcmp(argv[option.argi], "--help") || (!option.lgopt && !strcmp(argv[option.argi], "-h"))) {
      usage();
      continue;
    }

    // set info mode [setup]
    if (!strcmp(argv[option.argi], "--info") || (!option.lgopt && !strcmp(argv[option.argi], "-i"))) {
      debug("info mode on");
      logmsg_config.level = inform_level;
      logmsg_config.locate = 0;
      continue;
    }

    // set keep number [setup]
    if (!strcmp(argv[option.argi], "--keep") || (!option.lgopt && !strcmp(argv[option.argi], "-k"))) {
      option.keep = strtoul(argv[++option.argi],0,0);
      debug("keep set to %d", option.keep);
      continue;
    }

    // enable left result [setup]
    if (!strcmp(argv[option.argi], "--lhsrec")) {
      debug("recycling left file enabled");
      option.recycle = ndiff_recycle_left;
      continue;
    }

    // enable left result [setup]
    if (!strcmp(argv[option.argi], "--lhsres")) {
      debug("left results enabled");
      option.lhs_res = 1;
      continue;
    }

    // set list mode [setup]
    if (!strcmp(argv[option.argi], "--list") || (!option.lgopt && !strcmp(argv[option.argi], "-l"))) {
      debug("list mode on");
      option.list = 1;
      continue;
    }

    // disable short options [setup]
    if (!strcmp(argv[option.argi], "--long")) {
      debug("short options disabled");
      option.lgopt = 1;
      continue;
    }


    // disable color [setup]
    if (!strcmp(argv[option.argi], "--nocolor")) {
      debug("color output disabled");
      fail_str = "FAIL";
      pass_str = "PASS";
      continue;
    }

    // disable location trace [setup]
    if (!strcmp(argv[option.argi], "--noloc")) {
      debug("trace of location disabled");
      logmsg_config.locate = 0;
      continue;
    }

    // disable warnings [setup]
    if (!strcmp(argv[option.argi], "--nowarn")) {
      debug("no warning mode on");
      logmsg_config.level = error_level;
      logmsg_config.locate = 0;
      option.nowarn = 1;
      continue;
    }

    // set number of registers [setup]
    if (!strcmp(argv[option.argi], "--nregs")) {
      option.nregs = imin(REG_MAX, strtoul(argv[++option.argi],0,0));
      debug("number of registers set to %d", option.nregs);
      continue;
    }

    // set output extension [setup]
    if (!strcmp(argv[option.argi], "--outext")) {
      option.out_e = argv[++option.argi]; 
      debug("output extension set to '%s'", option.out_e);
      continue;
    }

    // set punctuation characters [setup]
    if (!strcmp(argv[option.argi], "--punct")) {
      option.pchr = argv[++option.argi]; 
      debug("punctuation characters set to '%s'", option.pchr);
      continue;
    }

    // set quiet mode [setup]
    if (!strcmp(argv[option.argi], "--quiet") || (!option.lgopt && !strcmp(argv[option.argi], "-q"))) {
      debug("quiet mode on");
      logmsg_config.level = warning_level;
      logmsg_config.locate = 0;
      continue;
    }

    // set reference extension [setup]
    if (!strcmp(argv[option.argi], "--refext")) {
      option.ref_e = argv[++option.argi]; 
      debug("reference extension set to '%s'", option.ref_e);
      continue;
    }

    // set register format [setup]
    if (!strcmp(argv[option.argi], "--regfmt")) {
      option.rfmt = argv[++option.argi];
      debug("register format set to '%s'", option.rfmt);
      continue;
    }

    // set result extension [setup]
    if (!strcmp(argv[option.argi], "--resext")) {
      option.res_e = argv[++option.argi]; 
      debug("result extension set to '%s'", option.res_e);
      continue;
    }

    // reset accumulation information [action]
    if (!strcmp(argv[option.argi], "--reset") || (!option.lgopt && !strcmp(argv[option.argi], "-r"))) {
      ensure(option.accum, "no accumulation file specified");
      debug("reseting file '%s'", option.accum);
      option.reset = 1;
      accum_summary(0, 0, 0, 0);
      continue;
    }

    // enable right result [setup]
    if (!strcmp(argv[option.argi], "--rhsrec")) {
      debug("recycling right file enabled");
      option.recycle = ndiff_recycle_right;
      continue;
    }

    // enable right result [setup]
    if (!strcmp(argv[option.argi], "--rhsres")) {
      debug("right results enabled");
      option.rhs_res = 1;
      continue;
    }

    // set serie mode [setup]
    if (!strcmp(argv[option.argi], "--serie") || (!option.lgopt && !strcmp(argv[option.argi], "-n"))) {
      debug("serie mode on");
      option.serie = 1;
      continue;
    }

    // set serie format [setup]
    if (!strcmp(argv[option.argi], "--seriefmt")) {
      option.fmt = argv[++option.argi];
      debug("serie format set to '%s'", option.fmt);
      continue;
    }

    // set suite name [setup]
    if (!strcmp(argv[option.argi], "--suite") || (!option.lgopt && !strcmp(argv[option.argi], "-s"))) {
      option.suite = argv[++option.argi];
      debug("suite name set to '%s'", option.suite);
      continue;
    }

    // set suite format [setup]
    if (!strcmp(argv[option.argi], "--suitefmt")) {
      option.sfmt = argv[++option.argi];
      debug("suite format set to '%s'", option.sfmt);
      continue;
    }

    // set test name [setup]
    if (!strcmp(argv[option.argi], "--test") || (!option.lgopt && !strcmp(argv[option.argi], "-t"))) {
      option.test = argv[++option.argi];
      debug("test name set to '%s'", option.test);
      continue;
    }

    // set trace mode [setup]
    if (!strcmp(argv[option.argi], "--trace")) {
      logmsg_config.level = trace_level;
      logmsg_config.locate = 1;
      debug("trace mode on");
      option.debug = 1;
      option.check = 1;
      continue;
    }

    // enable truncation [setup]
    if (!strcmp(argv[option.argi], "--trunc")) {
      debug("premature truncation allowed");
      option.trunc = 1;
      continue;
    }

    // run utests [action]
    if (!strcmp(argv[option.argi], "--utest")) {
      run_utest();
      option.utest += 1;
      continue;
    }

    // set check mode [setup]
    if (!strcmp(argv[option.argi], "--xcheck") || (!option.lgopt && !strcmp(argv[option.argi], "-x"))) {
      debug("check mode on");
      option.check = 1;
      continue;
    }

// ---- [decompression]

    // set tertiary unzip command [setup]
    if (!strcmp(argv[option.argi], "--bzip2")) {
      option.unzip[2] = argv[++option.argi]; 
      debug("bzip2 command set to '%s'", option.unzip[2]);
      continue;
    }

    // set secondary unzip command [setup]
    if (!strcmp(argv[option.argi], "--gzip")) {
      option.unzip[1] = argv[++option.argi]; 
      debug("gzip command set to '%s'", option.unzip[1]);
      continue;
    }

    // set primary unzip command [setup]
    if (!strcmp(argv[option.argi], "--unzip")) {
      option.unzip[0] = argv[++option.argi]; 
      debug("unzip command set to '%s'", option.unzip[0]);
      continue;
    }

// ---- [unknown]
    invalid_option(argv[option.argi]);
  }

  exit(EXIT_SUCCESS);
}

