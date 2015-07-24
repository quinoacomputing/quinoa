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
     provides utilities
 
 o---------------------------------------------------------------------o
*/

#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "error.h"
#include "utils.h"
#include "args.h"

#ifdef _WIN32
#ifndef popen
#define popen(cmd, mode)  _popen(cmd, mode)
#endif
#ifndef pclose
#define pclose(fp)        _pclose(fp)
#endif
#endif

// macros

#if defined(_WIN32) || defined(NOCOLORS)
#define CSTR_RED(s) s
#define CSTR_GREEN(s) s
#else
#define CSTR_RED(s)   "\033[31m" s "\033[0m"
#define CSTR_GREEN(s) "\033[32m" s "\033[0m"
#endif

// constants

const char* fail_str = CSTR_RED  ("FAIL");
const char* pass_str = CSTR_GREEN("PASS");

static struct {
  const char *ext;
  int         cmd;
}
const zipext[] = {
  { ""        ,  0 },     // no compression
  { ".zip"    ,  1 },     // unzip      (unzip[0])
  { ".gz"     ,  2 },     // gzip       (unzip[1])
  { ".z"      ,  2 },     // gzip       (unzip[1])
  { ".Z"      ,  2 },     // gzip       (unzip[1])
  { ".tgz"    ,  2 },     // gzip       (unzip[1])
  { ".taz"    ,  2 },     // gzip       (unzip[1])
  { ".taZ"    ,  2 },     // gzip       (unzip[1])
  { ".bz"     ,  3 },     // bzip2      (unzip[2])
  { ".bz2"    ,  3 },     // bzip2      (unzip[2])
  { ".tbz"    ,  3 },     // bzip2      (unzip[2])
  { ".tbz2"   ,  3 },     // bzip2      (unzip[2])
  { ""        ,  0 },     // no compression, for clean error reporting
  {  NULL     ,  0 }      // end marker
};

// functions 

static bool
is_zipext(const char *str)
{
  for (int i = 0; zipext[i].ext; i++)
    if (!strcmp(zipext[i].ext, str))
      return i;

  return 0;
}

void
close_file(FILE *fp, int zip)
{
  if (fp && fp != stdin && fp != stdout) {
    if (zip) pclose(fp);
    else     fclose(fp);
  }
}

FILE*
open_file(const char* str, FILE **res_fp, int *idx, const char *ext, int optext, int required)
{
  char  buf[FILENAME_MAX+100];
  char rbuf[FILENAME_MAX+100];
  const char *dot = 0, *zdot = 0;
  int pos = 0, zpos = 0, zid = 0;
  FILE *fp = 0;

  assert(ext);

  // no file
  if (!str) return 0;

  // stdin
  if (str[0] == '-' && str[1] == 0) {
    fp = stdin;
    strncpy(buf, str, sizeof buf-1);
    goto filename;
  }

retry:

  for (int i=0; zipext[i].ext; i++) {

    // copy filename
    sprintf(buf, "%s%s", str, zipext[i].ext);

    // find and save zip extension, if any
    zdot = strrchr(buf, '.');
    if (zdot && (zid = is_zipext(zdot))) {
      zpos = zdot-buf;
      buf[zpos] = 0;
    } else {
      zpos = (int)strlen(buf);
      zdot = 0;
    }

    // find and save expected extension, if any
    dot = strrchr(buf, '.');
    if (dot && !strcmp(dot, ext)) {
      pos = dot-buf; 
      buf[pos] = 0;
    } else {
      pos = (int)strlen(buf);
      dot = 0;
    }

    // add formatted index, if in serie
    if (option.serie && idx && *idx > 0) {
      int n = sprintf(buf+pos, option.fmt, *idx);
      pos += n; zpos += n;
    }

    // add extensions (always for first attempt: procedure is safer)
    strncat(buf+pos , ext            , sizeof buf -  pos);
    strncat(buf+zpos, zipext[zid].ext, sizeof buf - zpos);

    // try to open, try again upon failure if extension is optional
    trace("trying to open file '%s' for reading", buf);
    fp = fopen(buf, "r");
    if (!fp && optext) {
      buf[pos] = 0;
      strncat(buf+pos, zipext[zid].ext, sizeof buf - pos);
      trace("trying to open file '%s' for reading", buf);
      fp = fopen(buf, "r");
    }

    if (fp || !option.list) break;
  }

  // allow failure on first non-numbered file for serie
  if (!fp) {
    if (option.serie && idx && *idx == 0) { ++*idx; goto retry; }
    ensure(!required, "unable to find file '%s.ref' or any %scompressed variants", str,
            idx && *idx ? "numbered/" : "");
  }

  // debug information
  if (!fp) {
    trace("<-open_file: unable to open file '%s' for reading", buf);
    return 0;
  }

  // close file if zipped, reopen through popen
  if (zid) {
    char zbuf[2*(FILENAME_MAX+100)];
    fclose(fp);
    sprintf(zbuf, "%s %s", option.unzip[zid-1], buf);
    debug("trying to reopen compressed file '%s' for reading", zbuf);
    fp = popen(zbuf, "r");
    ensure(fp, "failed to execute '%s'", zbuf);
  }

  // resize buffer for faster read
  if (BUFSIZ < 65536 && setvbuf(fp, 0, _IOFBF, 65536)) {
    close_file(fp, zid);
    error("unable to resize the stream buffer size");
  }

filename:

  // copy filenames into option for further reporting
  if (ext == option.out_e) {
    option.lhs_zip = zid;
    strncpy(option.lhs_file, buf, sizeof option.lhs_file-1);
    option.lhs_file[sizeof option.lhs_file-1] = 0;
  }

  if (ext == option.ref_e) {
    option.rhs_zip = zid;
    strncpy(option.rhs_file, buf, sizeof option.rhs_file-1);
    option.rhs_file[sizeof option.rhs_file-1] = 0;
    ensure(strcmp(option.lhs_file, option.rhs_file), "lhs and rhs files have same name");
  }

  if (ext == option.cfg_e) {
    option.cfg_zip = zid;
    strncpy(option.cfg_file, buf, sizeof option.cfg_file-1);
    option.cfg_file[sizeof option.cfg_file-1] = 0;
  }

  // open result file, if any
  if (res_fp) {
    strncpy(rbuf, buf         , sizeof rbuf-1); rbuf[sizeof rbuf-1]=0;
    strncat(rbuf, option.res_e, sizeof rbuf-1); rbuf[sizeof rbuf-1]=0;
    *res_fp = fopen(rbuf, "w");
  }

  // debug information
  debug("file %s open for reading", buf);
  if (res_fp && *res_fp) debug("file %s open for writing", rbuf);

  return fp;
}

void
accum_summary(int total, int failed, long lines, long numbers)
{
  if (!option.accum) return;

  int n;
  double tz;
  struct tm tm;
  time_t now = time(0);
  int total_tests=0, total_passed=0, total_failed=0;
  long total_lines=0, total_numbers=0;
  double total_ndtime=0;

  FILE *fp;
  if (!option.reset && (fp = fopen(option.accum, "r+"))) {
    // read time stamps
    n = fscanf(fp, " = tests summary (started at %d.%d.%d %d:%d:%d %d)\n",
                   &tm.tm_year, &tm.tm_mon, &tm.tm_mday,
                   &tm.tm_hour, &tm.tm_min, &tm.tm_sec, &tm.tm_isdst);
    ensure(n == 7, "invalid summary file format %s", option.accum);
    tm.tm_year -= 1900;
    tm.tm_mon  -= 1;

    // correct for TZ shift (i.e. emulate non-standard timegm)
    struct tm tm2 = tm;
    option.dat_t0 = mktime(&tm2);
    tm2 = *gmtime(&now);
    tz = difftime(now, mktime(&tm2));

    // read diff time, line and number counts
    n = fscanf(fp, "   total diff time %lf s  -  total lines %ld  -  total numbers %ld\n",
                   &total_ndtime, &total_lines, &total_numbers);
    ensure(n == 3, "invalid summary file format %s (lines)", option.accum);
    // read tests counts
    n = fscanf(fp, "   total run  time %*f s  -  total files %d - PASSED %d - FAILED %d\n",
                   &total_tests, &total_passed, &total_failed);
    ensure(n == 3, "invalid summary file format %s (files)", option.accum);
    ensure(total_tests == total_passed+total_failed, "invalid summary count in file %s", option.accum);

    // reset file
    rewind(fp);
  }
  else {
    // create the file
    fp = fopen(option.accum, "w+");
    ensure(fp, "failed to create or read summary file %s", option.accum);
    option.reset = 0;

    // init the time stamp
    tm = *localtime(&option.dat_t0);
    tz = 0;
  }

  double total_time = difftime(now, option.dat_t0) + tz - tm.tm_isdst*3600;
  total_ndtime  += (option.clk_t1 - option.clk_t0) / CLOCKS_PER_SEC;
  total_lines   += lines;
  total_numbers += numbers;

  fprintf(fp, " = tests summary (started at %04d.%02d.%02d %02d:%02d:%02d %+d)\n",
          tm.tm_year+1900, tm.tm_mon+1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, tm.tm_isdst);

  fprintf(fp, "   total diff time %6.2lf s  -  total lines %6ld  -  total numbers %8ld\n",
              total_ndtime, total_lines, total_numbers);

  fprintf(fp, "   total run  time %6.0f s  -  total files %6d  -  PASSED %4d  -  FAILED %4d\n",
              total_time, total_tests+total, total_passed+(total-failed), total_failed+failed);

  fclose(fp);
}

static const double pow10_tbl[2*99+1] = { 
  1e-99, 1e-98, 1e-97, 1e-96, 1e-95, 1e-94, 1e-93, 1e-92, 1e-91, 1e-90,
  1e-89, 1e-88, 1e-87, 1e-86, 1e-85, 1e-84, 1e-83, 1e-82, 1e-81, 1e-80,
  1e-79, 1e-78, 1e-77, 1e-76, 1e-75, 1e-74, 1e-73, 1e-72, 1e-71, 1e-70,
  1e-69, 1e-68, 1e-67, 1e-66, 1e-65, 1e-64, 1e-63, 1e-62, 1e-61, 1e-60,
  1e-59, 1e-58, 1e-57, 1e-56, 1e-55, 1e-54, 1e-53, 1e-52, 1e-51, 1e-50,
  1e-49, 1e-48, 1e-47, 1e-46, 1e-45, 1e-44, 1e-43, 1e-42, 1e-41, 1e-40,
  1e-39, 1e-38, 1e-37, 1e-36, 1e-35, 1e-34, 1e-33, 1e-32, 1e-31, 1e-30,
  1e-29, 1e-28, 1e-27, 1e-26, 1e-25, 1e-24, 1e-23, 1e-22, 1e-21, 1e-20,
  1e-19, 1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10,
  1e-09, 1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01,

  1e+00, 1e+01, 1e+02, 1e+03, 1e+04, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09,
  1e+10, 1e+11, 1e+12, 1e+13, 1e+14, 1e+15, 1e+16, 1e+17, 1e+18, 1e+19,
  1e+20, 1e+21, 1e+22, 1e+23, 1e+24, 1e+25, 1e+26, 1e+27, 1e+28, 1e+29,
  1e+30, 1e+31, 1e+32, 1e+33, 1e+34, 1e+35, 1e+36, 1e+37, 1e+38, 1e+39,
  1e+40, 1e+41, 1e+42, 1e+43, 1e+44, 1e+45, 1e+46, 1e+47, 1e+48, 1e+49,
  1e+50, 1e+51, 1e+52, 1e+53, 1e+54, 1e+55, 1e+56, 1e+57, 1e+58, 1e+59,
  1e+60, 1e+61, 1e+62, 1e+63, 1e+64, 1e+65, 1e+66, 1e+67, 1e+68, 1e+69,
  1e+70, 1e+71, 1e+72, 1e+73, 1e+74, 1e+75, 1e+76, 1e+77, 1e+78, 1e+79,
  1e+80, 1e+81, 1e+82, 1e+83, 1e+84, 1e+85, 1e+86, 1e+87, 1e+88, 1e+89,
  1e+90, 1e+91, 1e+92, 1e+93, 1e+94, 1e+95, 1e+96, 1e+97, 1e+98, 1e+99,
};

const double *const pow10_table99 = &pow10_tbl[99];

