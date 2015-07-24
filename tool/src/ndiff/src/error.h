#ifndef ERROR_H
#define ERROR_H

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
     manage error level and log error messages
 
 o---------------------------------------------------------------------o
*/

// ----- constants

enum {
  trace_level, debug_level, inform_level, warning_level, error_level, abort_level
};

// ----- interface

#define trace(...)        logmsg_if(trace_level  , __VA_ARGS__)
#define debug(...)        logmsg_if(debug_level  , __VA_ARGS__)
#define inform(...)       logmsg_if(inform_level , __VA_ARGS__)
#define warning(...)      logmsg_if(warning_level, __VA_ARGS__)
#define error(...)        logmsg_if(error_level  , __VA_ARGS__) // exit
#define abort(...)        logmsg_if(abort_level  , __VA_ARGS__) // abort

#define ensure(cond, ...) ((void)(!(cond) && (error(__VA_ARGS__),0)))

// low-level
void logmsg(unsigned level, const char *file, int line, const char *fmt, ...);

// ----- configuration

struct logmsg_config {
  unsigned level;  // level >= logmsg_config.level are active,  default warning_level
  int      locate; // print (file, line) information,           default zero
  int      flush;  // flush stdout before output if non-zero,   default non-zero
};

extern struct logmsg_config logmsg_config;

// ----- enable/disable call to function logmsg

// compile-time logmsg level, override with -DLOGMSG_LEVEL=xxx_level or in files
#ifndef NDEBUG
#define LOGMSG_LEVEL trace_level
#else
#define LOGMSG_LEVEL inform_level
#endif

// runtime logmsg level, configure with logmsg_config.level=xxx_level
#define logmsg_if(LEVEL, ...) \
  ((void)(LEVEL >= LOGMSG_LEVEL &&   /* compile-time */ \
          LEVEL >= logmsg_config.level && /* runtime */ \
          (logmsg(LEVEL, __FILE__, __LINE__, __VA_ARGS__),0)))

#endif
