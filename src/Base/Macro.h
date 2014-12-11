//******************************************************************************
/*!
  \file      src/Base/Macro.h
  \author    J. Bakosi
  \date      Thu 11 Dec 2014 07:54:11 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Macro definitions
  \details   Macro definitions for various utility functionality.
*/
//******************************************************************************
#ifndef Macro_h
#define Macro_h

#include <sys/time.h>
#include <ctime>
#include <cstdlib>

namespace tk {

//! This macro can be used to suppress compiler warning on unused variable
#define IGNORE(expr) (static_cast<void>(expr))

//! Start-time macro for fine-grained profiling. Put this in the beginning of
//! the section of code to be profiled.
//! \author J. Bakosi
#define STARTTIME \
struct timeval START_TIME, END_TIME; \
int total_usecs; \
gettimeofday(&START_TIME, (struct timezone*)0);

//! End-time macro for fine-grained profiling. Put this at the end of the
//! section of code to be profiled.
//! \author J. Bakosi
#define ENDTIME \
gettimeofday(&END_TIME, (struct timezone*)0); \
total_usecs = (END_TIME.tv_sec-START_TIME.tv_sec) * 1000000 + (END_TIME.tv_usec-START_TIME.tv_usec); \
printf("Total time was %d uSec.\n", total_usecs);

} // tk::

#endif // Macro_h
