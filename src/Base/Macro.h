//******************************************************************************
/*!
  \file      src/Base/Macro.h
  \author    J. Bakosi
  \date      Mon Oct  7 08:08:07 2013
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Macro definitions
  \details   Macro definitions
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

//! Start-time macro for fine-grained profiling
#define STARTTIME \
struct timeval START_TIME, END_TIME; \
int total_usecs; \
gettimeofday(&START_TIME, (struct timezone*)0);

//! End-time macro for fine-grained profiling
#define ENDTIME \
gettimeofday(&END_TIME, (struct timezone*)0); \
total_usecs = (END_TIME.tv_sec-START_TIME.tv_sec) * 1000000 + (END_TIME.tv_usec-START_TIME.tv_usec); \
printf("Total time was %d uSec.\n", total_usecs);

} // tk::

#endif // Macro_h
