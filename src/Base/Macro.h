// *****************************************************************************
/*!
  \file      src/Base/Macro.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Macro definitions
  \details   Macro definitions for various utility functionality.
*/
// *****************************************************************************
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

//! Macro to detect strictly gcc.
//! \details __GNUC__ and __GNUG__ were intended to indicate the GNU compilers.
//! However, they're also defined by Clang/LLVM and Intel compilers to indicate
//! compatibility. This macro can be used to detect strictly gcc and not clang
//! or icc.
#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
  #define STRICT_GNUC
#endif

} // tk::

#endif // Macro_h
