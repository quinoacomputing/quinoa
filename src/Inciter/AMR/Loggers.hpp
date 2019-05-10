#ifndef Loggers_h
#define Loggers_h

#include <iostream>

//#define ENABLE_TRACE 1

#ifdef ENABLE_TRACE
#define trace_out std::cout << "__TRACE: "
#else
#define trace_out while(0) std::cout
#endif /* ENABLE_TRACE */

#ifdef ENABLE_DEBUG
#define debug_out std::cout << "__DEBUG: "
#else
#define debug_out while(0) std::cout
#endif /* ENABLE_DEBUG */

#endif
