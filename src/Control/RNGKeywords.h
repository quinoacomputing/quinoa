//******************************************************************************
/*!
  \file      src/Control/MKLRNGKeywords.h
  \author    J. Bakosi
  \date      Sun 22 Jun 2014 10:34:30 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     RNG keywords for Intel's MKL
  \details   Random number generator selector keywords for those generators
  available in Intel's Math Kernel Library's (MKL) Vector Statistical Library
  (VSL). The keywords are defined by specializing struct 'keyword', defined in
  Control/Keyword.h. Introducing a new keyword requires a more human readable
  (bust still short) name as well as a short, few-line, help-like description.
*/
//******************************************************************************
#ifndef Keywords
#error "MKLRNGKeywords.h should only be included within a *Keywords.h"
#endif

namespace tk {
namespace kw {

using namespace pegtl::ascii;

// Keyword 'seed'
struct seed_info {
  static const char* name() { return "seed"; }
  static const char* help() { return
    "This keyword is used to denote the random number generator seed, if "
    "applicable. Example:\n"
    "\trngsse_gm55\n"
    "\t  seed 1234\n"
    "\tend";
  }
};
using seed = keyword< seed_info, s,e,e,d >;

} // kw::
} // tk::
