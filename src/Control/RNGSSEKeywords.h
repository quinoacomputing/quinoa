//******************************************************************************
/*!
  \file      src/Control/RNGSSEKeywords.h
  \author    J. Bakosi
  \date      Fri 13 Dec 2013 07:30:44 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNG keywords for RNGSSE lib
  \details   Random number generator selector keywords for those generators
  available in the RNGSSE library's. The keywords are defined by specializing
  struct 'keyword', defined in Control/Keyword.h. Introducing a new keyword
  requires a more human readable (but still short) name as well as a short,
  few-line, help-like description.
*/
//******************************************************************************
#ifndef Keywords
#error "RNGSSEKeywords.h should only be included within a *Keywords.h"
#endif

namespace tk {
namespace kw {

using namespace pegtl::ascii;

// Keyword 'rngsse_mrg32k3a'
struct rngsse_mrg32k3a_info {
  static const char* name() { return "RNGSSE MRG32K3A"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's MRG32K3A generator, a "
    "combined multiple recursive random number generator with two components "
    "of order 3";
  }
};
using rngsse_mrg32k3a =
  keyword< rngsse_mrg32k3a_info, r,n,g,s,s,e,'_',m,r,g,'3','2',k,'3',a >;

} // kw::
} // tk::
