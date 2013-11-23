//******************************************************************************
/*!
  \file      src/Control/MKLRNGKeywords.h
  \author    J. Bakosi
  \date      Fri 22 Nov 2013 05:19:33 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNG keywords for Intel's MKL
  \details   Random number generator selector keywords for those generators
  available in Intel's Math Kernel Library's (MKL) Vector Statistical Library
  (VSL). The keywords are defined by specializing struct 'keyword', defined in
  Control/Keyword.h. Introducing a new keyword requires a more human readable
  (bust still short) name as well as a short, few-line, help-like description.
*/
//******************************************************************************
#ifndef Keywords
#error "RNGKeywords.h should only be included within a *Keywords.h"
#endif

namespace tk {
namespace kw {

using namespace pegtl::ascii;

// Random number generator seed
using seed = keyword<kw::undefined_info, s,e,e,d >;

} // kw::
} // tk::
