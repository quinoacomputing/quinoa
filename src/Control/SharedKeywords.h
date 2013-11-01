//******************************************************************************
/*!
  \file      src/Control/SharedKeywords.h
  \author    J. Bakosi
  \date      Thu 31 Oct 2013 08:21:16 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Shared keywords are shared among different executables
  \details   Shared keywords are shared among different executables
*/
//******************************************************************************
#ifndef SharedKeywords_h
#define SharedKeywords_h

namespace tk {
namespace kw {

// Base keywords recognized by all input deck parsers
#include <BaseKeywords.h>

// Keywords common to random number generators
#include <RNGKeywords.h>

// Intel's MKL's RNG keywords
#include <MKLRNGKeywords.h>

} // kw::
} // tk::

#endif // SharedKeywords_h
