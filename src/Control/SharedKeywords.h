//******************************************************************************
/*!
  \file      src/Control/SharedKeywords.h
  \author    J. Bakosi
  \date      Fri 22 Nov 2013 05:18:52 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Shared keywords are shared among different executables
  \details   Shared keywords are shared among different executables
*/
//******************************************************************************
#ifndef SharedKeywords_h
#define SharedKeywords_h

// Base keywords recognized by all input deck parsers
#include <BaseKeywords.h>

// Keywords common to random number generators
#include <RNGKeywords.h>

// Intel's MKL's RNG keywords
#include <MKLRNGKeywords.h>

#endif // SharedKeywords_h
