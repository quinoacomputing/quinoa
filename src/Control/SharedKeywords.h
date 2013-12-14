//******************************************************************************
/*!
  \file      src/Control/SharedKeywords.h
  \author    J. Bakosi
  \date      Fri 13 Dec 2013 06:23:24 PM MST
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

// RNGSSE's keywords
#include <RNGSSEKeywords.h>

#endif // SharedKeywords_h
