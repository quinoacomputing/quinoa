//******************************************************************************
/*!
  \file      src/Control/SharedKeywords.h
  \author    J. Bakosi
  \date      Sat 12 Jul 2014 06:39:43 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Shared keywords are shared among different executables
  \details   Shared keywords are shared among different executables
*/
//******************************************************************************
#ifndef SharedKeywords_h
#define SharedKeywords_h

// Base keywords recognized by all input deck parsers
#include <InputDeckBaseKeywords.h>

// Keywords common to random number generators
#include <RNGKeywords.h>

// Intel's MKL's RNG keywords
#include <MKLRNGKeywords.h>

// RNGSSE's keywords
#include <RNGSSEKeywords.h>

#endif // SharedKeywords_h
