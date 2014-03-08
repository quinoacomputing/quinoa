//******************************************************************************
/*!
  \file      src/Control/tkTypes.h
  \author    J. Bakosi
  \date      Sat 08 Mar 2014 06:38:34 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for tk control
  \details   Types for tk control
*/
//******************************************************************************
#ifndef tkTypes_h
#define tkTypes_h

#include <TaggedTuple.h>
#include <tkTags.h>
#include <Options/RNG.h>
#include <Options/RNGSSESeqLen.h>

#ifdef HAS_MKL
#include <Options/MKLUniformMethod.h>
#include <Options/MKLGaussianMethod.h>
#endif

namespace tk {
namespace ctr {

//! RNGSSE random number generator parameters storage
using RNGSSEParam = tuple::tagged_tuple<
  tag::seed,          unsigned int,              //!< seed
  tag::seqlen,        RNGSSESeqLenType           //!< sequence length type
>;
//! RNGSSE parameters bundle
using RNGSSEParameters = std::map< RNGType, RNGSSEParam >;

#ifdef HAS_MKL
//! MKL random number generator parameters storage
using MKLRNGParam = tuple::tagged_tuple<
  tag::seed,             unsigned int,              //!< seed
  tag::uniform_method,   MKLUniformMethodType,      //!< uniform method type
  tag::gaussian_method,  MKLGaussianMethodType      //!< Gaussian method type
>;
//! MKL RNG parameters bundle
using MKLRNGParameters = std::map< RNGType, MKLRNGParam >;
#endif

} // ctr::
} // tk::

#endif // tkTypes_h
