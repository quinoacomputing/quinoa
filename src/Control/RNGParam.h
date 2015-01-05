//******************************************************************************
/*!
  \file      src/Control/RNGParam.h
  \author    J. Bakosi
  \date      Wed 14 Jan 2015 01:40:28 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Types for storing parameters of random number generators
  \details   Types for storing parameters of random number generators.
*/
//******************************************************************************
#ifndef RNGParam_h
#define RNGParam_h

#include <TaggedTuple.h>
#include <Tags.h>
#include <Options/RNG.h>
#include <Options/RNGSSESeqLen.h>

#ifdef HAS_MKL
#include <Options/MKLUniformMethod.h>
#include <Options/MKLGaussianMethod.h>
#endif

namespace tk {
namespace ctr {

//! RNGSSE random number generator parameters storage
using RNGSSEParam = tk::tuple::tagged_tuple<
  tag::seed,          kw::seed::info::expect::type,  //!< seed
  tag::seqlen,        RNGSSESeqLenType               //!< sequence length type
>;
//! RNGSSE parameters bundle associating RNG types and their parameters
using RNGSSEParameters = std::map< RNGType, RNGSSEParam >;

#ifdef HAS_MKL
//! MKL random number generator parameters storage
using RNGMKLParam = tk::tuple::tagged_tuple<
  tag::seed,            kw::seed::info::expect::type, //!< seed
  tag::uniform_method,  MKLUniformMethodType,         //!< uniform method type
  tag::gaussian_method, MKLGaussianMethodType         //!< Gaussian method type
>;
//! MKL RNG parameters bundle associating RNG types and their parameters
using RNGMKLParameters = std::map< RNGType, RNGMKLParam >;
#endif

} // ctr::
} // tk::

#endif // RNGParam_h
