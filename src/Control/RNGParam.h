// *****************************************************************************
/*!
  \file      src/Control/RNGParam.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Types for storing parameters of random number generators
  \details   Types for storing parameters of random number generators.
*/
// *****************************************************************************
#ifndef RNGParam_h
#define RNGParam_h

#include "TaggedTuple.h"
#include "Tags.h"
#include "Options/RNG.h"
#include "Options/RNGSSESeqLen.h"
#include "QuinoaConfig.h"

#ifdef HAS_MKL
  #include "Options/MKLUniformMethod.h"
  #include "Options/MKLGaussianMethod.h"
  #include "Options/MKLBetaMethod.h"
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

//! Random123 random number generator parameters storage
using RNGRandom123Param = tk::tuple::tagged_tuple<
  tag::seed,          kw::seed::info::expect::type   //!< seed
>;
//! Random123 parameters bundle associating RNG types and their parameters
using RNGRandom123Parameters = std::map< RNGType, RNGRandom123Param >;

#ifdef HAS_MKL
//! MKL random number generator parameters storage
using RNGMKLParam = tk::tuple::tagged_tuple<
  tag::seed,            kw::seed::info::expect::type, //!< seed
  tag::uniform_method,  MKLUniformMethodType,         //!< uniform method type
  tag::gaussian_method, MKLGaussianMethodType,        //!< Gaussian method type
  tag::beta_method,     MKLBetaMethodType             //!< beta method type
>;
//! MKL RNG parameters bundle associating RNG types and their parameters
using RNGMKLParameters = std::map< RNGType, RNGMKLParam >;
#endif

} // ctr::
} // tk::

#endif // RNGParam_h
