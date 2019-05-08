// *****************************************************************************
/*!
  \file      src/Control/RNGParam.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Types for storing parameters of random number generators
  \details   Types for storing parameters of random number generators.
*/
// *****************************************************************************
#ifndef RNGParam_h
#define RNGParam_h

#include "TaggedTuple.hpp"
#include "Tags.hpp"
#include "Options/RNG.hpp"
#include "Options/RNGSSESeqLen.hpp"
#include "QuinoaConfig.hpp"

#ifdef HAS_MKL
  #include "Options/MKLUniformMethod.hpp"
  #include "Options/MKLGaussianMethod.hpp"
  #include "Options/MKLGaussianMVMethod.hpp"
  #include "Options/MKLBetaMethod.hpp"
  #include "Options/MKLGammaMethod.hpp"
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
  tag::seed,              kw::seed::info::expect::type, //!< seed
  tag::uniform_method,    MKLUniformMethodType,         //!< uniform method type
  tag::gaussian_method,   MKLGaussianMethodType,        //!< Gaussian method type
  //! multi-variate Gaussian method type
  tag::gaussianmv_method, MKLGaussianMVMethodType,
  tag::beta_method,       MKLBetaMethodType,            //!< beta method type
  tag::gamma_method,      MKLGammaMethodType            //!< gamma method type
>;
//! MKL RNG parameters bundle associating RNG types and their parameters
using RNGMKLParameters = std::map< RNGType, RNGMKLParam >;
#endif

} // ctr::
} // tk::

#endif // RNGParam_h
