// *****************************************************************************
/*!
  \file      src/RNG/RNGStack.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Stack of random number generators
  \details   This file defines class RNGStack, which implements various
    functionality related to registering and instantiating random number
    generators interfacing to multiple RNG libraries. Registration and
    instantiation use a random number generator factory, which is a std::map (an
    associative container), associating unique RNG keys to their constructor
    calls. For more details, see the in-code documentation of the constructor.
*/
// *****************************************************************************

#include <iterator>
#include <utility>

#include "NoWarning/charm.h"
#include "NoWarning/threefry.h"
#include "NoWarning/philox.h"

#include "Tags.h"
#include "Factory.h"
#include "Exception.h"
#include "RNGStack.h"
#include "RNGSSE.h"
#include "Options/RNGSSESeqLen.h"
#include "Random123.h"
#include "QuinoaConfig.h"

#ifdef HAS_RNGSSE2
  #include <gm19.h>
  #include <gm29.h>
  #include <gm31.h>
  #include <gm55.h>
  #include <gm61.h>
  #include <gq58x1.h>
  #include <gq58x3.h>
  #include <gq58x4.h>
  #include <mt19937.h>
  #include <lfsr113.h>
  #include <mrg32k3a.h>
#endif

#ifdef HAS_MKL
  #include "MKLRNG.h"
  #include "Options/MKLBetaMethod.h"
  #include "Options/MKLGaussianMethod.h"
  #include "Options/MKLUniformMethod.h"
#endif

using tk::RNGStack;

RNGStack::RNGStack(
                    #ifdef HAS_MKL
                    const tk::ctr::RNGMKLParameters& mklparam,
                    #endif
                    #ifdef HAS_RNGSSE2
                    const tk::ctr::RNGSSEParameters& rngsseparam,
                    #endif
                    const tk::ctr::RNGRandom123Parameters& r123param )
 : m_factory()
// *****************************************************************************
//  Constructor: register generators into factory for each supported library
//! \param[in] mklparam MKL RNG parameters to use to configure MKL RNGs
//! \param[in] rngsseparam RNGSSE RNG parameters to use to configure RNGSSE RNGs
//! \param[in] rngr123param Random123 RNG parameters to use to configure
//!   Random123 RNGs
//! \author J. Bakosi
// *****************************************************************************
{
  #ifdef HAS_MKL
  regMKL( CkNumPes(), mklparam );
  #endif
  #ifdef HAS_RNGSSE2
  regRNGSSE( CkNumPes(), rngsseparam );
  #endif
  regRandom123( CkNumPes(), r123param );
}

std::map< tk::ctr::RawRNGType, tk::RNG >
RNGStack::selected( const std::vector< tk::ctr::RNGType >& sel ) const
// *****************************************************************************
//  Instantiate selected RNGs from factory and place them in map
//! \param[in] sel Vector of selected RNGs to instantiate (selected by the user)
//! \return A std::map of keys and their associated instantiated RNG objects
//! \author  J. Bakosi
// *****************************************************************************
{
  using tk::ctr::RawRNGType;
  std::map< RawRNGType, tk::RNG > rng;
  for (const auto& s : sel) {
    const auto r = m_factory.find(s);
    if (r != end(m_factory))
      rng.emplace( static_cast< RawRNGType >( s ), r->second() );
    else Throw( "RNG not found in factory" );
  }
  return rng;
}

#ifdef HAS_MKL
void
RNGStack::regMKL( int nstreams, const tk::ctr::RNGMKLParameters& param )
// *****************************************************************************
//  Register MKL random number generators into factory
//! \details Note that registering these entries in the map does not
//!   invoke the constructors. The mapped value simply stores how the
//!   constructors should be invoked at a later time. At some point later,
//!   based on user input, we then instantiate only the RNGs (correctly
//!   configured by the pre-bound constructor arguments) that are requested by
//!   the user.
//! \param[in] nstreams Register MKL RNG using this many independent streams
//! \param[in] param MKL RNG parameters to use to configure the RNGs
//! \author J. Bakosi
// *****************************************************************************
{
  using tk::ctr::RNGType;
  using tk::ctr::MKLUniformMethodType;
  using tk::ctr::MKLGaussianMethodType;
  using tk::ctr::MKLBetaMethodType;
  using tag::uniform_method;
  using tag::gaussian_method;
  using tag::beta_method;

  // Defaults for MKL RNGs
  unsigned int s_def = 0;
  MKLUniformMethodType u_def = MKLUniformMethodType::STANDARD;
  MKLGaussianMethodType g_def = MKLGaussianMethodType::BOXMULLER;
  MKLBetaMethodType b_def = MKLBetaMethodType::CJA;

  tk::ctr::RNG opt;
  tk::ctr::MKLUniformMethod um_opt;
  tk::ctr::MKLGaussianMethod gm_opt;
  tk::ctr::MKLBetaMethod bm_opt;

  //! Lambda to register a MKL random number generator into factory
  auto regMKLRNG = [&]( RNGType rng ) {
    recordModel< tk::RNG, tk::MKLRNG >
      ( m_factory, rng,
        nstreams,
        opt.param( rng ),
        opt.param< tag::seed >( rng, s_def, param ),
        um_opt.param( opt.param< uniform_method >( rng, u_def, param ) ),
        gm_opt.param( opt.param< gaussian_method >( rng, g_def, param) ),
        bm_opt.param( opt.param< beta_method >( rng, b_def, param) ) );
  };

  // Register MKL RNGs
  regMKLRNG( RNGType::MKL_MCG31 );
  regMKLRNG( RNGType::MKL_R250 );
  regMKLRNG( RNGType::MKL_MRG32K3A );
  regMKLRNG( RNGType::MKL_MCG59 );
  regMKLRNG( RNGType::MKL_WH );
  regMKLRNG( RNGType::MKL_MT19937 );
  regMKLRNG( RNGType::MKL_MT2203 );
  regMKLRNG( RNGType::MKL_SFMT19937 );
  regMKLRNG( RNGType::MKL_SOBOL );
  regMKLRNG( RNGType::MKL_NIEDERR );
  //regMKLRNG( RNGType::MKL_IABSTRACT );
  //regMKLRNG( RNGType::MKL_DABSTRACT );
  //regMKLRNG( RNGType::MKL_SABSTRACT );
  regMKLRNG( RNGType::MKL_NONDETERM );
}
#endif

#ifdef HAS_RNGSSE2
void
RNGStack::regRNGSSE( int nstreams, const tk::ctr::RNGSSEParameters& param )
// *****************************************************************************
//  Register RNGSSE random number generators into factory
//! \details Note that registering these entries in the map does not
//!   invoke the constructors. The mapped value simply stores how the
//!   constructors should be invoked at a later time. At some point later,
//!   based on user input, we then instantiate only the RNGs (correctly
//!   configured by the pre-bound constructor arguments) that are requested by
//!   the user.
//! \param[in] nstreams Register RNGSSE RNG using this many independent streams
//! \param[in] param RNGSSE RNG parameters to use to configure the RNGs
//! \author J. Bakosi
// *****************************************************************************
{
  using tk::RNG;
  using tk::RNGSSE;
  using tk::ctr::RNGType;
  using tk::ctr::RNGSSESeqLenType;
  using tag::seqlen;

  tk::ctr::RNG opt;

  // Defaults for RNGSSE RNGs
  RNGSSESeqLenType l_def = RNGSSESeqLenType::SHORT;

  // Register RNGSSE RNGs
  recordModel< RNG, RNGSSE< gm19_state, unsigned, gm19_generate_ > >
    ( m_factory, RNGType::RNGSSE_GM19,
      nstreams,
      &gm19_init_sequence_ );

  recordModel< RNG, RNGSSE< gm29_state, unsigned, &gm29_generate_ > >
    ( m_factory, RNGType::RNGSSE_GM29,
      nstreams,
      &gm29_init_short_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GM29, l_def, param ),
      &gm29_init_long_sequence_,
      &gm29_init_medium_sequence_ );

  recordModel< RNG, RNGSSE< gm31_state, unsigned, gm31_generate_ > >
    ( m_factory, RNGType::RNGSSE_GM31,
      nstreams,
      &gm31_init_short_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GM31, l_def, param ),
      &gm31_init_long_sequence_,
      &gm31_init_medium_sequence_ );

  recordModel< RNG, RNGSSE< gm55_state, unsigned long long, gm55_generate_ > >
    ( m_factory, RNGType::RNGSSE_GM55,
      nstreams,
      &gm55_init_short_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GM55, l_def, param ),
      &gm55_init_long_sequence_ );

  recordModel< RNG, RNGSSE< gm61_state, unsigned long long, gm61_generate_ > >
    ( m_factory, RNGType::RNGSSE_GM61,
      nstreams,
      &gm61_init_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GM61, l_def, param ),
      &gm61_init_long_sequence_ );

  recordModel< RNG, RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ > >
    ( m_factory, RNGType::RNGSSE_GQ581,
      nstreams,
      &gq58x1_init_short_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GQ581, l_def, param ),
      &gq58x1_init_long_sequence_,
      &gq58x1_init_medium_sequence_ );

  recordModel< RNG, RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ > >
    ( m_factory, RNGType::RNGSSE_GQ583,
      nstreams,
      &gq58x3_init_short_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GQ583, l_def, param ),
      &gq58x3_init_long_sequence_,
      &gq58x3_init_medium_sequence_ );

  recordModel< RNG, RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ > >
    ( m_factory, RNGType::RNGSSE_GQ584,
      nstreams,
      &gq58x4_init_short_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GQ584, l_def, param ),
      &gq58x4_init_long_sequence_,
      &gq58x4_init_medium_sequence_ );

  recordModel< RNG,
               RNGSSE< mt19937_state, unsigned long long, mt19937_generate_ > >
    ( m_factory, RNGType::RNGSSE_MT19937,
      nstreams,
      &mt19937_init_sequence_ );

  recordModel< RNG,
               RNGSSE< lfsr113_state, unsigned long long, lfsr113_generate_ > >
    ( m_factory, RNGType::RNGSSE_LFSR113,
      nstreams,
      &lfsr113_init_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_LFSR113, l_def, param ),
      &lfsr113_init_long_sequence_ );

  recordModel< RNG,
               RNGSSE<mrg32k3a_state, unsigned long long, mrg32k3a_generate_> >
    ( m_factory, RNGType::RNGSSE_MRG32K3A,
      nstreams,
      &mrg32k3a_init_sequence_ );
}
#endif

void
RNGStack::regRandom123( int nstreams,
                        const tk::ctr::RNGRandom123Parameters& param )
// *****************************************************************************
//  Register Random123 random number generators into factory
//! \details Note that registering these entries in the map does not
//!   invoke the constructors. The mapped value simply stores how the
//!   constructors should be invoked at a later time. At some point later,
//!   based on user input, we then instantiate only the RNGs (correctly
//!   configured by the pre-bound constructor arguments) that are requested by
//!   the user.
//! \param[in] nstreams Register Randomer123 RNG using this many independent
//!   streams
//! \param[in] param Random123 RNG parameters to use to configure the RNGs
//! \author J. Bakosi
// *****************************************************************************
{
  using tk::ctr::RNGType;

  // Defaults for MKL RNGs
  uint32_t s_def = 0;

  tk::ctr::RNG opt;

  // Register Random123 RNGs
  recordModel< tk::RNG, tk::Random123< r123::Threefry2x64 > >
             ( m_factory, RNGType::R123_THREEFRY,
               nstreams,
               opt.param< tag::seed >( RNGType::R123_THREEFRY, s_def, param ) );

  recordModel< tk::RNG, tk::Random123< r123::Philox2x64 > >
             ( m_factory, RNGType::R123_PHILOX,
               nstreams,
               opt.param< tag::seed >( RNGType::R123_PHILOX, s_def, param ) );
}
