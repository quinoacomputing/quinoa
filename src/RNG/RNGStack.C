//******************************************************************************
/*!
  \file      src/Main/RNGStack.C
  \author    J. Bakosi
  \date      Sun 08 Jun 2014 01:47:39 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Stack of random number generators
  \details   Stack of random number generators
*/
//******************************************************************************

extern "C" {
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
}

#include <RNGStack.h>
#include <Factory.h>
#include <RNGSSE.h>

#ifdef HAS_MKL
  #include <MKLRNG.h>
#endif

using tk::RNGStack;

void
RNGStack::initFactory( tk::RNGFactory& factory,
                       int nthreads,
                       #ifdef HAS_MKL
                       const tk::ctr::RNGMKLParameters& mklparam,
                       #endif
                       const tk::ctr::RNGSSEParameters& rngsseparam )
//******************************************************************************
//  Register random number generators into factory for each supported library
//! \author  J. Bakosi
//******************************************************************************
{
  #ifdef HAS_MKL
  regMKL( factory, nthreads, mklparam );
  #endif
  regRNGSSE( factory, nthreads, rngsseparam );
}

std::map< tk::ctr::RawRNGType, tk::RNG >
RNGStack::createSelected( const tk::RNGFactory& factory,
                          const std::vector< tk::ctr::RNGType >& selected )
//******************************************************************************
//  Instantiate selected RNGs from factory and place them in map
//! \author  J. Bakosi
//******************************************************************************
{
  using tk::ctr::RawRNGType;
  std::map< RawRNGType, tk::RNG > rng;
  for (const auto& s : selected) {
    const auto r = factory.find(s);
    if (r != end(factory)) rng.emplace(static_cast<RawRNGType>(s), r->second());
    else Throw( tk::ExceptType::FATAL, "RNG not found in factory" );
  }
  return rng;
}

#ifdef HAS_MKL
void
RNGStack::regMKL( tk::RNGFactory& factory,
                  int nthreads,
                  const tk::ctr::RNGMKLParameters& param )
//******************************************************************************
//  Register MKL random number generators into factory
//! \author  J. Bakosi
//******************************************************************************
{
  using tk::ctr::RNGType;
  using tk::ctr::MKLUniformMethodType;
  using tk::ctr::MKLGaussianMethodType;
  using tk::tag::uniform_method;
  using tk::tag::gaussian_method;

  // Defaults for MKL RNGs
  unsigned int s_def = 0;
  MKLUniformMethodType u_def = MKLUniformMethodType::STANDARD;
  MKLGaussianMethodType g_def = MKLGaussianMethodType::BOXMULLER;

  tk::ctr::RNG opt;
  tk::ctr::MKLUniformMethod um_opt;
  tk::ctr::MKLGaussianMethod gm_opt;

  //! Lambda to register a MKL random number generator into factory
  auto regMKLRNG = [&]( RNGType rng ) {
    recordModel< tk::RNG, tk::MKLRNG >
      ( factory, rng,
        nthreads,
        opt.param( rng ),
        opt.param< tk::tag::seed >( rng, s_def, param ),
        um_opt.param( opt.param< uniform_method >( rng, u_def, param ) ),
        gm_opt.param( opt.param< gaussian_method >( rng, g_def, param) ) );
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

void
RNGStack::regRNGSSE( tk::RNGFactory& factory,
                     int nthreads,
                     const tk::ctr::RNGSSEParameters& param )
//******************************************************************************
//  Register RNGSSE random number generators into factory
//! \author  J. Bakosi
//******************************************************************************
{
  using tk::RNG;
  using tk::RNGSSE;
  using tk::ctr::RNGType;
  using tk::ctr::RNGSSESeqLenType;
  using tk::tag::seqlen;

  tk::ctr::RNG opt;

  // Defaults for RNGSSE RNGs
  RNGSSESeqLenType l_def = RNGSSESeqLenType::SHORT;

  // Register RNGSSE RNGs
  recordModel< RNG, RNGSSE< gm19_state, unsigned, gm19_generate_ > >
    ( factory, RNGType::RNGSSE_GM19,
      nthreads,
      &gm19_init_sequence_ );

  recordModel< RNG, RNGSSE< gm29_state, unsigned, &gm29_generate_ > >
    ( factory, RNGType::RNGSSE_GM29,
      nthreads,
      &gm29_init_short_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GM29, l_def, param ),
      &gm29_init_long_sequence_,
      &gm29_init_medium_sequence_ );

  recordModel< RNG, RNGSSE< gm31_state, unsigned, gm31_generate_ > >
    ( factory, RNGType::RNGSSE_GM31,
      nthreads,
      &gm31_init_short_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GM31, l_def, param ),
      &gm31_init_long_sequence_,
      &gm31_init_medium_sequence_ );

  recordModel< RNG, RNGSSE< gm55_state, unsigned long long, gm55_generate_ > >
    ( factory, RNGType::RNGSSE_GM55,
      nthreads,
      &gm55_init_short_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GM55, l_def, param ),
      &gm55_init_long_sequence_ );

  recordModel< RNG, RNGSSE< gm61_state, unsigned long long, gm61_generate_ > >
    ( factory, RNGType::RNGSSE_GM61,
      nthreads,
      &gm61_init_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GM61, l_def, param ),
      &gm61_init_long_sequence_ );

  recordModel< RNG, RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ > >
    ( factory, RNGType::RNGSSE_GQ581,
      nthreads,
      &gq58x1_init_short_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GQ581, l_def, param ),
      &gq58x1_init_long_sequence_,
      &gq58x1_init_medium_sequence_ );

  recordModel< RNG, RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ > >
    ( factory, RNGType::RNGSSE_GQ583,
      nthreads,
      &gq58x3_init_short_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GQ583, l_def, param ),
      &gq58x3_init_long_sequence_,
      &gq58x3_init_medium_sequence_ );

  recordModel< RNG, RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ > >
    ( factory, RNGType::RNGSSE_GQ584,
      nthreads,
      &gq58x4_init_short_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_GQ584, l_def, param ),
      &gq58x4_init_long_sequence_,
      &gq58x4_init_medium_sequence_ );

  recordModel< RNG,
               RNGSSE< mt19937_state, unsigned long long, mt19937_generate_ > >
    ( factory, RNGType::RNGSSE_MT19937,
      nthreads,
      &mt19937_init_sequence_ );

  recordModel< RNG,
               RNGSSE< lfsr113_state, unsigned long long, lfsr113_generate_ > >
    ( factory, RNGType::RNGSSE_LFSR113,
      nthreads,
      &lfsr113_init_sequence_,
      opt.param< seqlen >( RNGType::RNGSSE_LFSR113, l_def, param ),
      &lfsr113_init_long_sequence_ );

  recordModel< RNG,
               RNGSSE<mrg32k3a_state, unsigned long long, mrg32k3a_generate_> >
    ( factory, RNGType::RNGSSE_MRG32K3A,
      nthreads,
      &mrg32k3a_init_sequence_ );
}
