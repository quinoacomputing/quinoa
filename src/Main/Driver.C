//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Fri 27 Dec 2013 07:48:45 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base
  \details   Driver base
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

#include <Driver.h>
#include <MKLRNG.h>

using tk::Driver;

void
Driver::initRNGFactory( tk::RNGFactory& factory,
                        const quinoa::ctr::RNG& opt,
                        std::list< quinoa::ctr::RNGType >& reg,
                        int nthreads,
                        const quinoa::ctr::MKLRNGParameters& mklparam,
                        const quinoa::ctr::RNGSSEParameters& rngsseparam )
//******************************************************************************
//  Register random number generators into factory
//! \author  J. Bakosi
//******************************************************************************
{
  #ifdef HAS_MKL
  regMKL( factory, opt, reg, nthreads, mklparam );
  #endif
  regRNGSSE( factory, opt, reg, nthreads, rngsseparam );
}

#ifdef HAS_MKL
void
Driver::regMKL( tk::RNGFactory& factory,
                const quinoa::ctr::RNG& opt,
                std::list< quinoa::ctr::RNGType >& reg,
                int nthreads,
                const quinoa::ctr::MKLRNGParameters& param )
//******************************************************************************
//  Register MKL random number generators into factory
//! \author  J. Bakosi
//******************************************************************************
{
  using quinoa::ctr::RNGType;
  using quinoa::ctr::uniform_method;
  using quinoa::ctr::gaussian_method;
  using quinoa::ctr::MKLUniformMethodType;
  using quinoa::ctr::MKLGaussianMethodType;

  // Defaults for MKL RNGs
  unsigned int s_def = 0;
  MKLUniformMethodType u_def = MKLUniformMethodType::STANDARD;
  MKLGaussianMethodType g_def = MKLGaussianMethodType::BOXMULLER;

  quinoa::ctr::MKLUniformMethod um_opt;
  quinoa::ctr::MKLGaussianMethod gm_opt;

  //! Lambda to register a MKL random number generator into factory
  auto regMKLRNG = [&]( RNGType rng ) {
    add< quinoa::MKLRNG >
       ( factory, reg, opt, rng,
         nthreads,
         opt.param( rng ),
         opt.param< quinoa::ctr::seed >( rng, s_def, param ),
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
  regMKLRNG( RNGType::MKL_IABSTRACT );
  regMKLRNG( RNGType::MKL_DABSTRACT );
  regMKLRNG( RNGType::MKL_SABSTRACT );
  regMKLRNG( RNGType::MKL_NONDETERM );
}
#endif

void
Driver::regRNGSSE( tk::RNGFactory& factory,
                   const quinoa::ctr::RNG& opt,
                   std::list< quinoa::ctr::RNGType >& reg,
                   int nthreads,
                   const quinoa::ctr::RNGSSEParameters& param )
//******************************************************************************
//  Register RNGSSE random number generators into factory
//! \author  J. Bakosi
//******************************************************************************
{
  using quinoa::ctr::RNGType;
//  using quinoa::ctr::seqlen;
//  using quinoa::ctr::RNGSSESeqLenType;

  // Defaults for RNGSSE RNGs
  unsigned int s_def = 0;
//  RNGSSESeqLenType l_def = RNGSSESeqLenType::SHORT;

  // Register RNGSSE RNGs

  // Get sequence length
//   if ( opt.param< seqlen >( RNGType::RNGSSE_GM19, l_def, param ) == l_def ) {
//     std::cout << "\n>>> WARNING: RNGSSE_GM19 does not take a sequence length parameter. Remove the 'seqlen' parameter from the rngsse_gm19 ... end block from the control file.\n\n";
//   }
  add< quinoa::RNGSSE< gm19_state,
                       unsigned int,
                       gm19_init_sequence_,
                       gm19_generate_ > >
     ( factory, reg, opt, RNGType::RNGSSE_GM19,
       nthreads,
       opt.param< quinoa::ctr::seed >( RNGType::RNGSSE_GM19, s_def, param ) );

  add< quinoa::RNGSSE< gm29_state,
                       unsigned,
                       gm29_init_short_sequence_,
                       gm29_generate_ > >
     ( factory, reg, opt, RNGType::RNGSSE_GM29,
       nthreads,
       opt.param< quinoa::ctr::seed >( RNGType::RNGSSE_GM29, s_def, param ) );

  add< quinoa::RNGSSE< gm31_state,
                       unsigned,
                       gm31_init_short_sequence_,
                       gm31_generate_ > >
     ( factory, reg, opt, RNGType::RNGSSE_GM31,
       nthreads,
       opt.param< quinoa::ctr::seed >( RNGType::RNGSSE_GM31, s_def, param ) );

  add< quinoa::RNGSSE< gm55_state,
                       unsigned long long,
                       gm55_init_short_sequence_,
                       gm55_generate_ > >
     ( factory, reg, opt, RNGType::RNGSSE_GM55,
       nthreads,
       opt.param< quinoa::ctr::seed >( RNGType::RNGSSE_GM55, s_def, param ) );

  add< quinoa::RNGSSE< gm61_state,
                       unsigned long long,
                       gm61_init_sequence_,
                       gm61_generate_ > >
     ( factory, reg, opt, RNGType::RNGSSE_GM61,
       nthreads,
       opt.param< quinoa::ctr::seed >( RNGType::RNGSSE_GM61, s_def, param ) );

  add< quinoa::RNGSSE< gq58x1_state,
                       unsigned,
                       gq58x1_init_short_sequence_,
                       gq58x1_generate_ > >
     ( factory, reg, opt, RNGType::RNGSSE_GQ581,
       nthreads,
       opt.param< quinoa::ctr::seed >( RNGType::RNGSSE_GQ581, s_def, param ) );

  add< quinoa::RNGSSE< gq58x3_state,
                       unsigned,
                       gq58x3_init_short_sequence_,
                       gq58x3_generate_ > >
     ( factory, reg, opt, RNGType::RNGSSE_GQ583,
       nthreads,
       opt.param< quinoa::ctr::seed >( RNGType::RNGSSE_GQ583, s_def, param ) );

  add< quinoa::RNGSSE< gq58x4_state,
                       unsigned,
                       gq58x4_init_short_sequence_,
                       gq58x4_generate_ > >
     ( factory, reg, opt, RNGType::RNGSSE_GQ584,
       nthreads,
       opt.param< quinoa::ctr::seed >( RNGType::RNGSSE_GQ584, s_def, param ) );

  add< quinoa::RNGSSE< mt19937_state,
                       unsigned long long,
                       mt19937_init_sequence_,
                       mt19937_generate_ > >
     ( factory, reg, opt, RNGType::RNGSSE_MT19937,
       nthreads,
       opt.param< quinoa::ctr::seed >( RNGType::RNGSSE_MT19937, s_def, param ) );

  add< quinoa::RNGSSE< lfsr113_state,
                       unsigned long long,
                       lfsr113_init_sequence_,
                       lfsr113_generate_ > >
     ( factory, reg, opt, RNGType::RNGSSE_LFSR113,
       nthreads,
       opt.param< quinoa::ctr::seed >( RNGType::RNGSSE_LFSR113, s_def, param ) );

  add< quinoa::RNGSSE< mrg32k3a_state,
                       unsigned long long,
                       mrg32k3a_init_sequence_,
                       mrg32k3a_generate_ > >
     ( factory, reg, opt, RNGType::RNGSSE_MRG32K3A,
       nthreads,
       opt.param< quinoa::ctr::seed >( RNGType::RNGSSE_MRG32K3A, s_def, param ) );
}
