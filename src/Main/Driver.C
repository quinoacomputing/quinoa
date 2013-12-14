//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Sat 14 Dec 2013 11:43:47 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base
  \details   Driver base
*/
//******************************************************************************

extern "C" {
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
  regMKL( factory, opt, reg, nthreads, mklparam );
  regRNGSSE( factory, opt, reg, nthreads, rngsseparam );
}

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
         opt.param< quinoa::ctr::seed >( rng, 0, param ),
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

  // Register RNGSSE RNGs
  add< quinoa::RNGSSE< mrg32k3a_state,
                       mrg32k3a_init_sequence_,
                       mrg32k3a_generate_ > >
     ( factory, reg, opt, RNGType::RNGSSE_MRG32K3A,
       nthreads,
       opt.param< quinoa::ctr::seed >( RNGType::RNGSSE_MRG32K3A, 0, param ) );
}
