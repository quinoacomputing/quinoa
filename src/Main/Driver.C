//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Thu 21 Nov 2013 03:15:10 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base
  \details   Driver base
*/
//******************************************************************************

#include <Driver.h>
#include <MKLRNG.h>

using tk::Driver;

void
Driver::initRNGFactory( tk::RNGFactory& factory,
                        const quinoa::ctr::RNG& opt,
                        std::list< quinoa::ctr::RNGType >& reg,
                        int nthreads,
                        const quinoa::ctr::MKLRNGParameters& mklparam )
//******************************************************************************
//  Register random number generators into factory
//! \author  J. Bakosi
//******************************************************************************
{
  using quinoa::ctr::RNGType;
  using quinoa::ctr::seed;
  using quinoa::ctr::uniform_method;

  quinoa::ctr::MKLUniformMethod um_opt;
  quinoa::ctr::MKLGaussianMethod gm_opt;

  //! Lambda to register a MKL random number generator into factory
  auto regMKLRNG = [&]( RNGType rng ) {
    add< quinoa::MKLRNG >
       ( factory, reg, opt, rng,
         nthreads,
         opt.param( rng ),
         opt.mkl_seed( rng, mklparam ),
         um_opt.param( opt.mkl_uniform_method( rng, mklparam ) ),
         gm_opt.param( opt.mkl_gaussian_method( rng, mklparam) ) );
  };

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
