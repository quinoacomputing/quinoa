//******************************************************************************
/*!
  \file      src/RNG/TestU01Suite.C
  \author    J. Bakosi
  \date      Wed 04 Dec 2013 12:13:17 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 suite
  \details   TestU01 suite
*/
//******************************************************************************

extern "C" {
  #include <smarsa.h>
  #include <sknuth.h>
  #include <svaria.h>
}

#include <TestU01Suite.h>

namespace rngtest {

std::vector< std::unique_ptr< tk::RNG > > g_rng;  //!< Global pointers to RNGs
int g_tid;                                        //!< Global thread id
#ifdef _OPENMP
#pragma omp threadprivate(g_tid)
#endif

template< int id >
static double uniform(void*, void*)
//******************************************************************************
//  TestU01 uniform RNG wrapper
//! \author  J. Bakosi
//******************************************************************************
{
  double r;
  g_rng[id]->uniform( g_tid, 1, &r );
  return r;
}

template< int id >
static unsigned long uniform_bits(void*, void*)
//******************************************************************************
//  TestU01 uniform RNG bits wrapper
//! \author  J. Bakosi
//******************************************************************************
{
  double r;
  g_rng[id]->uniform( g_tid, 1, &r );
  return static_cast<unsigned long>(r * unif01_NORM32);
}

} // rngtest::

using rngtest::TestU01Suite;
using Pvals = rngtest::StatTest::Pvals;

TestU01Suite::TestU01Suite( const Base& base ) : Battery(base)
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  using quinoa::ctr::RNGType;

  // I know this is lame for at least two reasons:
  //   * The order is important, and must start from zero
  //   * The ids must be literals as this is done at compile-time
  // This is due to two opposing requirements:
  //   * The wrappers, uniform() and uniform_bits(), this way can be templates
  //   * They must be in global scope as they are passed to TestU01
  addRNG< 0>( RNGType::MKL_MCG31,     uniform< 0>, uniform_bits< 0> );
  addRNG< 1>( RNGType::MKL_R250,      uniform< 1>, uniform_bits< 1> );
  addRNG< 2>( RNGType::MKL_MRG32K3A,  uniform< 2>, uniform_bits< 2> );
  addRNG< 3>( RNGType::MKL_MCG59,     uniform< 3>, uniform_bits< 3> );
  addRNG< 4>( RNGType::MKL_WH,        uniform< 4>, uniform_bits< 4> );
  addRNG< 5>( RNGType::MKL_MT19937,   uniform< 5>, uniform_bits< 5> );
  addRNG< 6>( RNGType::MKL_MT2203,    uniform< 6>, uniform_bits< 6> );
  addRNG< 7>( RNGType::MKL_SFMT19937, uniform< 7>, uniform_bits< 7> );
  addRNG< 8>( RNGType::MKL_SOBOL,     uniform< 8>, uniform_bits< 8> );
  addRNG< 9>( RNGType::MKL_NIEDERR,   uniform< 9>, uniform_bits< 9> );
  //addRNG<10>( RNGType::MKL_IABSTRACT, uniform<10>, uniform_bits<10> );
  //addRNG<11>( RNGType::MKL_DABSTRACT, uniform<11>, uniform_bits<11> );
  //addRNG<12>( RNGType::MKL_SABSTRACT, uniform<12>, uniform_bits<12> );
  //addRNG<13>( RNGType::MKL_NONDETERM, uniform<13>, uniform_bits<13> );
}

Pvals
TestU01Suite::BirthdaySpacings( unif01_Gen* gen, sres_Poisson* res )
//******************************************************************************
//  Run Marsdaglia's BirthdaySpacings test
//! \author  J. Bakosi
//******************************************************************************
{
#ifdef USE_LONGLONG
  smarsa_BirthdaySpacings( gen, res, 1, 5 * MILLION, 0, 1073741824, 2, 1 );
#else
  smarsa_BirthdaySpacings( gen, res, 10, MILLION / 2, 0, 67108864, 2, 1 );
#endif
  return Pvals( { res->pVal2 } );
}

Pvals
TestU01Suite::Collision( unif01_Gen* gen, sknuth_Res2* res )
//******************************************************************************
//  Run Knuth's Collision test
//! \author  J. Bakosi
//******************************************************************************
{
  sknuth_Collision( gen, res, 1, 5 * MILLION, 0, 65536, 2 );
  return Pvals( { res->Pois->pVal2 } );
}

Pvals
TestU01Suite::Gap( unif01_Gen* gen, sres_Chi2* res )
//******************************************************************************
//  Run Knuth's Gap test
//! \author  J. Bakosi
//******************************************************************************
{
  sknuth_Gap( gen, res, 1, MILLION / 5, 22, 0.0, .00390625 );
  return Pvals( { res->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::SimpPoker( unif01_Gen* gen, sres_Chi2* res )
//******************************************************************************
//  Run Knuth's Simplified Poker test
//! \author  J. Bakosi
//******************************************************************************
{
  sknuth_SimpPoker( gen, res, 1, 2 * MILLION / 5, 24, 64, 64 );
  return Pvals( { res->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::CouponCollector( unif01_Gen* gen, sres_Chi2* res )
//******************************************************************************
//  Run Knuth's Coupon Collector test
//! \author  J. Bakosi
//******************************************************************************
{
  sknuth_CouponCollector( gen, res, 1, MILLION / 2, 26, 16 );
  return Pvals( { res->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::MaxOft( unif01_Gen* gen, sknuth_Res1* res)
//******************************************************************************
//  Run Knuth's Maximum-of-t test
//! \author  J. Bakosi
//******************************************************************************
{
  sknuth_MaxOft( gen, res, 1, 2 * MILLION, 0, MILLION / 10, 6 );
  return Pvals( { res->Chi->pVal2[gofw_Mean],
                  res->Bas->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::WeightDistrib( unif01_Gen* gen, sres_Chi2* res)
//******************************************************************************
//  Run Matsumoto's Weight Distribution test
//! \author  J. Bakosi
//******************************************************************************
{
  svaria_WeightDistrib( gen, res, 1, MILLION / 5, 27, 256, 0.0, 0.125 );
  return Pvals( { res->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::MatrixRank( unif01_Gen* gen, sres_Chi2* res )
//******************************************************************************
//  Run Marsdaglia's Matrix Rank test
//! \author  J. Bakosi
//******************************************************************************
{
  smarsa_MatrixRank( gen, res, 1, 20 * THOUSAND, 20, 10, 60, 60 );
  return Pvals( { res->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::HammingIndep( unif01_Gen* gen, sstring_Res* res )
//******************************************************************************
//  Run L'Ecuyer's Hamming Independence test
//! \author  J. Bakosi
//******************************************************************************
{
  sstring_HammingIndep( gen, res, 1, MILLION/2, 20, 10, 300, 0 );
  return Pvals( { res->Bas->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::RandomWalk1( unif01_Gen* gen, swalk_Res* res )
//******************************************************************************
//  Run Random Walk 1 test
//! \author  J. Bakosi
//******************************************************************************
{
  swalk_RandomWalk1( gen, res, 1, MILLION, 0, 30, 150, 150 );
  return Pvals( { res->H[0]->pVal2[gofw_Mean],
                  res->M[0]->pVal2[gofw_Mean],
                  res->J[0]->pVal2[gofw_Mean],
                  res->R[0]->pVal2[gofw_Mean],
                  res->C[0]->pVal2[gofw_Mean] } );
}
