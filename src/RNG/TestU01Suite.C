//******************************************************************************
/*!
  \file      src/RNG/TestU01Suite.C
  \author    J. Bakosi
  \date      Fri 29 Nov 2013 08:16:24 AM MST
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

using rngtest::TestU01Suite;
using Pvals = rngtest::StatTest::Pvals;

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
