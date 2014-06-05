//******************************************************************************
/*!
  \file      src/RNGTest/TestU01SuitePolicy.C
  \author    J. Bakosi
  \date      Wed 21 May 2014 03:33:16 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 suite policy base
  \details   TestU01 suite policy base
*/
//******************************************************************************

extern "C" {
  #include <svaria.h>
  #include <swrite.h>
  #include <bbattery.h>
}

#include <TestU01SuitePolicy.h>

using rngtest::TestU01SuitePolicy;

std::vector< double >
TestU01SuitePolicy::BirthdaySpacings( unif01_Gen* gen, sres_Poisson* res,
  const std::tuple<long, long, int, long, int, int>& xargs )
//******************************************************************************
//  Run Marsaglia's BirthdaySpacings test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  smarsa_BirthdaySpacings( gen, res, get<0>(xargs), get<1>(xargs),
                           get<2>(xargs), get<3>(xargs), get<4>(xargs),
                           get<5>(xargs) );
  return { res->pVal2 };
}

std::vector< double >
TestU01SuitePolicy::Collision( unif01_Gen* gen, sknuth_Res2* res,
  const std::tuple<long, long, int, long, int>& xargs )
//******************************************************************************
//  Run Knuth's Collision test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sknuth_Collision( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                    get<3>(xargs), get<4>(xargs) );
  return { res->Pois->pVal2 };
}

std::vector< double >
TestU01SuitePolicy::Gap( unif01_Gen* gen, sres_Chi2* res,
  const std::tuple<long, long, int, double, double>& xargs )
//******************************************************************************
//  Run Knuth's Gap test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sknuth_Gap( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
              get<3>(xargs), get<4>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::SimpPoker( unif01_Gen* gen, sres_Chi2* res,
  const std::tuple<long, long, int, int, int>& xargs )
//******************************************************************************
//  Run Knuth's Simplified Poker test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sknuth_SimpPoker( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                    get<3>(xargs), get<4>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::CouponCollector( unif01_Gen* gen, sres_Chi2* res,
  const std::tuple<long, long, int, int>& xargs )
//******************************************************************************
//  Run Knuth's Coupon Collector test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sknuth_CouponCollector( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                          get<3>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::MaxOft( unif01_Gen* gen, sknuth_Res1* res,
  const std::tuple<long, long, int, int, int, gofw_TestType,
                   gofw_TestType>& xargs )
//******************************************************************************
//  Run Knuth's Maximum-of-t test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sknuth_MaxOft( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                 get<3>(xargs), get<4>(xargs) );
  return { res->Chi->pVal2[get<5>(xargs)],
           res->Bas->pVal2[get<6>(xargs)] };
}

std::vector< double >
TestU01SuitePolicy::WeightDistrib( unif01_Gen* gen, sres_Chi2* res,
  const std::tuple<long, long, int, long, double, double>& xargs )
//******************************************************************************
//  Run Matsumoto's Weight Distribution test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  svaria_WeightDistrib( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                        get<3>(xargs), get<4>(xargs), get<5>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::MatrixRank( unif01_Gen* gen, sres_Chi2* res,
  const std::tuple<long, long, int, int, int, int>& xargs )
//******************************************************************************
//  Run Marsaglia's Matrix Rank test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  smarsa_MatrixRank( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                     get<3>(xargs), get<4>(xargs), get<5>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::HammingIndep( unif01_Gen* gen, sstring_Res* res,
  const std::tuple<long, long, int, int, int, int>& xargs )
//******************************************************************************
//  Run L'Ecuyer's Hamming Independence test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sstring_HammingIndep( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                        get<3>(xargs), get<4>(xargs), get<5>(xargs) );
  return { res->Bas->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::RandomWalk1( unif01_Gen* gen, swalk_Res* res,
  const std::tuple<long, long, int, int, long, long>& xargs )
//******************************************************************************
//  Run Random Walk 1 test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  swalk_RandomWalk1( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                     get<3>(xargs), get<4>(xargs), get<5>(xargs) );
  return { res->H[0]->pVal2[gofw_Mean],
           res->M[0]->pVal2[gofw_Mean],
           res->J[0]->pVal2[gofw_Mean],
           res->R[0]->pVal2[gofw_Mean],
           res->C[0]->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::SerialOver( unif01_Gen* gen, sres_Basic* res,
  const std::tuple<long, long, int, long, int>& xargs )
//******************************************************************************
//  Run Marsaglia's Serial Over test, t = 2
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  smarsa_SerialOver( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                     get<3>(xargs), get<4>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::CollisionOver( unif01_Gen* gen, smarsa_Res* res,
  const std::tuple<long, long, int, long, int>& xargs )
//******************************************************************************
//  Run Marsaglia's Serial Over test, t = 2
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  smarsa_CollisionOver( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                        get<3>(xargs), get<4>(xargs) );
  return { res->Pois->pVal2 };
}

std::vector< double >
TestU01SuitePolicy::ClosePairs( unif01_Gen* gen, snpair_Res* res,
  const std::tuple<long, long, int, int, int, int, int>& xargs )
//******************************************************************************
//  Run the close-pairs test, t = 2
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  snpair_ClosePairs( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                        get<3>(xargs), get<4>(xargs), get<5>(xargs) );
  if (get<6>(xargs)) {
    return { res->pVal[snpair_NP],
             res->pVal[snpair_mNP],
             res->pVal[snpair_mNP1],
             res->pVal[snpair_mNP2],
             res->pVal[snpair_NJumps],
             res->pVal[snpair_mNP2S] };
  } else {
    return { res->pVal[snpair_NP],
             res->pVal[snpair_mNP],
             res->pVal[snpair_mNP1],
             res->pVal[snpair_mNP2],
             res->pVal[snpair_NJumps] };
  }
}

std::vector< double >
TestU01SuitePolicy::ClosePairsBitMatch( unif01_Gen* gen, snpair_Res* res,
  const std::tuple<long, long, int, int>& xargs )
//******************************************************************************
//  Run the close-pairs test using bit match distance, t = 2
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  snpair_ClosePairsBitMatch( gen, res, get<0>(xargs), get<1>(xargs),
                             get<2>(xargs), get<3>(xargs) );
  return { res->pVal[snpair_BM] };
}

std::vector< double >
TestU01SuitePolicy::Run( unif01_Gen* gen, sres_Chi2* res,
  const std::tuple<long, long, int, int>& xargs )
//******************************************************************************
//  Run Knuth's Run test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sknuth_Run( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
              get<3>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::Permutation( unif01_Gen* gen, sres_Chi2* res,
  const std::tuple<long, long, int, int>& xargs )
//******************************************************************************
//  Run Knuth's Permutation test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sknuth_Permutation( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                      get<3>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::CollisionPermut( unif01_Gen* gen, sknuth_Res2* res,
  const std::tuple<long, long, int, int>& xargs )
//******************************************************************************
//  Run Knuth's Collision test with permutations
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sknuth_CollisionPermut( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                          get<3>(xargs) );
  return { res->Pois->pVal2 };
}

std::vector< double >
TestU01SuitePolicy::SampleProd( unif01_Gen* gen, sres_Basic* res,
  const std::tuple<long, long, int, int>& xargs )
//******************************************************************************
//  Run Sample Products test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  svaria_SampleProd( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                     get<3>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::SampleMean( unif01_Gen* gen, sres_Basic* res,
  const std::tuple<long, long, int>& xargs )
//******************************************************************************
//  Run Sample Mean test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  svaria_SampleMean( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs) );
  return { res->pVal2[gofw_AD] };
}

std::vector< double >
TestU01SuitePolicy::SampleCorr( unif01_Gen* gen, sres_Basic* res,
  const std::tuple<long, long, int, int>& xargs )
//******************************************************************************
//  Run Sample Autocorrelation test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  svaria_SampleCorr( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                     get<3>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::AppearanceSpacings( unif01_Gen* gen, sres_Basic* res,
  const std::tuple<long, long, long, int, int, int>& xargs )
//******************************************************************************
//  Run Maurer's "universal" test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  svaria_AppearanceSpacings( gen, res, get<0>(xargs), get<1>(xargs),
                             get<2>(xargs), get<3>(xargs), get<4>(xargs),
                             get<5>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::SumCollector( unif01_Gen* gen, sres_Chi2* res,
  const std::tuple<long, long, int, double>& xargs )
//******************************************************************************
//  Run Sum Collector test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  svaria_SumCollector( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                       get<3>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::Savir2( unif01_Gen* gen, sres_Chi2* res,
  const std::tuple<long, long, int, long, int>& xargs )
//******************************************************************************
//  Run Marsaglia's modified Savir test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  smarsa_Savir2( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                 get<3>(xargs), get<4>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::GCD( unif01_Gen* gen, smarsa_Res2* res,
  const std::tuple<long, long, int, int>& xargs )
//******************************************************************************
//  Run Marsaglia's greatest common divisor test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  smarsa_GCD( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
              get<3>(xargs) );
  return { res->GCD->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::LinearComp( unif01_Gen* gen, scomp_Res* res,
  const std::tuple<long, long, int, int>& xargs )
//******************************************************************************
//  Run Linear Complexity test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  scomp_LinearComp( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                    get<3>(xargs) );
  return { res->JumpNum->pVal2[gofw_Mean],
           res->JumpSize->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::LempelZiv( unif01_Gen* gen, sres_Basic* res,
  const std::tuple<long, int, int, int>& xargs )
//******************************************************************************
//  Run Lempel-Ziv Compressibility test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  scomp_LempelZiv( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                   get<3>(xargs) );
  return { res->pVal2[gofw_Sum] };
}

std::vector< double >
TestU01SuitePolicy::Fourier3( unif01_Gen* gen, sspectral_Res* res,
  const std::tuple<long, int, int, int>& xargs )
//******************************************************************************
//  Run Fourier3 test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sspectral_Fourier3( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                      get<3>(xargs) );
  return { res->Bas->pVal2[gofw_AD] };
}

std::vector< double >
TestU01SuitePolicy::LongestHeadRun( unif01_Gen* gen, sstring_Res2* res,
  const std::tuple<long, long, int, int, long>& xargs )
//******************************************************************************
//  Run Longest Head Run test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sstring_LongestHeadRun( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                          get<3>(xargs), get<4>(xargs) );
  return { res->Chi->pVal2[gofw_Mean],
           res->Disc->pVal2 };
}

std::vector< double >
TestU01SuitePolicy::PeriodsInStrings( unif01_Gen* gen, sres_Chi2* res,
  const std::tuple<long, long, int, int>& xargs )
//******************************************************************************
//  Run Periods In Strings test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sstring_PeriodsInStrings( gen, res, get<0>(xargs), get<1>(xargs),
                            get<2>(xargs), get<3>(xargs) );
  return { res->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::HammingWeight2( unif01_Gen* gen, sres_Basic* res,
  const std::tuple<long, long, int, int, long>& xargs )
//******************************************************************************
//  Run Hamming Weight 2 test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sstring_HammingWeight2( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                          get<3>(xargs), get<4>(xargs) );
  return { res->pVal2[gofw_Sum] };
}

std::vector< double >
TestU01SuitePolicy::HammingCorr( unif01_Gen* gen, sstring_Res* res,
  const std::tuple<long, long, int, int, int>& xargs )
//******************************************************************************
//  Run Hamming Weight Correlation test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sstring_HammingCorr( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                       get<3>(xargs), get<4>(xargs) );
  return { res->Bas->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::StringRun( unif01_Gen* gen, sstring_Res3* res,
  const std::tuple<long, long, int, int>& xargs )
//******************************************************************************
//  Run String Run test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sstring_Run( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
               get<3>(xargs) );
  return { res->NRuns->pVal2[gofw_Mean],
           res->NBits->pVal2[gofw_Mean] };
}

std::vector< double >
TestU01SuitePolicy::AutoCor( unif01_Gen* gen, sres_Basic* res,
  const std::tuple<long, long, int, int, int>& xargs )
//******************************************************************************
//  Run Autocorrelation test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sstring_AutoCor( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                   get<3>(xargs), get<4>(xargs) );
  return { res->pVal2[gofw_Sum] };
}
