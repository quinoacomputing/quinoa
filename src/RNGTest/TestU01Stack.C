//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Stack.C
  \author    J. Bakosi
  \date      Sun 29 Jun 2014 08:28:40 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Stack of TestU01 tests
  \details   Stack of TestU01 tests
*/
//******************************************************************************

extern "C" {
  #include <svaria.h>
  #include <swrite.h>
  #include <bbattery.h>
}

#include <TestU01Stack.h>
#include <RNG.h>
#include <RNGTest/InputDeck/InputDeck.h>

namespace rngtest {

extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;
extern ctr::InputDeck g_inputdeck;

//! Global thread id
static int g_tid = 0;
#ifdef _OPENMP
#pragma omp threadprivate(g_tid)
#endif

template< tk::ctr::RawRNGType id >
static double uniform( void*, void* )
//******************************************************************************
//  Global-scope TestU01 uniform RNG wrapper
//! \author  J. Bakosi
//******************************************************************************
{
  double r = 0.0;
  const auto rng = g_rng.find( id );
  if (rng != end(g_rng))
    rng->second.uniform( g_tid, 1, &r );
  else
    Throw( tk::ExceptType::FATAL, "RNG not found" );
  return r;
}

template< tk::ctr::RawRNGType id >
static unsigned long uniform_bits( void*, void* )
//******************************************************************************
//  Global-scope TestU01 uniform RNG bits wrapper
//! \author  J. Bakosi
//******************************************************************************
{
  double r = 0.0;
  const auto rng = g_rng.find( id );
  if (rng != end(g_rng))
    rng->second.uniform( g_tid, 1, &r );
  else
    Throw( tk::ExceptType::FATAL, "RNG not found" );
  return static_cast<unsigned long>(r * unif01_NORM32);
}

template< tk::ctr::RawRNGType id >
inline unif01_Gen* createTestU01Gen( const std::string& name )
//******************************************************************************
//  Global-scope create TestU01 external generator
//! \author  J. Bakosi
//******************************************************************************
{
  return unif01_CreateExternGen01( const_cast<char*>(name.c_str()),
                                   uniform< id >, uniform_bits< id > );
}

} // rngtest::

using rngtest::TestU01Stack;

TestU01Stack::TestU01Stack()
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Associate RNGs to global-scope wrappers. Admittedly, the code below is ugly
  // and looks stupid at first sight. However, this is a translation of runtime
  // information (user-selected RNGs) to compile-time information: associating
  // RNG ids from an enum class, tk::ctr::RNGType::value, to a compile-time
  // constant, underlying_type value, facilitating a different pair of
  // global-scope RNG wrappers (uniform and uniform_bits) with code reuse. Note
  // that uniform and uniform_bits wrappers must be global-scope as they are
  // used as external generators to TestU01. Templating them on the id enables
  // the compiler generate a different wrapper for a different RNG facilitating
  // simultaneous calls to any or all wrappers as they are unique functions.
  for (const auto& r : g_inputdeck.get< tag::selected, tk::tag::rng >()) {
    using tk::ctr::RNGType;
    using tk::ctr::raw;
    #ifdef HAS_MKL
    if (r == RNGType::MKL_MCG31)
      addRNG< raw(RNGType::MKL_MCG31) >( r );
    else if (r == RNGType::MKL_R250)
      addRNG< raw(RNGType::MKL_R250) >( r );
    else if (r == RNGType::MKL_MRG32K3A)
      addRNG< raw(RNGType::MKL_MRG32K3A) >( r );
    else if (r == RNGType::MKL_MCG59)
      addRNG< raw(RNGType::MKL_MCG59) >( r );
    else if (r == RNGType::MKL_WH)
      addRNG< raw(RNGType::MKL_WH) >( r );
    else if (r == RNGType::MKL_MT19937)
      addRNG< raw(RNGType::MKL_MT19937) >( r );
    else if (r == RNGType::MKL_MT2203)
      addRNG< raw(RNGType::MKL_MT2203) >( r );
    else if (r == RNGType::MKL_SFMT19937)
      addRNG< raw(RNGType::MKL_SFMT19937) >( r );
    else if (r == RNGType::MKL_SOBOL)
      addRNG< raw(RNGType::MKL_SOBOL) >( r );
    else if (r == RNGType::MKL_NIEDERR)
      addRNG< raw(RNGType::MKL_NIEDERR) >( r );
    //else if (r == RNGType::MKL_IABSTRACT)
    //  addRNG< raw(RNGType::MKL_IABSTRACT) >( r );
    //else if (r == RNGType::MKL_DABSTRACT)
    //  addRNG< raw(RNGType::MKL_DABSTRACT) >( r );
    //else if (r == RNGType::MKL_SABSTRACT)
    //  addRNG< raw(RNGType::MKL_SABSTRACT) >( r );
    else if (r == RNGType::MKL_NONDETERM)
      addRNG< raw(RNGType::MKL_NONDETERM) >( r );
    else
    #endif
    if (r == RNGType::RNGSSE_GM19)         
      addRNG< raw(RNGType::RNGSSE_GM19) >( r );
    else if (r == RNGType::RNGSSE_GM29)    
      addRNG< raw(RNGType::RNGSSE_GM29) >( r );
    else if (r == RNGType::RNGSSE_GM31)    
      addRNG< raw(RNGType::RNGSSE_GM31) >( r );
    else if (r == RNGType::RNGSSE_GM55)    
      addRNG< raw(RNGType::RNGSSE_GM55) >( r );
    else if (r == RNGType::RNGSSE_GM61)    
      addRNG< raw(RNGType::RNGSSE_GM61) >( r );
    else if (r == RNGType::RNGSSE_GQ581)   
      addRNG< raw(RNGType::RNGSSE_GQ581) >( r );
    else if (r == RNGType::RNGSSE_GQ583)   
      addRNG< raw(RNGType::RNGSSE_GQ583) >( r );
    else if (r == RNGType::RNGSSE_GQ584)   
      addRNG< raw(RNGType::RNGSSE_GQ584) >( r );
    else if (r == RNGType::RNGSSE_MT19937) 
      addRNG< raw(RNGType::RNGSSE_MT19937) >( r );
    else if (r == RNGType::RNGSSE_LFSR113) 
      addRNG< raw(RNGType::RNGSSE_LFSR113) >( r );
    else if (r == RNGType::RNGSSE_MRG32K3A)
      addRNG< raw(RNGType::RNGSSE_MRG32K3A) >( r );
  }
}

template< tk::ctr::RawRNGType id >
void TestU01Stack::addRNG( tk::ctr::RNGType r )
//******************************************************************************
//! Create TestU01 RNG wrapper
//! \author  J. Bakosi
//******************************************************************************
{
  m_rngprops.emplace( r, createTestU01Gen< id >( tk::ctr::RNG().name(r) ) );
}

const TestU01Stack::RNGProps&
TestU01Stack::rngprops( tk::ctr::RNGType r ) const
//******************************************************************************
//! Find RNG properties based on RNG id
//! \author  J. Bakosi
//******************************************************************************
{
  const auto prop = m_rngprops.find( r );
  if (prop == end(m_rngprops)) Throw( tk::ExceptType::FATAL, "RNG not found" );
  return prop->second;
}

std::vector< double >
TestU01Stack::BirthdaySpacings( unif01_Gen* gen, sres_Poisson* res,
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
TestU01Stack::Collision( unif01_Gen* gen, sknuth_Res2* res,
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
TestU01Stack::Gap( unif01_Gen* gen, sres_Chi2* res,
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
TestU01Stack::SimpPoker( unif01_Gen* gen, sres_Chi2* res,
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
TestU01Stack::CouponCollector( unif01_Gen* gen, sres_Chi2* res,
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
TestU01Stack::MaxOft( unif01_Gen* gen, sknuth_Res1* res,
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
TestU01Stack::WeightDistrib( unif01_Gen* gen, sres_Chi2* res,
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
TestU01Stack::MatrixRank( unif01_Gen* gen, sres_Chi2* res,
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
TestU01Stack::HammingIndep( unif01_Gen* gen, sstring_Res* res,
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
TestU01Stack::RandomWalk1( unif01_Gen* gen, swalk_Res* res,
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
TestU01Stack::SerialOver( unif01_Gen* gen, sres_Basic* res,
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
TestU01Stack::CollisionOver( unif01_Gen* gen, smarsa_Res* res,
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
TestU01Stack::ClosePairs( unif01_Gen* gen, snpair_Res* res,
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
TestU01Stack::ClosePairsBitMatch( unif01_Gen* gen, snpair_Res* res,
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
TestU01Stack::Run( unif01_Gen* gen, sres_Chi2* res,
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
TestU01Stack::Permutation( unif01_Gen* gen, sres_Chi2* res,
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
TestU01Stack::CollisionPermut( unif01_Gen* gen, sknuth_Res2* res,
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
TestU01Stack::SampleProd( unif01_Gen* gen, sres_Basic* res,
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
TestU01Stack::SampleMean( unif01_Gen* gen, sres_Basic* res,
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
TestU01Stack::SampleCorr( unif01_Gen* gen, sres_Basic* res,
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
TestU01Stack::AppearanceSpacings( unif01_Gen* gen, sres_Basic* res,
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
TestU01Stack::SumCollector( unif01_Gen* gen, sres_Chi2* res,
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
TestU01Stack::Savir2( unif01_Gen* gen, sres_Chi2* res,
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
TestU01Stack::GCD( unif01_Gen* gen, smarsa_Res2* res,
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
TestU01Stack::LinearComp( unif01_Gen* gen, scomp_Res* res,
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
TestU01Stack::LempelZiv( unif01_Gen* gen, sres_Basic* res,
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
TestU01Stack::Fourier3( unif01_Gen* gen, sspectral_Res* res,
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
TestU01Stack::LongestHeadRun( unif01_Gen* gen, sstring_Res2* res,
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
TestU01Stack::PeriodsInStrings( unif01_Gen* gen, sres_Chi2* res,
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
TestU01Stack::HammingWeight2( unif01_Gen* gen, sres_Basic* res,
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
TestU01Stack::HammingCorr( unif01_Gen* gen, sstring_Res* res,
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
TestU01Stack::StringRun( unif01_Gen* gen, sstring_Res3* res,
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
TestU01Stack::AutoCor( unif01_Gen* gen, sres_Basic* res,
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
