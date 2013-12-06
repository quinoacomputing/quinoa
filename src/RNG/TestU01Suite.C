//******************************************************************************
/*!
  \file      src/RNG/TestU01Suite.C
  \author    J. Bakosi
  \date      Fri 06 Dec 2013 12:56:55 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 suite
  \details   TestU01 suite
*/
//******************************************************************************

extern "C" {
  #include <smarsa.h>
  #include <svaria.h>
  #include <swrite.h>
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

TestU01Suite::TestU01Suite( const Base& base, const std::string& name )
  : Battery( base ), m_name( name )
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
  //   * The wrappers, uniform() and uniform_bits(), this way can be templates,
  //     which facilitates code-reuse
  //   * The wrappers must be in global scope as they are passed to TestU01
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

void
TestU01Suite::total()
//******************************************************************************
//  Count up total number of statistics expected
//! \author  J. Bakosi
//******************************************************************************
{
  m_npval = 0;
  for (const auto& t : m_pvals) {
    m_npval += t.size();
  }
}


rngtest::StatTest::Pvals::size_type
TestU01Suite::failed()
//******************************************************************************
//  Count up number of failed tests
//! \author  J. Bakosi
//******************************************************************************
{
  using Pval = StatTest::Pvals::value_type;

  StatTest::Pvals::size_type failed = 0;
  for (const auto& t : m_pvals) {
    for (const auto& p : t) {
      if (fabs(p+1.0) > std::numeric_limits< Pval >::epsilon()) {
        ++failed;
      }
    }
  }
  return failed;
}

void
TestU01Suite::run()
//******************************************************************************
//  Run battery
//! \author  J. Bakosi
//******************************************************************************
{
  const RNGTestPrint& print = m_base.print;

  print.part( m_name );
  print.section( "Statistical tests finished" );

  swrite_Basic = FALSE;         // Want screen no putput from TestU01

//   g_tid = 0;
//   bbattery_SmallCrush( m_rng[2].get() );

  using Psize = StatTest::Pvals::size_type;
  using Tsize = TestContainer::size_type;

  Tsize ntest = m_tests.size();
  Psize n = 0;
  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef _OPENMP
    g_tid = omp_get_thread_num();
    #else
    g_tid = 0;
    #endif

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (Tsize i=0; i<ntest; ++i) {

      // Run test
      StatTest::Pvals pvals = m_tests[i]->run();

      // Evaluate test
      Psize npval = pvals.size();
      for (Psize p=0; p<npval; ++p) {

        // Increase number tests completed
        #ifdef _OPENMP
        #pragma omp atomic
        #endif
        ++n;

        // Output one-liner
        print.test< StatTest, TestContainer >
                  ( n, m_npval, m_tests[i], p, pvals[p] );

        // Save suspect p-value
        if ((pvals[p] <= gofw_Suspectp) || (pvals[p] >= 1.0 - gofw_Suspectp)) {
          m_pvals[i][p] = pvals[p];
        }
      }
    }
  }

  // Count up number of filed tests
  StatTest::Pvals::size_type nfail = failed();

  // Output summary of failed tests
  if (nfail) {
    print.failed< StatTest >("Failed tests", m_npval, nfail, m_pvals, m_tests);
  } else {
    print.note("All tests passed");
  }

  print.endpart();
}

void
TestU01Suite::print() const
//******************************************************************************
//  Print list of registered statistical tests
//! \author  J. Bakosi
//******************************************************************************
{
  // Output test names (only for the first RNG tested, the rest are repeats)
  m_base.print.names< StatTest >( m_tests, m_tests.size()/m_testRNGs.size() );
}

Pvals
TestU01Suite::BirthdaySpacings( unif01_Gen* gen, sres_Poisson* res,
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
  return Pvals( { res->pVal2 } );
}

Pvals
TestU01Suite::Collision( unif01_Gen* gen, sknuth_Res2* res,
                         const std::tuple<long, long, int, long, int>& xargs )
//******************************************************************************
//  Run Knuth's Collision test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sknuth_Collision( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                    get<3>(xargs), get<4>(xargs) );
  return Pvals( { res->Pois->pVal2 } );
}

Pvals
TestU01Suite::Gap( unif01_Gen* gen, sres_Chi2* res,
                   const std::tuple<long, long, int, double, double>& xargs )
//******************************************************************************
//  Run Knuth's Gap test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sknuth_Gap( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
              get<3>(xargs), get<4>(xargs) );
  return Pvals( { res->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::SimpPoker( unif01_Gen* gen, sres_Chi2* res,
                         const std::tuple<long, long, int, int, int>& xargs )
//******************************************************************************
//  Run Knuth's Simplified Poker test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sknuth_SimpPoker( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                    get<3>(xargs), get<4>(xargs) );
  return Pvals( { res->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::CouponCollector( unif01_Gen* gen, sres_Chi2* res,
                               const std::tuple<long, long, int, int>& xargs )
//******************************************************************************
//  Run Knuth's Coupon Collector test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sknuth_CouponCollector( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                          get<2>(xargs) );
  return Pvals( { res->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::MaxOft( unif01_Gen* gen, sknuth_Res1* res,
                      const std::tuple<long, long, int, int, int>& xargs )
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
TestU01Suite::WeightDistrib( unif01_Gen* gen, sres_Chi2* res,
                const std::tuple<long, long, int, long, double, double>& xargs )
//******************************************************************************
//  Run Matsumoto's Weight Distribution test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  svaria_WeightDistrib( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                        get<3>(xargs), get<4>(xargs), get<5>(xargs) );
  return Pvals( { res->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::MatrixRank( unif01_Gen* gen, sres_Chi2* res,
                const std::tuple<long, long, int, int, int, int>& xargs )
//******************************************************************************
//  Run Marsaglia's Matrix Rank test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  smarsa_MatrixRank( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                     get<3>(xargs), get<4>(xargs), get<5>(xargs) );
  return Pvals( { res->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::HammingIndep( unif01_Gen* gen, sstring_Res* res,
                const std::tuple<long, long, int, int, int, int>& xargs )
//******************************************************************************
//  Run L'Ecuyer's Hamming Independence test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  sstring_HammingIndep( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                        get<3>(xargs), get<4>(xargs), get<5>(xargs) );
  return Pvals( { res->Bas->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::RandomWalk1( unif01_Gen* gen, swalk_Res* res,
                const std::tuple<long, long, int, int, long, long>& xargs )
//******************************************************************************
//  Run Random Walk 1 test
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  swalk_RandomWalk1( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                     get<3>(xargs), get<4>(xargs), get<5>(xargs) );
  return Pvals( { res->H[0]->pVal2[gofw_Mean],
                  res->M[0]->pVal2[gofw_Mean],
                  res->J[0]->pVal2[gofw_Mean],
                  res->R[0]->pVal2[gofw_Mean],
                  res->C[0]->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::SerialOver( unif01_Gen* gen, sres_Basic* res,
                          const std::tuple<long, long, int, long, int>& xargs )
//******************************************************************************
//  Run Marsaglia's Serial Over test, t = 2
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  smarsa_SerialOver( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                     get<3>(xargs), get<4>(xargs) );
  return Pvals( { res->pVal2[gofw_Mean] } );
}

Pvals
TestU01Suite::CollisionOver( unif01_Gen* gen, smarsa_Res* res,
                const std::tuple<long, long, int, long, int>& xargs )
//******************************************************************************
//  Run Marsaglia's Serial Over test, t = 2
//! \author  J. Bakosi
//******************************************************************************
{
  using std::get;
  smarsa_CollisionOver( gen, res, get<0>(xargs), get<1>(xargs), get<2>(xargs),
                        get<3>(xargs), get<4>(xargs) );
  return Pvals( { res->Pois->pVal2 } );
}
