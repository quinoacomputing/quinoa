// *****************************************************************************
/*!
  \file      src/UnitTest/tests/RNG/TestRNG.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for RNG/RNG.h
  \details   Unit tests for RNG/RNG.h
*/
// *****************************************************************************
#ifndef test_RNG_h
#define test_RNG_h

#include <functional>

#include <boost/functional/value_factory.hpp>
#include "NoWarning/tut.h"
#include "QuinoaConfig.h"

#include "NoWarning/threefry.h"
#include "NoWarning/philox.h"

#ifdef HAS_MKL
  #include <mkl_vsl_types.h>
  #include "MKLRNG.h"
#endif

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
  #include "RNGSSE.h"
#endif

#include "RNG.h"
#include "Random123.h"

namespace tut {

//! All tests in group inherited from this base
struct RNG_common {

  //! Constructor: create a vector of RNGs to be tested by most tests
  RNG_common() : rngs()
  {
    #ifdef HAS_MKL
    rngs.emplace_back( tk::MKLRNG( 4, VSL_BRNG_MCG31 ) );
    #endif
    #ifdef HAS_RNGSSE2
    rngs.emplace_back( tk::RNGSSE< gm19_state, unsigned, gm19_generate_ >
                                 ( 4, gm19_init_sequence_ ) );
    rngs.emplace_back( tk::RNGSSE< gm29_state, unsigned, gm29_generate_ >
                                 ( 4, gm29_init_short_sequence_ ) );
    rngs.emplace_back( tk::RNGSSE< gm31_state, unsigned, gm31_generate_ >
                                 ( 4, gm31_init_short_sequence_ ) );
    rngs.emplace_back( tk::RNGSSE< gm55_state, unsigned long long,
                                   gm55_generate_ >
                                 ( 4, gm55_init_short_sequence_ ) );
    rngs.emplace_back( tk::RNGSSE< gm61_state, unsigned long long,
                                   gm61_generate_ >
                                 ( 4, gm61_init_sequence_ ) );
    rngs.emplace_back( tk::RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ >
                                 ( 4, gq58x1_init_short_sequence_ ) );
    rngs.emplace_back( tk::RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ >
                                 ( 4, gq58x3_init_short_sequence_ ) );
    rngs.emplace_back( tk::RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ >
                                 ( 4, gq58x4_init_short_sequence_ ) );
    rngs.emplace_back( tk::RNGSSE< mt19937_state, unsigned long long,
                                   mt19937_generate_ >
                                 ( 4, mt19937_init_sequence_ ) );
    rngs.emplace_back( tk::RNGSSE< lfsr113_state, unsigned long long,
                                   lfsr113_generate_ >
                                 ( 4, lfsr113_init_long_sequence_ ) );
    rngs.emplace_back( tk::RNGSSE< mrg32k3a_state, unsigned long long,
                                   mrg32k3a_generate_ >
                                 ( 4, mrg32k3a_init_sequence_ ) );
    #endif
    rngs.emplace_back( tk::Random123< r123::Threefry2x64 >( 4 ) );
    rngs.emplace_back( tk::Random123< r123::Philox2x64 >( 4 ) );
  }

  //! \brief Add a model constructor bound to its arguments to a vector of
  //!   objects modeling tk::RNG
  template< class ModelCtor, typename... Args >
  void add( std::vector< tk::RNG >& v, Args&&... args ) {
    // Bind model constructor to its arguments
    std::function< ModelCtor() > c =
      std::bind( boost::value_factory< ModelCtor >(),
                 std::forward< Args >( args )... );
    // Add to vector
    v.emplace_back( std::move(c) );
  }

  //! Test uniform distribution of a random number generator
  //! \param[in] r RNG to test
  //! \details In real code the member function used to generate random numbers
  //!   are called by different threads, but here we pass a different thread ID
  //!   (first argument to the member function called) will exercise multiple
  //!   streams in serial, i.e., emulate multi-threaded RNG on a single thread.
  //!   This is done here this way because the unit test harness deals out each
  //!   test to a PE, thus calling the generators from real threads would create
  //!   contention.
  template< class rng >
  static void test_uniform( const rng& r ) {
    std::size_t num = 100000;
    std::vector< double > numbers( num );
    auto n = r.nthreads();
    for (std::size_t i=0; i<n; ++i)
      r.uniform( static_cast<int>(i), num/n, &numbers[i*num/n] );
    for (auto m : numbers) ensure( "sample space incorrect", 0.0<m && m<1.0 );
    test_stats( numbers, 0.5, 1.0/12.0, 0.0, -6.0/5.0 );
  }

  //! Test Gaussian distribution of a random number generator
  //! \param[in] r RNG to test
  //! \details In real code the member function used to generate random numbers
  //!   are called by different threads, but here we pass a different thread ID
  //!   (first argument to the member function called) will exercise multiple
  //!   streams in serial, i.e., emulate multi-threaded RNG on a single thread.
  //!   This is done here this way because the unit test harness deals out each
  //!   test to a PE, thus calling the generators from real threads would create
  //!   contention.
  template< class rng >
  static void test_gaussian( const rng& r ) {
    std::size_t num = 100000;
    std::vector< double > numbers( num );
    auto n = r.nthreads();
    for (std::size_t i=0; i<n; ++i)
      r.gaussian( static_cast<int>(i), num/n, &numbers[i*num/n] );
    test_stats( numbers, 0.0, 1.0, 0.0, 0.0 );
  }

  //! Test beta distribution of a random number generator
  //! \param[in] r RNG to test
  //! \details In real code the member function used to generate random numbers
  //!   are called by different threads, but here we pass a different thread ID
  //!   (first argument to the member function called) will exercise multiple
  //!   streams in serial, i.e., emulate multi-threaded RNG on a single thread.
  //!   This is done here this way because the unit test harness deals out each
  //!   test to a PE, thus calling the generators from real threads would create
  //!   contention.
  template< class rng >
  static void test_beta( const rng& r ) {
    std::size_t num = 100000;
    double a = 0.2;
    double b = 0.6;
    std::vector< double > numbers( num );
    auto n = r.nthreads();
    for (std::size_t i=0; i<n; ++i)
      r.beta( static_cast<int>(i), num/n, a, b, 0.0, 1.0, &numbers[i*num/n] );
    // test sample space
    for (auto m : numbers)
      ensure( "sample space incorrect", -eps<m && m<1.0+eps );
    // test first four moments
    test_stats( numbers, a/(a+b), a*b/(a+b)/(a+b)/(a+b+1.0),
                2.0*(b-a)*std::sqrt(a+b+1.0)/(a+b+2.0)/std::sqrt(a*b),
            6.0*((a-b)*(a-b)*(a+b+1.0)-a*b*(a+b+2.0))/a/b/(a+b+2.0)/(a+b+3.0) );
    // test scaled sample space
    for (std::size_t i=0; i<n; ++i)
      r.beta( static_cast<int>(i), num/n, a, b, 0.5, 2.5, &numbers[i*num/n] );
    for (auto m : numbers)
      ensure( "sample space incorrect", 0.5-eps<m && m<3.0+eps );
  }

  //! Test copy constructor of a random number generator
  //! \param[in] r RNG to test
  template< class rng >
  static void test_copy_ctor( const rng& r ) {
    std::vector< rng > v;
    v.push_back( r );
    test_gaussian( r );         // test that source of copy still works
    test_gaussian( v[0] );      // test that the copy works
  }

  //! Test move constructor of a random number generator
  //! \param[in] r RNG to test
  template< class rng >
  static void test_move_ctor( const rng& r ) {
    std::vector< rng > v;
    auto p = r;
    v.emplace_back( std::move(p) );
    test_gaussian( v[0] );      // test that the newly moved RNG works
  }

  //! Test copy assignment of a random number generator
  //! \param[in] r RNG to test
  template< class rng >
  static void test_copy_assignment( const rng& r ) {
    auto v = r;
    test_gaussian( r );         // test that source of copy still works
    test_gaussian( v );         // test that the copy works
  }

  //! Test move assignment of a random number generator
  //! \param[in] r RNG to test
  template< class rng >
  static void test_move_assignment( const rng& r ) {
    auto p = r;
    auto v = std::move(p);
    test_gaussian( v );        // test that the newly moved RNG works
  }

  // Test the first four moments of random numbers passed in
  static void test_stats( const std::vector< double >& numbers,
                          double correct_mean,
                          double correct_variance,
                          double correct_skewness,
                          double correct_excess_kurtosis )
  {
    const double precision = 0.05;
    std::array< double, 4 > s{{ 0.0, 0.0, 0.0, 0.0 }};
    for (auto n : numbers) s[0] += n;
    s[0] /= static_cast< double >( numbers.size() );
    ensure_equals( "mean inaccurate", s[0], correct_mean, precision );
    for (auto n : numbers) {
      s[1] += (n-s[0])*(n-s[0]);
      s[2] += (n-s[0])*(n-s[0])*(n-s[0]);
      s[3] += (n-s[0])*(n-s[0])*(n-s[0])*(n-s[0]);
    }
    s[1] /= static_cast< double >( numbers.size() );
    s[2] /= static_cast< double >( numbers.size() );
    s[3] /= static_cast< double >( numbers.size() );
    ensure_equals( "variance inaccurate", s[1], correct_variance, precision );
    ensure_equals( "skewness inaccurate",
                   s[2]/std::pow(s[1],1.5), correct_skewness, precision );
    ensure_equals( "excess kurtosis inaccurate",
                   s[3]/s[1]/s[1]-3.0, correct_excess_kurtosis, precision );
  }

  static constexpr double eps = 1.0e-9;
  std::vector< tk::RNG > rngs;
};

//! Test group shortcuts
using RNG_group = test_group< RNG_common, MAX_TESTS_IN_GROUP >;
using RNG_object = RNG_group::object;

//! Define test group
static RNG_group RNG( "RNG/RNG" );

//! Test definitions for group

//! Test constructor taking an object modeling Concept in tk::RNG
//! \author J. Bakosi
template<> template<>
void RNG_object::test< 1 >() {
  set_test_name( "ctor( rng() ) & nthreads()" );
  for (const auto& r : rngs)
    ensure_equals( "nthreads() via polymorphic tk::RNG call incorrect",
                   r.nthreads(), 4 );
}

//! \brief Test constructor taking a function pointer to a constructor of an
//!   object modeling Concept in tk::RNG
//! \author J. Bakosi
template<> template<>
void RNG_object::test< 2 >() {
  set_test_name( "ctor( std::function<rng>(rng) )" );

  std::vector< tk::RNG > v;

  #ifdef HAS_MKL
  add< tk::MKLRNG >( v, 4, VSL_BRNG_MCG31 );
  #endif

  #ifdef HAS_RNGSSE2
  add< tk::RNGSSE< gm19_state, unsigned, gm19_generate_ > >
     ( v, 4, gm19_init_sequence_ );
  add< tk::RNGSSE< gm29_state, unsigned, gm29_generate_ > >
     ( v, 4, gm29_init_short_sequence_ );
  add< tk::RNGSSE< gm31_state, unsigned, gm31_generate_ > >
     ( v, 4, gm31_init_short_sequence_ );
  add< tk::RNGSSE< gm55_state, unsigned long long, gm55_generate_ > >
     ( v, 4, gm55_init_short_sequence_ );
  add< tk::RNGSSE< gm61_state, unsigned long long, gm61_generate_ > >
     ( v, 4, gm61_init_sequence_ );
  add< tk::RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ > >
     ( v, 4, gq58x1_init_short_sequence_ );
  add< tk::RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ > >
     ( v, 4, gq58x3_init_short_sequence_ );
  add< tk::RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ > >
     ( v, 4, gq58x4_init_short_sequence_ );
  add< tk::RNGSSE< mt19937_state, unsigned long long, mt19937_generate_ > >
     ( v, 4, mt19937_init_sequence_ );
  add< tk::RNGSSE< lfsr113_state, unsigned long long, lfsr113_generate_ > >
     ( v, 4, lfsr113_init_long_sequence_ );
  add< tk::RNGSSE< mrg32k3a_state, unsigned long long, mrg32k3a_generate_ > >
     ( v, 4, mrg32k3a_init_sequence_ );
  #endif

  add< tk::Random123< r123::Threefry2x64 > >( v, 4 );
  add< tk::Random123< r123::Philox2x64 > >( v, 4 );

  for (const auto& r : v)
    ensure_equals( "nthreads() via polymorphic tk::RNG call incorrect",
                   r.nthreads(), 4 );
}

//! Test Gaussian generator statistics via polymorphic call in tk::RNG
//! \author J. Bakosi
template<> template<>
void RNG_object::test< 3 >() {
  set_test_name( "Gaussian from 4 emulated streams" );
  for (const auto& r : rngs) test_gaussian( r );
}

//! Test beta generator statistics via polymorphic call in tk::RNG
//! \author J. Bakosi
template<> template<>
void RNG_object::test< 4 >() {
  set_test_name( "beta from 4 emulated streams" );
  for (const auto& r : rngs) test_beta( r );
}

//! Test uniform generator statistics via polymorphic call in tk::RNG
//! \author J. Bakosi
template<> template<>
void RNG_object::test< 5 >() {
  set_test_name( "uniform from 4 emulated streams" );
  for (const auto& r : rngs) test_uniform( r );
}

//! Test copy constructor
//! \author J. Bakosi
template<> template<>
void RNG_object::test< 6 >() {
  set_test_name( "copy constructor" );
  for (const auto& r : rngs) test_copy_ctor( r );
}

//! Test move constructor
//! \author J. Bakosi
template<> template<>
void RNG_object::test< 7 >() {
  set_test_name( "move constructor" );
  for (const auto& r : rngs) test_move_ctor( r );
}

//! Test copy assignment
//! \author J. Bakosi
template<> template<>
void RNG_object::test< 8 >() {
  set_test_name( "copy assignment" );
  for (const auto& r : rngs) test_copy_assignment( r );
}

//! Test move assignment
//! \author J. Bakosi
template<> template<>
void RNG_object::test< 9 >() {
  set_test_name( "move assignment" );
  for (const auto& r : rngs) test_move_assignment( r );
}

} // tut::

#endif // test_RNG_h
