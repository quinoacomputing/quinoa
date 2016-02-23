//******************************************************************************
/*!
  \file      src/UnitTest/tests/RNG/TestRNG.h
  \author    J. Bakosi
  \date      Tue 23 Feb 2016 03:23:53 PM MST
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Unit tests for RNG/RNG.h
  \details   Unit tests for RNG/RNG.h
*/
//******************************************************************************
#ifndef test_RNG_h
#define test_RNG_h

#include <tut/tut.hpp>

#include "RNG.h"

namespace tut {

//! All tests in group inherited from this base
struct RNG_common {

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
    for (auto n : numbers) ensure( "sample space incorrect", 0.0<n && n<1.0 );
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
    for (auto n : numbers)
      ensure( "sample space incorrect", -eps<n && n<1.0+eps );
    // test first four moments
    test_stats( numbers, a/(a+b), a*b/(a+b)/(a+b)/(a+b+1.0),
                2.0*(b-a)*std::sqrt(a+b+1.0)/(a+b+2.0)/std::sqrt(a*b),
            6.0*((a-b)*(a-b)*(a+b+1.0)-a*b*(a+b+2.0))/a/b/(a+b+2.0)/(a+b+3.0) );
    // test scaled sample space
    for (std::size_t i=0; i<n; ++i)
      r.beta( static_cast<int>(i), num/n, a, b, 0.5, 2.5, &numbers[i*num/n] );
    for (auto n : numbers)
      ensure( "sample space incorrect", 0.5-eps<n && n<3.0+eps );
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
  static void test_move_ctor( rng& r ) {
    std::vector< rng > v;
    v.emplace_back( std::move(r) );
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
  static void test_move_assignment( rng& r ) {
    auto v = std::move(r);
    test_gaussian( v );        // test that the newly moved RNG works
  }

  // Test the first the four first moments of random numbers passed in
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
};

//! Test group shortcuts
using RNG_group = test_group< RNG_common, MAX_TESTS_IN_GROUP >;
using RNG_object = RNG_group::object;

//! Define test group
RNG_group RNG( "RNG/RNG" );

//! Test definitions for group



} // tut::

#endif // test_RNG_h
