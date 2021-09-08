// *****************************************************************************
/*!
  \file      tests/unit/RNG/TestRNG.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Reused test code for unit testing random number generators
  \details   Reused test code for unit testing random number generators
*/
// *****************************************************************************

#include <functional>

#include "NoWarning/value_factory.hpp"

#include "QuinoaConfig.hpp"

#include "NoWarning/threefry.hpp"
#include "NoWarning/philox.hpp"

#ifdef HAS_MKL
  #include <mkl_vsl_types.h>
  #include "MKLRNG.hpp"
#endif

#ifdef HAS_MKL
  #include <mkl_lapacke.h>
#else
  #include <lapacke.h>
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
  #include "RNGSSE.hpp"
#endif

#include "RNG.hpp"
#include "Random123.hpp"

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

  // Compute the Cholesky decomposition of a covariance matrix
  //! \tparam D Number of dimensions of covariance matrix
  //! \param[in] C Covariance matrix to decompose, only the upper diagonal and
  //!   the diagonal are stored as the matrix is assumed to be symmetric (and
  //!   positive definite)
  //! \return Choleksy decomposition (upper triangle only)
  template< std::size_t D >
  static std::array< double, D*(D+1)/2 >
  cholesky( const std::array< double, D*(D+1)/2 >& C ) {
    auto cov = C;
    lapack_int ndim = static_cast< lapack_int >( D );
    #ifndef NDEBUG
    lapack_int info =
    #endif
      LAPACKE_dpptrf( LAPACK_ROW_MAJOR, 'U', ndim, cov.data() );
    Assert( info == 0, "Error in Cholesky-decomposition" );
    return cov;
  }

  //! Test multi-variate Gaussian distribution of a random number generator
  //! \tparam D Number of dimensions of covariance matrix
  //! \tparam rng Random number generator class type whose instance to test
  //! \param[in] r RNG to test
  //! \param[in] mean Vector of means
  //! \param[in] C Upper triangle of the covariance matrix
  //! \details In real code the member function used to generate random numbers
  //!   are called by different threads, but here we pass a different thread ID
  //!   (first argument to the member function called) will exercise multiple
  //!   streams in serial, i.e., emulate multi-threaded RNG on a single thread.
  //!   This is done here this way because the unit test harness deals out each
  //!   test to a PE, thus calling the generators from real threads would create
  //!   contention.
  template< std::size_t D, class rng >
  static void test_gaussianmv( const rng& r,
                               const std::array< double, D >& mean,
                               const std::array< double, D*(D+1)/2 >& C )
  {
    const auto n = r.nthreads();
    const std::size_t num = 1000000;
    // Compute Cholesky decomposition of covariance matrix
    auto cov = cholesky< D >( C );
    // Generate joint Gaussian random vectors of D dimension
    std::vector< double > numbers( num*D );
    for (std::size_t i=0; i<n; ++i)
      r.gaussianmv( static_cast<int>(i), num/n, D, mean.data(), cov.data(),
                    &numbers[i*num*D/n] );
    // Test first two moments
    test_statsmv< D >( numbers, mean, C );
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
    rng v( r );
    v = r;
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
  //! \param[in] numbers Random numbers to test
  //! \param[in] correct_mean Baseline mean to compare to
  //! \param[in] correct_variance Baseline variance to compare to
  //! \param[in] correct_skewness Baseline skewness to compare to
  //! \param[in] correct_excess_kurtosis Baseline excess kurtosis to compare to
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

  // Test the first the two moments of multi-variate random numbers passed in
  //! \tparam D Number of dimensions of sample space in multi-variate
  //!    distribution
  //! \param[in] numbers Random vectors to test, size: num*D, vector size: D
  //! \param[in] correct_mean Baseline mean vector to compare to
  //! \param[in] correct_cov Baseline covariance matrix to compare to
  template< std::size_t D >
  static void test_statsmv(
     const std::vector< double >& numbers,
     const std::array< double, D >& correct_mean,
     const std::array< double, D*(D+1)/2 >& correct_cov )
  {
    const double precision = 0.05;
    const auto N = numbers.size() / D;

    // Compute and test means
    std::array< double, D > mean;
    for (std::size_t i=0; i<D; ++i) {
      mean[i] = 0.0;
      for (std::size_t j=0; j<N; ++j) mean[i] += numbers[j*D+i];
      mean[i] /= static_cast<tk::real>(N);
      ensure_equals( "mean inaccurate", mean[i], correct_mean[i],
                     std::abs(correct_mean[i])*precision );
    }

    // Compute and test covariance matrix
    std::array< double, D*(D+1)/2 > cov;
    cov.fill( 0.0 );
    std::size_t e = 0;
    for (std::size_t r=0; r<D; ++r)
      for (std::size_t c=0; c<D; ++c) {
        if (r<=c) {
          for (std::size_t j=0; j<N; ++j)
            cov[e] += (numbers[j*D+r] - mean[r]) * (numbers[j*D+c] - mean[c]);
          ++e;
        }
      }
    for (std::size_t c=0; c<cov.size(); ++c) {
      cov[c] /= static_cast<tk::real>(N);
      ensure_equals( "covariance matrix entry " + std::to_string(c) +
                     " inaccurate", cov[c], correct_cov[c],
                     std::abs(correct_cov[c])*precision );
    }
  }

  static constexpr double eps = 1.0e-9;
  std::vector< tk::RNG > rngs;
};

} // tut::
