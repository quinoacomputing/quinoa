// *****************************************************************************
/*!
  \file      src/RNG/Random123.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Interface to Random123 random number generators
  \details   Interface to Random123 random number generators
*/
// *****************************************************************************
#ifndef Random123_h
#define Random123_h

#include <cstring>
#include <random>
#include <limits>
#include <array>

#include "NoWarning/uniform.h"
#include "NoWarning/beta_distribution.h"

#include "Make_unique.h"
#include "Exception.h"
#include "Keywords.h"
#include "Macro.h"

namespace tk {

//! Random123-based random number generator used polymorphically with tk::RNG
template< class CBRNG >
class Random123 {

  private:
    static const std::size_t CBRNG_DATA_SIZE = 3;
    using ncomp_t = kw::ncomp::info::expect::type;    
    using ctr_type = typename CBRNG::ctr_type;
    using key_type = typename CBRNG::key_type;
    using value_type = typename CBRNG::ctr_type::value_type;
    using arg_type = std::vector< std::array< value_type, CBRNG_DATA_SIZE > >;

    //! Adaptor to use a std distribution with the Random123 generator
    //! \see C++ concepts: UniformRandomNumberGenerator
    struct Adaptor {
      using result_type = unsigned long;
      Adaptor( CBRNG& r, arg_type& d, int t ) : rng(r), data(d), tid(t) {}
      static constexpr result_type min() { return 0u; }
      static constexpr result_type max() {
        return std::numeric_limits< result_type >::max();
      }
      result_type operator()()
      {
        auto& d = data[ static_cast< std::size_t >( tid ) ];
        d[2] = static_cast< result_type >( tid );
        ctr_type ctr = {{ d[0], d[1] }};      // assemble counter
        key_type key = {{ d[2] }};            // assemble key
        auto res = rng( ctr, key );           // generate
        ctr.incr();
        d[0] = ctr[0];
        d[1] = ctr[1];
        return res[0];
      }
      CBRNG& rng;
      arg_type& data;
      int tid;
    };

  public:
    //! Constructor
    //! \param[in] n Initialize RNG using this many independent streams
    //! \param[in] seed RNG seed
    explicit Random123( uint64_t n = 1, uint32_t seed = 0 ) {
      Assert( n > 0, "Need at least one thread" );
      m_data.resize( n, {{ 0, static_cast< uint64_t >( seed ) << 32, 0 }} );
    }

    //! Uniform RNG: Generate uniform random numbers
    //! \param[in] tid Thread (or more precisely) stream ID
    //! \param[in] num Number of RNGs to generate
    //! \param[in,out] r Pointer to memory to write the random numbers to
    void uniform( int tid, ncomp_t num, double* r ) const {
      auto& d = m_data[ static_cast< std::size_t >( tid ) ];
      for (ncomp_t i=0; i<num; ++i) {
        d[2] = static_cast< unsigned long >( tid );
        ctr_type ctr = {{ d[0], d[1] }};      // assemble counter
        key_type key = {{ d[2] }};            // assemble key
        auto res = m_rng( ctr, key );         // generate
        r[i] = r123::u01fixedpt< double, value_type >( res[0] );
        ctr.incr();
        d[0] = ctr[0];
        d[1] = ctr[1];
      }
    }

    //! Gaussian RNG: Generate Gaussian random numbers
    //! \param[in] tid Thread (or more precisely stream) ID
    //! \param[in] num Number of RNGs to generate
    //! \param[in,out] r Pointer to memory to write the random numbers to
    //! \details Generating Gaussian random numbers is implemented via an
    //!   adaptor, modeling std::UniformRandomNumberGenerator, outsourcing the
    //!   transformation of uniform random numbers to Gaussian ones, to the
    //!   standard library. The adaptor is instantiated here because a standard
    //!   distribution, such as e.g., std::normal_distribution, generates
    //!   numbers using operator() with no arguments, thus the RNG state and the
    //!   thread ID (this latter only known here) must be stored in the adaptor
    //!   functor's state. Even though creating the adaptor seems like a
    //!   potentially costly operation for every call, using the standard
    //!   library implementation is still faster than a hand-coded
    //!   implementation of the Box-Muller algorithm. Note that libc++ uses a
    //!   cache, as Box-Muller, implemented using the polar algorithm generates
    //!   2 Gaussian numbers for each pair of uniform ones, caching every 2nd.
    void gaussian( int tid, ncomp_t num, double* r ) const {
      Adaptor generator( m_rng, m_data, tid );
      std::normal_distribution<> gauss_dist( 0.0, 1.0 );
      for (ncomp_t i=0; i<num; ++i) r[i] = gauss_dist( generator );
    }

    //! Beta RNG: Generate beta random numbers
    //! \param[in] tid Thread (or more precisely stream) ID
    //! \param[in] num Number of RNGs to generate
    //! \param[in] p First beta shape parameter
    //! \param[in] q Second beta shape parameter
    //! \param[in] a Beta displacement parameter
    //! \param[in] b Beta scale factor
    //! \param[in,out] r Pointer to memory to write the random numbers to
    //! \details Generating beta-distributed random numbers is implemented via
    //!   an adaptor, modeling boost::UniformRandomNumberGenerator, outsourcing
    //!   the transformation of uniform random numbers to beta-distributed ones,
    //!   to boost::random. The adaptor is instantiated here because a boost
    //!   random number distribution, such as e.g.,
    //!   boost::random::beta_distribution, generates numbers using operator()
    //!   with no arguments, thus the RNG state and the thread ID (this latter
    //!   only known here) must be stored in the adaptor functor's state.
    void beta( int tid, ncomp_t num, double p, double q, double a, double b,
               double* r ) const {
      Adaptor generator( m_rng, m_data, tid );
      boost::random::beta_distribution<> beta_dist( p, q );
      for (ncomp_t i=0; i<num; ++i) r[i] = beta_dist( generator ) * b + a;
    }

    //! Accessor to the number of threads we operate on
    uint64_t nthreads() const noexcept { return m_data.size(); }

  private:
    mutable CBRNG m_rng;        //!< Random123 RNG object
    mutable arg_type m_data;    //!< RNG arguments
};

} // tk::

#endif // Random123_h
