// *****************************************************************************
/*!
  \file      src/RNG/RNGSSE.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Interface to RNGSSE random number generators
  \details   Interface to RNGSSE random number generators
*/
// *****************************************************************************
#ifndef RNGSSE_h
#define RNGSSE_h

#include <cstring>
#include <random>

#include "NoWarning/beta_distribution.h"

#include "Make_unique.h"
#include "Exception.h"
#include "Macro.h"
#include "Options/RNGSSESeqLen.h"

namespace tk {

//! RNGSSE-based random number generator used polymorphically with tk::RNG
template< class State, typename SeqNumType, unsigned int (*Generate)(State*) >
class RNGSSE {

  private:
    using InitFn = void (*)( State*, SeqNumType );
    using ncomp_t = kw::ncomp::info::expect::type;    

    //! Adaptor to use a std distribution with the RNGSSE generator
    //! \see C++ concepts: UniformRandomNumberGenerator
    struct Adaptor {
      using result_type = unsigned int;
      Adaptor( const std::unique_ptr< State[] >& s, int t ) : str(s), tid(t) {}
      static constexpr result_type min() { return 0u; }
      static constexpr result_type max() { return 4294967295u; }
      result_type operator()()
      { return Generate( &str[ static_cast<std::size_t>(tid) ] ); }
      const std::unique_ptr< State[] >& str;
      int tid;
    };

  public:
    //! Constructor
    //! \param[in] n Initialize RNG using this many independent streams
    //! \param[in] fnShort RNG initializer function for short streams
    //! \param[in] seqlen Sequence length enum: short, medium or long
    //! \param[in] fnLong RNG initializer function for long streams
    //! \param[in] fnMed RNG initializer function for medium streams
    explicit RNGSSE( SeqNumType n,
                     InitFn fnShort,
                     ctr::RNGSSESeqLenType seqlen = ctr::RNGSSESeqLenType::SHORT,
                     InitFn fnLong = nullptr,
                     InitFn fnMed = nullptr) :
       m_nthreads( n ),
       m_init( seqlen == ctr::RNGSSESeqLenType::LONG ? fnLong :
               seqlen == ctr::RNGSSESeqLenType::MEDIUM ? fnMed : fnShort ),
       m_stream()
    {
      Assert( m_init != nullptr, "nullptr passed to RNGSSE constructor" );
      Assert( n > 0, "Need at least one thread" );
      // Allocate array of stream-pointers for threads
      m_stream = tk::make_unique< State[] >( n );
      // Initialize thread-streams
      for (SeqNumType i=0; i<n; ++i) m_init( &m_stream[i], i );
    }

    //! Uniform RNG: Generate uniform random numbers
    //! \param[in] tid Thread (or more precisely) stream ID
    //! \param[in] num Number of RNGs to generate
    //! \param[in,out] r Pointer to memory to write the random numbers to
    void uniform( int tid, ncomp_t num, double* r ) const {
      for (ncomp_t i=0; i<num; ++i)
        r[i] = static_cast<double>(
                 Generate( &m_stream[ static_cast<std::size_t>(tid) ] ) )
               / 4294967296.0;
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
      Adaptor generator( m_stream, tid );
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
      Adaptor generator( m_stream, tid );
      boost::random::beta_distribution<> beta_dist( p, q );
      for (ncomp_t i=0; i<num; ++i) r[i] = beta_dist( generator ) * b + a;
    }

    //! Copy assignment
    RNGSSE& operator=( const RNGSSE& x ) {
      m_nthreads = x.m_nthreads;
      m_init = x.m_init;
      m_stream = tk::make_unique< State[] >( x.m_nthreads );
      for (SeqNumType i=0; i<x.m_nthreads; ++i) m_init( &m_stream[i], i );
      return *this;
    }

    //! Copy constructor: in terms of copy assignment
    RNGSSE( const RNGSSE& x ) { operator=(x); }

    //! Move assignment
    RNGSSE& operator=( RNGSSE&& x ) {
      m_nthreads = x.m_nthreads;
      m_init = x.m_init;
      m_stream = tk::make_unique< State[] >( x.m_nthreads );
      for (SeqNumType i=0; i<x.m_nthreads; ++i) {
        m_stream[i] = x.m_stream[i];
        std::memset( &x.m_stream[i], 0, sizeof(x.m_stream[i]) );
      }
      x.m_nthreads = 0;
      x.m_init = nullptr;
      x.m_stream.reset( nullptr );
      return *this;
    }

    //! Move constructor: in terms of move assignment
    RNGSSE( RNGSSE&& x ) :
      m_nthreads( 0 ),
      m_init( nullptr ),
      m_stream( nullptr )
    { *this = std::move( x ); }

    //! Accessor to the number of threads we operate on
    SeqNumType nthreads() const noexcept { return m_nthreads; }

  private:
    SeqNumType m_nthreads;                 //!< Number of threads
    InitFn m_init;                         //!< Sequence length initializer
    std::unique_ptr< State[] > m_stream;   //!< Random number stream for threads
};

} // tk::

#endif // RNGSSE_h
