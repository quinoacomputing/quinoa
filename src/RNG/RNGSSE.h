//******************************************************************************
/*!
  \file      src/RNG/RNGSSE.h
  \author    J. Bakosi
  \date      Fri 01 May 2015 07:31:28 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Interface to RNGSSE random number generators
  \details   Interface to RNGSSE random number generators
*/
//******************************************************************************
#ifndef RNGSSE_h
#define RNGSSE_h

#include <cstring>

#include <make_unique.h>
#include <Exception.h>
#include <Macro.h>

namespace tk {

//! RNGSSE-based random number generator used polymorphically with tk::RNG
template< class State, typename SeqNumType, unsigned int (*Generate)(State*) >
class RNGSSE {

    using InitFn = void (*)( State*, SeqNumType );
    using ncomp_t = kw::ncomp::info::expect::type;    

  public:
    //! \brief Constructor
    //! \param[in] nthreads Initialize RNG using this many independent streams
    //! \param[in] fnShort RNG initializer function for short streams
    //! \param[in] seqlen Sequence length enum: short, medium or long
    //! \param[in] fnLong RNG initializer function for long streams
    //! \param[in] fnMed RNG initializer function for medium streams
    explicit RNGSSE( SeqNumType nthreads,
                     InitFn fnShort,
                     ctr::RNGSSESeqLenType seqlen = ctr::RNGSSESeqLenType::SHORT,
                     InitFn fnLong = nullptr,
                     InitFn fnMed = nullptr) :
       m_nthreads( nthreads ),
       m_init( seqlen == ctr::RNGSSESeqLenType::LONG ? fnLong :
               seqlen == ctr::RNGSSESeqLenType::MEDIUM ? fnMed : fnShort ) {
      Assert( m_init != nullptr, "nullptr passed to RNGSSE constructor" );
      Assert( nthreads > 0, "Need at least one thread" );
      // Allocate array of stream-pointers for threads
      m_stream = tk::make_unique< State[] >( nthreads );
      // Initialize thread-streams
      for (SeqNumType i=0; i<nthreads; ++i) m_init( &m_stream[i], i );
    }

    //! Uniform RNG: Generate uniform random numbers
    //! \param[in] tid Thread (or more precisely) stream ID
    //! \param[in] num Number of RNGs to generate
    //! \param[inout] r Pointer to memory to write the RNGs to
    void uniform( int tid, ncomp_t num, double* r ) const {
      for (int i=0; i<num; ++i)
        r[i] = static_cast<double>(
                 Generate( &m_stream[ static_cast<std::size_t>(tid) ] ) )
               / 4294967296.0;
    }

    //! Gaussian RNG: Generate Gaussian random numbers
    //! \param[in] tid Thread (or more precisely stream) ID
    //! \param[in] num Number of RNGs to generate
    //! \param[inout] r Pointer to memory to write the RNGs to
    // TODO: not yet implemented
    void gaussian( int tid, ncomp_t num, double* r ) const {
      Throw( "RNGSSE::gaussian undefined" );
      IGNORE(tid);
      IGNORE(num);
      IGNORE(r);
    }

    //! Beta RNG: Generate beta random numbers
    //! \param[in] tid Thread (or more precisely stream) ID
    //! \param[in] num Number of RNGs to generate
    //! \param[in] p First beta shape parameter
    //! \param[in] q Second beta shape parameter
    //! \param[in] a Beta displacement parameter
    //! \param[in] b Beta scale factor
    //! \param[inout] r Pointer to memory to write the RNGs to
    // TODO: not yet implemented
    void beta( int tid, ncomp_t num, double p, double q, double a, double b,
               double* r ) const {
      Throw( "RNGSSE::beta undefined" );
      IGNORE(tid);
      IGNORE(num);
      IGNORE(p);
      IGNORE(q);
      IGNORE(a);
      IGNORE(b);
      IGNORE(r);
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

  private:
    SeqNumType m_nthreads;                 //!< Number of threads
    InitFn m_init;                         //!< Sequence length initializer
    std::unique_ptr< State[] > m_stream;   //!< Random number stream for threads
};

} // tk::

#endif // RNGSSE_h
