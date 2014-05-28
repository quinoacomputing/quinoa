//******************************************************************************
/*!
  \file      src/RNG/RNGSSE.h
  \author    J. Bakosi
  \date      Wed 28 May 2014 12:17:03 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGSSE-based random number generator
  \details   RNGSSE-based random number generator
*/
//******************************************************************************
#ifndef RNGSSE_h
#define RNGSSE_h

#include <cstring>

#include <make_unique.h>
#include <Exception.h>
#include <Macro.h>

namespace tk {

//! RNGSSE-based random number generator used polymorphically with RNG
template< class State, typename SeqNumType, unsigned int (*Generate)(State*) >
class RNGSSE {

    using InitFn = void (*)(State*, SeqNumType);

  public:
    //! Constructor with sequence length option: short, long, medium
    explicit RNGSSE( SeqNumType nthreads,
                     InitFn fnShort,
                     ctr::RNGSSESeqLenType seqlen = ctr::RNGSSESeqLenType::SHORT,
                     InitFn fnLong = nullptr,
                     InitFn fnMed = nullptr) :
       m_nthreads( nthreads ),
       m_init( seqlen == ctr::RNGSSESeqLenType::LONG ? fnLong :
               seqlen == ctr::RNGSSESeqLenType::MEDIUM ? fnMed : fnShort ) {
      Assert( m_init != nullptr,
              ExceptType::FATAL, "nullptr passed to RNGSSE constructor" );
      Assert( nthreads > 0, ExceptType::FATAL, "Need at least one thread" );
      // Allocate array of stream-pointers for threads
      m_stream = make_unique< State[] >( nthreads );
      // Initialize thread-streams
      for (SeqNumType i=0; i<nthreads; ++i) m_init( &m_stream[i], i );
    }

    //! Uniform RNG
    void uniform( int tid, int num, double* r ) const {
      for (int i=0; i<num; ++i)
        r[i] = static_cast<double>( Generate( &m_stream[tid] ) ) / 4294967296.0;
    }

    //! Gaussian RNG
    void gaussian( int tid, int num, double* r ) const {
      Throw( ExceptType::WARNING, "RNGSSE::gaussian undefined" );
      IGNORE(tid);
      IGNORE(num);
      IGNORE(r);
    }

    //! Copy assignment
    RNGSSE& operator=( const RNGSSE& x ) {
      m_nthreads = x.m_nthreads;
      m_init = x.m_init;
      m_stream = make_unique< State[] >( x.m_nthreads );
      for (SeqNumType i=0; i<x.m_nthreads; ++i) m_init( &m_stream[i], i );
      return *this;
    }

    //! Copy constructor: in terms of copy assignment
    RNGSSE( const RNGSSE& x ) { operator=(x); }

    //! Move assignment
    RNGSSE& operator=( RNGSSE&& x ) {
      m_nthreads = x.m_nthreads;
      m_init = x.m_init;
      m_stream = make_unique< State[] >( x.m_nthreads );
      for (SeqNumType i=0; i<x.m_nthreads; ++i) {
        m_stream[i] = x.m_stream[i];
        memset( &x.m_stream[i], 0, sizeof(x.m_stream[i]) );
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
