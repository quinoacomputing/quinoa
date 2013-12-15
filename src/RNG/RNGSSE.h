//******************************************************************************
/*!
  \file      src/RNG/RNGSSE.h
  \author    J. Bakosi
  \date      Sat 14 Dec 2013 09:35:35 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGSSE-based random number generator
  \details   RNGSSE-based random number generator
*/
//******************************************************************************
#ifndef RNGSSE_h
#define RNGSSE_h

#include <RNG.h>
#include <Exception.h>

namespace quinoa {

//! RNGSSE-based random number generator
template< class State,
          typename SeqNumType,
          void (*Init)(State*, SeqNumType),
          unsigned int (*Generate)(State*) > 
class RNGSSE : public tk::RNG {

  public:
    //! Constructor
    explicit RNGSSE( SeqNumType nthreads, unsigned int seed )
    {
      // Throw if not NDEBUG and nthreads invalid
      Assert(nthreads > 0, tk::ExceptType::FATAL, "Need at least one thread");
      // Allocate array of stream-pointers for threads
      m_stream = std::unique_ptr< State[] >( new State [nthreads] );
      // Initialize thread-streams
      for ( SeqNumType i=0; i<nthreads; ++i) {
        Init( &m_stream[i], i );
      }
    }

    //! Destructor
    // ICC: should be '= default'
    ~RNGSSE() noexcept override {}

    //! Uniform RNG
    void uniform( int tid, int num, double* r) const override {
      r[0] = static_cast<double>( Generate( &m_stream[tid] ) ) / 4294967296.0;
    }

    //! Gaussian RNG
    void gaussian(int tid, int num, double* r) const override {}

  private:
    //! Don't permit copy constructor
    RNGSSE(const RNGSSE&) = delete;
    //! Don't permit copy assigment
    RNGSSE& operator=(const RNGSSE&) = delete;
    //! Don't permit move constructor
    RNGSSE(RNGSSE&&) = delete;
    //! Don't permit move assigment
    RNGSSE& operator=(RNGSSE&&) = delete;

    //! Random number stream for threads
    std::unique_ptr< State[] > m_stream;
};

} // quinoa::

#endif // RNGSSE_h
