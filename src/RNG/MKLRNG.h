//******************************************************************************
/*!
  \file      src/RNG/MKLRNG.h
  \author    J. Bakosi
  \date      Mon 15 Sep 2014 01:36:49 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************
#ifndef MKLRNG_h
#define MKLRNG_h

#include <mkl_vsl_types.h>

#include <Exception.h>

namespace tk {

//! MKL-based random number generator used polymorphically with RNG
class MKLRNG {

  public:
    //! Constructor
    explicit MKLRNG( int nthreads,
                     int brng,
                     unsigned int seed,
                     int uniform_method,
                     int gaussian_method ) :
      m_brng( brng ),
      m_seed( seed ),
      m_uniform_method( uniform_method ),
      m_gaussian_method( gaussian_method ),
      m_nthreads( nthreads ) {
      Assert( nthreads > 0, "Need at least one thread" );
      // Allocate array of stream-pointers for threads
      m_stream = tk::make_unique< VSLStreamStatePtr[] >( nthreads );
      // Initialize thread-streams for block-splitting. These MKL VSL functions
      // dynamically allocate memory, so these calls being in a constructor are
      // a potential memory leak hazard in the presence of exceptions. However,
      // thankfully, the MKL functions below only emit warnings if they
      // encounter errors and always continue. As a result, the constructor
      // finishes, the MKLRNG object gets created, so the destructor will also
      // get called when leaving scope.
      for (auto i=0; i<nthreads; ++i) {
        vslNewStream( &m_stream[i], brng, seed );
        vslLeapfrogStream( m_stream[i], i, nthreads );
      }
    }

    //! Destructor
    ~MKLRNG() { deleteStreams(); }

    //! Uniform RNG
    void uniform( int tid, int num, double* r ) const
    { vdRngUniform( m_uniform_method, m_stream[tid], num, r, 0.0, 1.0 ); }

    //! Gaussian RNG
    void gaussian( int tid, int num, double* r ) const
    { vdRngGaussian( m_gaussian_method, m_stream[tid], num, r, 0.0, 1.0 ); }

    //! Copy assignment
    MKLRNG& operator=( const MKLRNG& x ) {
      m_brng = x.m_brng;
      m_seed = x.m_seed;
      m_uniform_method = x.m_uniform_method;
      m_gaussian_method = x.m_gaussian_method;
      m_nthreads = x.m_nthreads;
      m_stream = tk::make_unique< VSLStreamStatePtr[] >( x.m_nthreads );
      for (auto i=0; i<x.m_nthreads; ++i) {
        vslNewStream( &m_stream[i], x.m_brng, x.m_seed );
        vslLeapfrogStream( m_stream[i], i, x.m_nthreads );
      }
      return *this;
    }

    //! Copy constructor: in terms of copy assignment
    MKLRNG( const MKLRNG& x ) { operator=(x); }

    //! Move assignment
    MKLRNG& operator=( MKLRNG&& x ) {
      deleteStreams();
      m_brng = x.m_brng;
      m_seed = x.m_seed;
      m_uniform_method = x.m_uniform_method;
      m_gaussian_method = x.m_gaussian_method;
      m_nthreads = x.m_nthreads;
      m_stream = tk::make_unique< VSLStreamStatePtr[] >( x.m_nthreads );
      for (auto i=0; i<x.m_nthreads; ++i) {
        m_stream[i] = x.m_stream[i];
        x.m_stream[i] = nullptr;
      }
      x.m_brng = 0;
      x.m_seed = 0;
      x.m_uniform_method = 0;
      x.m_gaussian_method = 0;
      x.m_nthreads = 0;
      x.m_stream.reset( nullptr );
      return *this;
    }

    //! Move constructor: in terms of move assignment
    MKLRNG( MKLRNG&& x ) :
      m_brng( 0 ),
      m_seed( 0 ),
      m_uniform_method( 0 ),
      m_gaussian_method( 0 ),
      m_nthreads( 0 ),
      m_stream( nullptr )
    { *this = std::move(x); }

  private:
    //! Delete all thread streams
    void deleteStreams() {
      for (auto i=0; i<m_nthreads; ++i)
        if (m_stream[i]) {
          vslDeleteStream( &m_stream[i] );
          m_stream[i] = nullptr;
        }
    }

    int m_brng;                                      //!< MKL RNG id
    unsigned int m_seed;                             //!< Seed
    int m_uniform_method;                            //!< Uniform method to use
    int m_gaussian_method;                           //!< Gaussian method to use
    int m_nthreads;                                  //!< Number of threads
    std::unique_ptr< VSLStreamStatePtr[] > m_stream; //! Random number streams
};

} // tk::

#endif // MKLRNG_h
