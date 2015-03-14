//******************************************************************************
/*!
  \file      src/RNG/MKLRNG.h
  \author    J. Bakosi
  \date      Thu 12 Mar 2015 09:58:40 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Interface to Intel MKL VSL random number generators
  \details   Interface to Intel MKL VSL random number generators.
*/
//******************************************************************************
#ifndef MKLRNG_h
#define MKLRNG_h

#include <mkl_vsl_types.h>

#include <Exception.h>

namespace tk {

//! MKL-based random number generator used polymorphically with tk::RNG
class MKLRNG {

    using ncomp_t = kw::ncomp::info::expect::type;    

  public:
    //! \brief Constructor
    //! \param[in] nthreads Initialize RNG using this many independent streams
    //! \param[in] brng Index of the basic generator to initialize the stream
    //! \param[in] seed RNG seed
    //! \param[in] uniform_method MKL ID of the method to use for uniform RNGs
    //! \param[in] gaussian_method MKL ID of the method to use for Gaussian RNGs
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
      m_stream = tk::make_unique< VSLStreamStatePtr[] >(
                   static_cast<std::size_t>(nthreads) );
      // Initialize thread-streams for block-splitting. These MKL VSL functions
      // dynamically allocate memory, so these calls being in a constructor are
      // a potential memory leak hazard in the presence of exceptions. However,
      // thankfully, the MKL functions below only emit warnings if they
      // encounter errors and always continue. As a result, the constructor
      // finishes, the MKLRNG object gets created, so the destructor will also
      // get called when leaving scope.
      for (int i=0; i<nthreads; ++i) {
        auto I = static_cast< std::size_t >( i );
        vslNewStream( &m_stream[ I ], brng, seed );
        vslLeapfrogStream( m_stream[ I ], i, nthreads );
      }
    }

    //! Destructor
    ~MKLRNG() { deleteStreams(); }

    //! Uniform RNG: Generate uniform random numbers
    //! \param[in] tid Thread (or more precisely) stream ID
    //! \param[in] num Number of RNGs to generate
    //! \param[inout r Pointer to memory to write the RNGs to
    void uniform( int tid, ncomp_t num, double* r ) const {
      vdRngUniform( m_uniform_method,
                    m_stream[ static_cast<std::size_t>(tid) ],
                    static_cast< long long >( num ),
                    r,
                    0.0,
                    1.0 );
    }

    //! Gaussian RNG: Generate Gaussian random numbers
    //! \param[in] tid Thread (or rather) stream ID
    //! \param[in] num Number of RNGs to generate
    //! \param[inout r Pointer to memory to write the RNGs to
    void gaussian( int tid, ncomp_t num, double* r ) const {
      vdRngGaussian( m_gaussian_method,
                     m_stream[ static_cast<std::size_t>(tid) ],
                     static_cast< long long >( num ),
                     r,
                     0.0,
                     1.0 );
    }

    //! Copy assignment
    MKLRNG& operator=( const MKLRNG& x ) {
      m_brng = x.m_brng;
      m_seed = x.m_seed;
      m_uniform_method = x.m_uniform_method;
      m_gaussian_method = x.m_gaussian_method;
      m_nthreads = x.m_nthreads;
      m_stream = tk::make_unique< VSLStreamStatePtr[] >(
                   static_cast<std::size_t>(x.m_nthreads) );
      for (int i=0; i<x.m_nthreads; ++i) {
        auto I = static_cast< std::size_t >( i );
        vslNewStream( &m_stream[ I ], x.m_brng, x.m_seed );
        vslLeapfrogStream( m_stream[ I ], i, x.m_nthreads );
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
      m_stream = tk::make_unique< VSLStreamStatePtr[] >(
                   static_cast<std::size_t>(x.m_nthreads) );
      for (int i=0; i<x.m_nthreads; ++i) {
        auto I = static_cast< std::size_t >( i );
        m_stream[ I ] = x.m_stream[ I ];
        x.m_stream[ I ] = nullptr;
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
      for (int i=0; i<m_nthreads; ++i) {
        auto I = static_cast< std::size_t >( i );
        if (m_stream[ I ]) {
          vslDeleteStream( &m_stream[ I ] );
          m_stream[ I ] = nullptr;
        }
      }
    }

    int m_brng;                                      //!< MKL RNG id
    unsigned int m_seed;                             //!< Seed
    int m_uniform_method;                            //!< Uniform method to use
    int m_gaussian_method;                           //!< Gaussian method to use
    int m_nthreads;                                  //!< Number of threads
    std::unique_ptr< VSLStreamStatePtr[] > m_stream; //!< Random number streams
};

} // tk::

#endif // MKLRNG_h
