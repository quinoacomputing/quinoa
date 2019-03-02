// *****************************************************************************
/*!
  \file      src/RNG/MKLRNG.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Interface to Intel MKL VSL random number generators
  \details   Interface to Intel MKL VSL random number generators.
*/
// *****************************************************************************
#ifndef MKLRNG_h
#define MKLRNG_h

#include <mkl_vsl.h>

#include "Exception.h"
#include "Make_unique.h"
#include "Keywords.h"

namespace tk {

//! MKL-based random number generator used polymorphically with tk::RNG
class MKLRNG {

    using ncomp_t = kw::ncomp::info::expect::type;    

  public:
    //! Constructor
    //! \param[in] n Initialize RNG using this many independent streams
    //! \param[in] brng Index of the basic generator to initialize the stream
    //! \param[in] seed RNG seed
    //! \param[in] uniform_method MKL ID of the method to use for uniform RNGs
    //! \param[in] gaussian_method MKL ID of the method to use for Gaussian RNGs
    //! \param[in] gaussianmv_method MKL ID of the method to use for
    //!    multi-variate Gaussian RNGs
    //! \param[in] beta_method MKL ID of the method to use for beta RNGs
    //! \param[in] gamma_method MKL ID of the method to use for gamma RNGs
    explicit MKLRNG(
      int n = 1,
      int brng = VSL_BRNG_MCG59,
      unsigned int seed = 0,
      int uniform_method = VSL_RNG_METHOD_UNIFORM_STD,
      int gaussian_method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
      int gaussianmv_method = VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2,
      int beta_method = VSL_RNG_METHOD_BETA_CJA,
      int gamma_method = VSL_RNG_METHOD_GAMMA_GNORM ) :
      m_brng( brng ),
      m_seed( seed ),
      m_uniform_method( uniform_method ),
      m_gaussian_method( gaussian_method ),
      m_gaussianmv_method( gaussianmv_method ),
      m_beta_method( beta_method ),
      m_gamma_method( gamma_method ),
      m_nthreads( n ),
      m_stream()
    {
      Assert( n > 0, "Need at least one thread" );
      Assert( brng > 0, "Basic RNG MKL parameter must be positive" );
      // Allocate array of stream-pointers for threads
      m_stream = tk::make_unique< VSLStreamStatePtr[] >(
                   static_cast<std::size_t>(n) );
      // Initialize thread-streams for block-splitting. These MKL VSL functions
      // dynamically allocate memory, so these calls being in a constructor are
      // a potential memory leak hazard in the presence of exceptions. However,
      // thankfully, the MKL functions below only emit warnings if they
      // encounter errors and always continue. As a result, the constructor
      // finishes, the MKLRNG object gets created, so the destructor will also
      // get called when leaving scope.
      if (n == 1)
        errchk( vslNewStream( &m_stream[0], brng, seed ) );
      else
        for (int i=0; i<n; ++i) {
          auto I = static_cast< std::size_t >( i );
          errchk( vslNewStream( &m_stream[I], brng, seed ) );
          errchk( vslLeapfrogStream( m_stream[I], i, n ) );
        }
    }

    //! Destructor
    ~MKLRNG() noexcept { deleteStreams(); }

    //! Uniform RNG: Generate uniform random numbers
    //! \param[in] tid Thread (or more precisely stream) ID
    //! \param[in] num Number of RNGs to generate
    //! \param[in,out] r Pointer to memory to write the random numbers to
    void uniform( int tid, ncomp_t num, double* r ) const {
      vdRngUniform( m_uniform_method,
                    m_stream[ static_cast<std::size_t>(tid) ],
                    static_cast< long long >( num ),
                    r,
                    0.0, 1.0 );
    }

    //! Gaussian RNG: Generate Gaussian random numbers
    //! \param[in] tid Thread (or more precisely stream) ID
    //! \param[in] num Number of RNGs to generate
    //! \param[in,out] r Pointer to memory to write the random numbers to
    void gaussian( int tid, ncomp_t num, double* r ) const {
      vdRngGaussian( m_gaussian_method,
                     m_stream[ static_cast<std::size_t>(tid) ],
                     static_cast< long long >( num ),
                     r,
                     0.0, 1.0 );
    }

    //! \brief Multi-variate Gaussian RNG: Generate multi-variate Gaussian
    //!    random numbers
    //! \param[in] tid Thread (or more precisely stream) ID
    //! \param[in] num Number of RNGs to generate
    //! \param[in] d Dimension d ( d ≥ 1) of output random vectors
    //! \param[in] mean Mean vector of dimension d
    //! \param[in] cov Lower triangle of covariance matrix, stored as a vector
    //!   of length d(d+1)/2
    //! \param[in,out] r Pointer to memory to write the random numbers to
    void gaussianmv( int tid, ncomp_t num, ncomp_t d, const double* const mean,
                     const double* const cov, double* r ) const
    {
      Assert( d > 0,
              "Dimension of multi-variate Gaussian RNGs must be positive" );
      vdRngGaussianMV( m_gaussianmv_method,
                       m_stream[ static_cast<std::size_t>(tid) ],
                       static_cast< long long >( num ),
                       r,
                       static_cast< int >( d ),
                       VSL_MATRIX_STORAGE_PACKED,
                       mean,
                       cov );
    }

    //! Beta RNG: Generate beta random numbers
    //! \param[in] tid Thread (or more precisely stream) ID
    //! \param[in] num Number of RNGs to generate
    //! \param[in] p First beta shape parameter
    //! \param[in] q Second beta shape parameter
    //! \param[in] a Beta displacement parameter
    //! \param[in] b Beta scale factor
    //! \param[in,out] r Pointer to memory to write the random numbers to
    void beta( int tid, ncomp_t num, double p, double q, double a, double b,
               double* r ) const
    {
      vdRngBeta( m_beta_method,
                 m_stream[ static_cast<std::size_t>(tid) ],
                 static_cast< long long >( num ),
                 r,
                 p, q, a, b );
    }

    //! Gamma RNG: Generate gamma random numbers
    //! \param[in] tid Thread (or more precisely stream) ID
    //! \param[in] num Number of RNGs to generate
    //! \param[in] a Gamma shape parameter
    //! \param[in] b Gamma scale factor
    //! \param[in,out] r Pointer to memory to write the random numbers to
    void gamma( int tid, ncomp_t num, double a, double b, double* r ) const
    {
      vdRngGamma( m_beta_method,
                  m_stream[ static_cast<std::size_t>(tid) ],
                  static_cast< long long >( num ),
                  r,
                  a, 0.0, b );  // displacement = 0.0
    }

    //! Copy assignment
    MKLRNG& operator=( const MKLRNG& x ) {
      m_brng = x.m_brng;
      m_seed = x.m_seed;
      m_uniform_method = x.m_uniform_method;
      m_gaussian_method = x.m_gaussian_method;
      m_gaussianmv_method = x.m_gaussianmv_method;
      m_beta_method = x.m_beta_method;
      m_gamma_method = x.m_gamma_method;
      m_nthreads = x.m_nthreads;
      m_stream = tk::make_unique< VSLStreamStatePtr[] >(
                   static_cast<std::size_t>(x.m_nthreads) );
      if (m_nthreads == 1)
        errchk( vslNewStream( &m_stream[0], x.m_brng, x.m_seed ) );
      else
        for (int i=0; i<x.m_nthreads; ++i) {
          auto I = static_cast< std::size_t >( i );
          errchk( vslNewStream( &m_stream[I], x.m_brng, x.m_seed ) );
          errchk( vslLeapfrogStream( m_stream[I], i, x.m_nthreads ) );
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
      m_gaussianmv_method = x.m_gaussianmv_method;
      m_beta_method = x.m_beta_method;
      m_gamma_method = x.m_gamma_method;
      m_nthreads = x.m_nthreads;
      m_stream = tk::make_unique< VSLStreamStatePtr[] >(
                   static_cast<std::size_t>(x.m_nthreads) );
      for (int i=0; i<x.m_nthreads; ++i) {
        auto I = static_cast< std::size_t >( i );
        m_stream[I] = x.m_stream[I];
        x.m_stream[I] = nullptr;
      }
      x.m_brng = 0;
      x.m_seed = 0;
      x.m_uniform_method = 0;
      x.m_gaussian_method = 0;
      x.m_gaussianmv_method = 0;
      x.m_beta_method = 0;
      x.m_gamma_method = 0;
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
      m_gaussianmv_method( 0 ),
      m_beta_method( 0 ),
      m_gamma_method( 0 ),
      m_nthreads( 0 ),
      m_stream( nullptr )
    { *this = std::move(x); }

    //! Accessor to the number of threads we operate on
    std::size_t nthreads() const noexcept
    { return static_cast< std::size_t >( m_nthreads); }

  private:
    //! Delete all thread streams
    void deleteStreams() {
      for (int i=0; i<m_nthreads; ++i) {
        auto I = static_cast< std::size_t >( i );
        if (m_stream[I]) {
          vslDeleteStream( &m_stream[I] );
          m_stream[I] = nullptr;
        }
      }
    }

    //! MKL VSL error check
    //! \param[in] err MKL VSL error code as returned from MKL VSL functions
    //! \details This calls ErrChk(), i.e., it is not compiled away in Release
    //!   mode as an error here can result due to user input incompatible with
    //!   the MKL library.
    void errchk( int err ) {
      ErrChk( err == VSL_STATUS_OK, "MKL VSL Error Code: " +
              std::to_string(err) + ", see mkl_vsl_defines.h for more info" );
    }

    int m_brng;                                      //!< MKL RNG id
    unsigned int m_seed;                             //!< Seed
    int m_uniform_method;                            //!< Uniform method to use
    int m_gaussian_method;                           //!< Gaussian method to use
    int m_gaussianmv_method;           //!< Multi-variate Gaussian method to use
    int m_beta_method;                               //!< Beta method to use
    int m_gamma_method;                              //!< Gamma method to use
    int m_nthreads;                                  //!< Number of threads
    std::unique_ptr< VSLStreamStatePtr[] > m_stream; //!< Random number streams
};

} // tk::

#endif // MKLRNG_h
