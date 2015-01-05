//******************************************************************************
/*!
  \file      src/DiffEq/WrightFisher.h
  \author    J. Bakosi
  \date      Tue 13 Jan 2015 11:00:25 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Wright-Fisher SDE
  \details   Wright-Fisher SDE, see
             http://www.sciencedirect.com/science/article/pii/S0040580912001013
*/
//******************************************************************************
#ifndef WrightFisher_h
#define WrightFisher_h

#include <numeric>

#ifdef HAS_MKL
  #include <mkl_lapacke.h>
#else
  #include <lapacke.h>
#endif

#include <InitPolicy.h>
#include <WFCoeffPolicy.h>
#include <RNG.h>

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! Wright-Fisher SDE used polymorphically with DiffEq
template< class Init, class Coefficients >
class WrightFisher {

  void print_matrix(const char *name, const double *mat, int n, int m) const {
    int i, j;
    printf("%s:\n", name);
    for(i=0; i<n; ++i) {
        for(j=0; j<m; ++j)
            printf("%10g ", mat[i*n+j]);
        printf("\n");
    }
    printf("\n");
  }

  public:
    //! Constructor: use coefficients policy to initialize coefficients
    explicit WrightFisher( unsigned int c ) :
      m_ncomp( g_inputdeck.get< tag::component >().
                           get< tag::wrightfisher >()[c] ),
      m_offset( g_inputdeck.get< tag::component >().
                            offset< tag::wrightfisher >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::wrightfisher, tag::rng >()[c] ) ) )
    {
      Throw( "Wright-Fisher diffusion matrix not yet implemented! See comments "
             "in code for details." );
      const auto& omega =
        g_inputdeck.get< tag::param, tag::wrightfisher, tag::omega >();
      ErrChk( omega.size() > c,
              "Indexing out of Wright-Fisher SDE parameters 'omega'");
      // Use coefficients policy to initialize coefficients
      Coefficients( m_ncomp, omega[c], m_omega );
    }

    //! Set initial conditions
    void initialize( tk::ParProps& particles ) const {
      //Init( { particles } );
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Initialize the first m_ncomp (N-1) scalars
        tk::ctr::ncomp_type i;
        for (i=0; i<m_ncomp-1; ++i) {
          particles( p, i, m_offset ) = (1.0+i)/m_ncomp;
        }
        // Initialize the (N-1)th scalar from unit-sum
        tk::real& par = particles( p, i, m_offset );
        par = 1.0 - particles( p, 0, m_offset );
        for (i=1; i<m_ncomp-1; ++i) {
          par -= particles( p, i, m_offset );
        }
      }
    }

    //! Advance particles
    void advance( tk::ParProps& particles, int stream, tk::real dt ) const {
      // Compute sum of coefficients
      const auto omega = std::accumulate( begin(m_omega), end(m_omega), 0.0 );
      const auto npar = particles.npar();

      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Need to build the square-root of the Wright-Fisher diffusion matrix:
        // B_ij = y_i * ( delta_ij - y_j ). If the matrix is positive definite,
        // the Cholesky decomposition would work, however, B_ij is only positive
        // semi-definite, so Cholesky may fail and considering floating-point
        // errors it will definitely fail at some point.
        //
        // A stable square-root for B_ij is not yet implemented. Here is some
        // details and points for further development:
        //
        // An interesting fact from http://math.stackexchange.com/a/332465, on
        // gauging the eigen values of a matrix:
        //
        // "Choosing each diagonal entry to be greater than the sum of the
        // absolute values of the other entries in the same row will immediately
        // imply that all of the eigenvalues of A are positive, and therefore
        // that A is positive definite."
        //
        // The WF diffusion matrix seems to be positive semi-definite with at
        // least one very small eigenvalue O(10e-18). The row and column sums
        // are zero. It does not appear that the matrix has negative
        // eigenvalues.
        //
        // Using Cholesky decomposition to compute the square-root may fail.
        // Other methods, such as diagonalizing the matrix or the LDL
        // decomposition, may work, but have not tried. The former requires
        // finding the eigenvalues, with e.g., LAPACK's DSPEVD, thereby
        // diagonalizing the matrix, and taking the square-root of the diagonal
        // elements to find the square-root of the matrix. See
        // http://en.wikipedia.org/wiki/Square_root_of_a_matrix.  This may be
        // okay, but probably overkill, as we don't really need the eigenvalues,
        // only the square-root of the matrix.
        //
        // The latter (LDL) is similar to Cholesky, but does not involve square
        // root. See http://en.wikipedia.org/wiki/Cholesky_decomposition. See
        // also Golub van Loan: Sec.4.1.2 LDL. It appears that LDL would compare
        // in perfromance to Cholesky and would work even on large indefinite
        // matrices. LAPACK probably has LDL.
        //
        // I'm putting this aside for now. The next step should probably be
        // idintifying the LAPACK for LDL and use it to advance the WF system
        // using LD^{1/2}.
        //
        // See also:
        // http://en.wikipedia.org/wiki/Positive-definite_matrix
        // file:///opt/intel/composerxe/Documentation/en_US/mkl/mklman/index.htm
        // http://icl.cs.utk.edu/lapack-forum/archives/lapack/msg01497.html
        // http://yarchive.net/comp/sqrtm.html

        tk::real B[m_ncomp][m_ncomp];
        tk::real Bo[m_ncomp][m_ncomp];
        for (tk::ctr::ncomp_type i=0; i<m_ncomp; ++i) {
          const tk::real& pari = particles( p, i, m_offset );
          for (tk::ctr::ncomp_type j=0; j<m_ncomp; ++j) {
            const tk::real& parj = particles( p, j, m_offset );
            if (i == j) {
              B[i][i] = std::abs( pari * (1.0 - pari) );
              if (B[i][i] < 1.0e-10) B[i][i] = 1.0;
            } else {
              B[i][j] = -pari*parj;
            }
          }
        }
        std::memcpy( Bo, B, m_ncomp*m_ncomp*sizeof(tk::real) );
        // Compute diffusion matrix (lower triangle of Cholesky-decomposition)
        lapack_int n = static_cast< lapack_int >( m_ncomp );
        lapack_int info = LAPACKE_dpotrf( LAPACK_ROW_MAJOR, 'L', n, B[0], n );
        if (info != 0 ) {
          print_matrix( "=======\nOriginal Matrix", Bo[0], n, n );
          std::cout << "info: " << info << std::endl;
          print_matrix( "Result of Cholesky factorization", B[0], n, n );
          for (tk::ctr::ncomp_type i=0; i<m_ncomp; ++i) {
            std::cout <<
              std::setprecision( std::numeric_limits< tk::real >::digits10 )
            << i << " par: " << particles( p, i, m_offset) << std::endl;
          }
        }
        //ErrChk( info == 0, "Wright-Fisher Cholesky decomposition unsuccessful, "
        //                   "info = " + std::to_string( info ) + ", particle: " +
        //                   std::to_string( p ) );

        // Advance the first m_ncomp (N-1) scalars
        if (info == 0) {
          tk::ctr::ncomp_type i = 0;
          for (i=0; i<m_ncomp-1; ++i) {
            tk::real& par = particles( p, i, m_offset );
            // Advance first m_ncomp (K=N-1) scalars due to drift
            par += 0.5*(m_omega[i] - omega*par)*dt;
            // Advance first m_ncomp (K=N-1) particles with Cholesky-decomposed
            // lower triangle (diffusion matrix)
            for (tk::ctr::ncomp_type j=0; j<m_ncomp-1; ++j)
              if (j<=i) {
                tk::real dW;
                m_rng.gaussian( stream, 1, &dW );
                par += B[i][j] * sqrt(dt) * dW;
              }
          }
          // Compute the (N-1)th scalar from unit-sum
          tk::real& par = particles( p, i, m_offset );
          par = 1.0 - particles( p, 0, m_offset );
          for (i=1; i<m_ncomp-1; ++i) {
            par -= particles( p, i, m_offset );
          }
        }
      }
    }

  private:
    const tk::ctr::ncomp_type m_ncomp;  //!< Number of components
    const int m_offset;                 //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator
    std::vector< kw::sde_omega::info::expect::type > m_omega;  //!< Coefficients
};

} // walker::

#endif // WrightFisher_h
