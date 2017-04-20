// *****************************************************************************
/*!
  \file      src/DiffEq/WrightFisher.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Wright-Fisher SDE
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs), whose invariant is the Dirichlet
    distribution. For more details on the Wright-Fisher SDE, see
    http://www.sciencedirect.com/science/article/pii/S0040580912001013.
*/
// *****************************************************************************
#ifndef WrightFisher_h
#define WrightFisher_h

#include <numeric>
#include <vector>
#include <iomanip>

#include "QuinoaConfig.h"

#ifdef HAS_MKL
  #include "NoWarning/mkl_lapacke.h"
#else
  #include "NoWarning/lapacke.h"
#endif

#include "Macro.h"
#include "InitPolicy.h"
#include "WrightFisherCoeffPolicy.h"
#include "RNG.h"
#include "Particles.h"

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief Wright-Fisher SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see DiffEq/WrightFisherCoeffPolicy.h
template< class Init, class Coefficients >
class WrightFisher {

    void print_matrix( const char *name, const double *mat, lapack_int n,
                       lapack_int m ) const
    {
      int i, j;
      printf("%s:\n", name);
      for(i=0; i<n; ++i) {
          for(j=0; j<m; ++j)
              printf("%10g ", mat[i*n+j]);
          printf("\n");
      }
      printf("\n");
    }

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of Wright-Fisher SDEs to
    //!   construct. There can be multiple wright-fisher ... end blocks in a
    //!   control file. This index specifies which Wright-Fisher SDE system to
    //!   instantiate. The index corresponds to the order in which the
    //!   wright-fisher ... end blocks are given the control file.
    //! \author J. Bakosi
    explicit WrightFisher( ncomp_t c ) :
      m_c( c ),
      m_depvar(
        g_inputdeck.get< tag::param, tag::wrightfisher, tag::depvar >().at(c) ),
      m_ncomp(
        g_inputdeck.get< tag::component >().get< tag::wrightfisher >().at(c) ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< tag::wrightfisher >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::wrightfisher, tag::rng >().at(c) ) ) ),
      m_omega(),
      coeff(
        m_ncomp,
        g_inputdeck.get< tag::param, tag::wrightfisher, tag::omega >().at(c),
        m_omega )
    {
      Throw( "Wright-Fisher diffusion matrix not yet implemented! See comments "
             "in code for details." );
    }
    //! Initalize SDE, prepare for time integration
    //! \param[in] stream Thread (or more precisely stream) ID 
    //! \param[in,out] particles Array of particle properties 
    //! \author J. Bakosi
    void initialize( int stream, tk::Particles& particles ) {
      IGNORE( stream );
      //! Set initial conditions using initialization policy
      //Init::template
      //  init< tag::wrightfisher >
      //      ( g_inputdeck, m_rng, stream, particles, m_c, m_ncomp, m_offset );

      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Initialize the first m_ncomp (N-1) scalars
        ncomp_t i;
        for (i=0; i<m_ncomp-1; ++i) {
          particles( p, i, m_offset ) =
            (1.0 + static_cast<tk::real>(i)) / static_cast<tk::real>(m_ncomp);
        }
        // Initialize the (N-1)th scalar from unit-sum
        tk::real& par = particles( p, i, m_offset );
        par = 1.0 - particles( p, 0, m_offset );
        for (i=1; i<m_ncomp-1; ++i) {
          par -= particles( p, i, m_offset );
        }
      }
    }

    //! \brief Advance particles according to the Wright-Fisher SDE
    //! \param[in,out] particles Array of particle properties
    //! \param[in] stream Thread (or more precisely stream) ID
    //! \param[in] dt Time step size
    //! \author J. Bakosi
    void advance( tk::Particles& particles,
                  int stream,
                  tk::real dt,
                  tk::real,
                  const std::map< tk::ctr::Product, tk::real >& )
    {
      // Compute sum of coefficients
      const auto omega = std::accumulate( begin(m_omega), end(m_omega), 0.0 );
      const auto npar = particles.nunk();

      #if defined(__clang__)
        #pragma clang diagnostic push
        #pragma clang diagnostic ignored "-Wvla"
        #pragma clang diagnostic ignored "-Wvla-extension"
      #elif defined(STRICT_GNUC)
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wvla"
      #endif

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
        for (ncomp_t i=0; i<m_ncomp; ++i) {
          const tk::real& pari = particles( p, i, m_offset );
          for (ncomp_t j=0; j<m_ncomp; ++j) {
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
          for (ncomp_t i=0; i<m_ncomp; ++i) {
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
          ncomp_t i = 0;
          for (i=0; i<m_ncomp-1; ++i) {
            tk::real& par = particles( p, i, m_offset );
            // Advance first m_ncomp (K=N-1) scalars due to drift
            par += 0.5*(m_omega[i] - omega*par)*dt;
            // Advance first m_ncomp (K=N-1) particles with Cholesky-decomposed
            // lower triangle (diffusion matrix)
            for (ncomp_t j=0; j<m_ncomp-1; ++j)
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

      #if defined(__clang__)
        #pragma clang diagnostic pop
      #elif defined(STRICT_GNUC)
        #pragma GCC diagnostic pop
      #endif
    }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const char m_depvar;                //!< Dependent variable
    const ncomp_t m_ncomp;              //!< Number of components
    const ncomp_t m_offset;             //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator

    //! Coefficients
    std::vector< kw::sde_omega::info::expect::type > m_omega;

    //! Coefficients policy
    Coefficients coeff;
};

} // walker::

#endif // WrightFisher_h
