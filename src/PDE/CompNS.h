// *****************************************************************************
/*!
  \file      src/PDE/CompNS
  \author    J. Bakosi
  \date      Mon 02 May 2016 04:58:30 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Navier-Stokes equations describing compressible flow
  \details   This file implements the time integration of the Navier-Stokes
    equations governing compressible fluid flow.
*/
// *****************************************************************************
#ifndef CompNS_h
#define CompNS_h

#include "Macro.h"
#include "CompNSProblem.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! \brief CompNS used polymorphically with tk::PDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Problem - problem configuration, see PDE/CompNSProblem.h
template< class Problem >
class CompNS {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! \brief Constructor
    //! \author J. Bakosi
    explicit CompNS( ncomp_t c ) :
      m_ncomp( 5 ),
      m_offset( 0 )
    {
      IGNORE(c);
    }

    //! Initalize the Navier-Stokes equations, prepare for time integration
    //! \param[in,out] unk Array of unknowns
    //! \author J. Bakosi
    void initialize( const std::array< std::vector< tk::real >, 3 >& coord,
                     tk::MeshNodes& unk,
                     tk::real t )
    {
      IGNORE(coord);
      IGNORE(unk);
      IGNORE(t);
      //! Set initial conditions using problem configuration policy
      //Problem::template init< tag::euler >( g_inputdeck, unk, m_offset );
    }

    //! Compute the left hand side sparse matrix
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] psup Linked lists storing IDs of points surrounding points
    //! \param[in,out] lhsd Diagonal of the sparse matrix storing nonzeros
    //! \param[in,out] lhso Off-diagonal of the sparse matrix storing nonzeros
    //! \details Sparse matrix storing the nonzero matrix values at rows and
    //!   columns given by psup. The format is similar to compressed row
    //!   storage, but the diagonal and off-diagonal data are stored in separate
    //!   vectors. For the off-diagonal data the local row and column indices,
    //!   at which values are nonzero, are stored by psup (psup1 and psup2,
    //!   where psup2 holds the indices at which psup1 holds the point ids
    //!   surrounding points, see also tk::genPsup()). Note that the number of
    //!   mesh points (our chunk) npoin = psup.second.size()-1.
    //! \author J. Bakosi
    void lhs( const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const std::pair< std::vector< std::size_t >,
                               std::vector< std::size_t > >& psup,
              tk::MeshNodes& lhsd,
              tk::MeshNodes& lhso )
    {
      Assert( psup.second.size()-1 == coord[0].size(),
              "Number of mesh points and number of global IDs unequal" );
      Assert( lhsd.nunk() == psup.second.size()-1, "Number of unknowns in "
              "diagonal sparse matrix storage incorrect" );
      Assert( lhso.nunk() == psup.first.size(), "Number of unknowns in "
              "off-diagonal sparse matrix storage incorrect" );

      // Lambda to compute the sparse matrix vector index for row and column
      // indices. Used only for off-diagonal entries.
      auto spidx = [ &psup ]( std::size_t r, std::size_t c ) -> std::size_t {
        Assert( r != c, "Only for computing the off-diagonal indices" );
        for (auto i=psup.second[r]+1; i<=psup.second[r+1]; ++i)
          if (c == psup.first[i]) return i;
        Throw( "Cannot find row, column: " + std::to_string(r) + ',' +
               std::to_string(c) + " in sparse matrix" );
      };

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // Zero matrix for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) {
        lhsd.fill( c, m_offset, 0.0 );
        lhso.fill( c, m_offset, 0.0 );
      }

      for (std::size_t e=0; e<inpoel.size()/4; ++e) {
        const auto A = inpoel[e*4+0];
        const auto B = inpoel[e*4+1];
        const auto C = inpoel[e*4+2];
        const auto D = inpoel[e*4+3];
        std::array< tk::real, 3 > ba{{ x[B]-x[A], y[B]-y[A], z[B]-z[A] }},
                                  ca{{ x[C]-x[A], y[C]-y[A], z[C]-z[A] }},
                                  da{{ x[D]-x[A], y[D]-y[A], z[D]-z[A] }};
        const auto J = tk::triple( ba, ca, da ) / 120.0;

        for (ncomp_t c=0; c<m_ncomp; ++c) {
          const auto r = lhsd.cptr( c, m_offset );
          lhsd.var( r, A ) += 2.0 * J;
          lhsd.var( r, B ) += 2.0 * J;
          lhsd.var( r, C ) += 2.0 * J;
          lhsd.var( r, D ) += 2.0 * J;

          const auto s = lhso.cptr( c, m_offset );
          lhso.var( s, spidx(A,B) ) += J;
          lhso.var( s, spidx(A,C) ) += J;
          lhso.var( s, spidx(A,D) ) += J;

          lhso.var( s, spidx(B,A) ) += J;
          lhso.var( s, spidx(B,C) ) += J;
          lhso.var( s, spidx(B,D) ) += J;

          lhso.var( s, spidx(C,A) ) += J;
          lhso.var( s, spidx(C,B) ) += J;
          lhso.var( s, spidx(C,D) ) += J;

          lhso.var( s, spidx(D,A) ) += J;
          lhso.var( s, spidx(D,B) ) += J;
          lhso.var( s, spidx(D,C) ) += J;
        }
      }
    }

    //! Compute right hand side
    //! \param[in] mult Multiplier differentiating the different stages in
    //!    multi-stage time stepping
    //! \param[in] dt Size of time step
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] U Solution vector at recent time step stage
    //! \param[in] Un Solution vector at previous time step
    //! \param[in,out] R Right-hand side vector computed
    //! \author J. Bakosi
    void rhs( tk::real mult,
              tk::real dt,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const tk::MeshNodes& U,
              const tk::MeshNodes& Un,
              tk::MeshNodes& R )
    {
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      Assert( Un.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at previous time step incorrect" );
      Assert( R.nunk() == coord[0].size(), "Number of unknowns in right-hand "
              "side vector incorrect" );

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // zero right hand side for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) R.fill( c, m_offset, 0.0 );

      for (std::size_t e=0; e<inpoel.size()/4; ++e) {
        const auto A = inpoel[e*4+0];
        const auto B = inpoel[e*4+1];
        const auto C = inpoel[e*4+2];
        const auto D = inpoel[e*4+3];

        // compute element Jacobi determinant
        const std::array< tk::real, 3 > ba{{ x[B]-x[A], y[B]-y[A], z[B]-z[A] }},
                                        ca{{ x[C]-x[A], y[C]-y[A], z[C]-z[A] }},
                                        da{{ x[D]-x[A], y[D]-y[A], z[D]-z[A] }};
        const auto J = tk::triple( ba, ca, da );

        // construct tetrahedron element-level matrices

        // consistent mass
        std::array< std::array< tk::real, 4 >, 4 > mass;  // nnode*nnode [4][4]
        // diagonal
        mass[0][0] = mass[1][1] = mass[2][2] = mass[3][3] = J/60.0;
        // off-diagonal
        mass[0][1] = mass[0][2] = mass[0][3] =
        mass[1][0] = mass[1][2] = mass[1][3] =
        mass[2][0] = mass[2][1] = mass[2][3] =
        mass[3][0] = mass[3][1] = mass[3][2] = J/120.0;

        // shape function derivatives
        std::array< std::array< tk::real, 3 >, 4 > grad;  // nnode*ndim [4][3]
        grad[1] = tk::crossdiv( ca, da, J );
        grad[2] = tk::crossdiv( da, ba, J );
        grad[3] = tk::crossdiv( ba, ca, J );
        for (std::size_t i=0; i<3; ++i)
          grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

        // access solution at element nodes at time n
        std::vector< std::array< tk::real, 4 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c)
          u[c] = Un.extract( c, m_offset, A, B, C, D );
        // access solution at element nodes at recent time step stage
        std::vector< std::array< tk::real, 4 > > s( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c)
          s[c] = U.extract( c, m_offset, A, B, C, D );
        // access pointer to right hand side at component and offset
        std::vector< const tk::real* > r( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) r[c] = R.cptr( c, m_offset );

        tk::real gamma = 1.4;
        tk::real mu = 0.1;
	tk::real cv = 1.005;
	// thermal conductivity
	tk::real kc = 0.029;

        // compute pressure
        std::array< tk::real, 4 > p;
        for (std::size_t i=0; i<4; ++i)
          p[i] = (gamma-1.0)*(s[4][i] - (s[1][i]*s[1][i] + s[2][i]*s[2][i] + s[3][i]*s[3][i])/2.0/s[0][i]);

        // compute temperature
        std::array< tk::real, 4 > T;
        for (std::size_t i=0; i<4; ++i)
          T[i] = cv*(s[4][i]/s[0][i] - (s[1][i]*s[1][i] + s[2][i]*s[2][i] + s[3][i]*s[3][i])/2.0/s[0][i]);
 
        // add mass contribution to right hand side for all equations
        for (ncomp_t c=0; c<m_ncomp; ++c)
          for (std::size_t j=0; j<4; ++j) {
            R.var(r[c],A) += mass[0][j] * u[c][j];
            R.var(r[c],B) += mass[1][j] * u[c][j];
            R.var(r[c],C) += mass[2][j] * u[c][j];
            R.var(r[c],D) += mass[3][j] * u[c][j];
          }

        // add advection contribution for conservation of mass to right hand side
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t k=0; k<4; ++k) {
            R.var(r[0],A) -= mult * dt * J/24.0 * grad[k][j] * s[j+1][k];
            R.var(r[0],B) -= mult * dt * J/24.0 * grad[k][j] * s[j+1][k];
            R.var(r[0],C) -= mult * dt * J/24.0 * grad[k][j] * s[j+1][k];
            R.var(r[0],D) -= mult * dt * J/24.0 * grad[k][j] * s[j+1][k];
          }

        // add advection contribution for conservation of momentum to right hand side
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t k=0; k<4; ++k) {
              R.var(r[i+1],A) -= mult * dt * J/24.0 * grad[k][j] * s[i+1][k]*s[j+1][k]/s[0][k];
              R.var(r[i+1],B) -= mult * dt * J/24.0 * grad[k][j] * s[i+1][k]*s[j+1][k]/s[0][k];
              R.var(r[i+1],C) -= mult * dt * J/24.0 * grad[k][j] * s[i+1][k]*s[j+1][k]/s[0][k];
              R.var(r[i+1],D) -= mult * dt * J/24.0 * grad[k][j] * s[i+1][k]*s[j+1][k]/s[0][k];
            }

        // add pressure gradient contribution for conservation of momentum to right hand side
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t k=0; k<4; ++k) {
              R.var(r[i+1],A) -= mult * dt * J/24.0 * grad[k][j] * p[k];
              R.var(r[i+1],B) -= mult * dt * J/24.0 * grad[k][j] * p[k];
              R.var(r[i+1],C) -= mult * dt * J/24.0 * grad[k][j] * p[k];
              R.var(r[i+1],D) -= mult * dt * J/24.0 * grad[k][j] * p[k];
            }

        // add deviatoric viscous stress contribution for conservation of momentum to right hand side
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t k=0; k<4; ++k) {
              R.var(r[i+1],A) -= mult * dt * J/6.0 *
                                 grad[0][j]*(grad[k][j]*s[i+1][k] + grad[k][i]*s[j+1][k])/s[0][k]*mu;
              R.var(r[i+1],B) -= mult * dt * J/6.0 *
                                 grad[1][j]*(grad[k][j]*s[i+1][k] + grad[k][i]*s[j+1][k])/s[0][k]*mu;
              R.var(r[i+1],C) -= mult * dt * J/6.0 *
                                 grad[2][j]*(grad[k][j]*s[i+1][k] + grad[k][i]*s[j+1][k])/s[0][k]*mu;
              R.var(r[i+1],D) -= mult * dt * J/6.0 *
                                 grad[3][j]*(grad[k][j]*s[i+1][k] + grad[k][i]*s[j+1][k])/s[0][k]*mu;
            }

        // add isotropic viscous stress contribution for conservation of momentum to right hand side
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t k=0; k<4; ++k) {
              R.var(r[i+1],A) += mult * dt * J/6.0 *
                                 grad[0][i]*grad[k][j]*s[j+1][k]/s[0][k]*2.0/3.0*mu;
              R.var(r[i+1],B) += mult * dt * J/6.0 *
                                 grad[1][i]*grad[k][j]*s[j+1][k]/s[0][k]*2.0/3.0*mu;
              R.var(r[i+1],C) += mult * dt * J/6.0 *
                                 grad[2][i]*grad[k][j]*s[j+1][k]/s[0][k]*2.0/3.0*mu;
              R.var(r[i+1],D) += mult * dt * J/6.0 *
                                 grad[3][i]*grad[k][j]*s[j+1][k]/s[0][k]*2.0/3.0*mu;
            }

        // add advection and pressure gradient contribution for conservation of energy to right hand side
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t k=0; k<4; ++k) {
            R.var(r[4],A) -= mult * dt * J/24.0 * grad[k][j] * (s[4][k]+p[k])*s[j+1][k]/s[0][k];
            R.var(r[4],B) -= mult * dt * J/24.0 * grad[k][j] * (s[4][k]+p[k])*s[j+1][k]/s[0][k];
            R.var(r[4],C) -= mult * dt * J/24.0 * grad[k][j] * (s[4][k]+p[k])*s[j+1][k]/s[0][k];
            R.var(r[4],D) -= mult * dt * J/24.0 * grad[k][j] * (s[4][k]+p[k])*s[j+1][k]/s[0][k];
          }

        // add deviatoric viscous stress contribution for conservation of energy to right hand side
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t k=0; k<4; ++k) {
              R.var(r[4],A) -= mult * dt * J/24.0 * s[i+1][k]/s[0][k] *
                                 grad[0][j]*(grad[k][j]*s[i+1][k] + grad[k][i]*s[j+1][k])/s[0][k]*mu;
              R.var(r[4],B) -= mult * dt * J/24.0 * s[i+1][k]/s[0][k] *
                                 grad[1][j]*(grad[k][j]*s[i+1][k] + grad[k][i]*s[j+1][k])/s[0][k]*mu;
              R.var(r[4],C) -= mult * dt * J/24.0 * s[i+1][k]/s[0][k] *
                                 grad[2][j]*(grad[k][j]*s[i+1][k] + grad[k][i]*s[j+1][k])/s[0][k]*mu;
              R.var(r[4],D) -= mult * dt * J/24.0 * s[i+1][k]/s[0][k] *
                                 grad[3][j]*(grad[k][j]*s[i+1][k] + grad[k][i]*s[j+1][k])/s[0][k]*mu;
            }

        // add isotropic viscous stress contribution for conservation of energy to right hand side
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t k=0; k<4; ++k) {
              R.var(r[4],A) += mult * dt * J/24.0 * s[i+1][k]/s[0][k] *
                                 grad[0][i]*grad[k][j]*s[j+1][k]/s[0][k]*2.0/3.0*mu;
              R.var(r[4],B) += mult * dt * J/24.0 * s[i+1][k]/s[0][k] *
                                 grad[1][i]*grad[k][j]*s[j+1][k]/s[0][k]*2.0/3.0*mu;
              R.var(r[4],C) += mult * dt * J/24.0 * s[i+1][k]/s[0][k] *
                                 grad[2][i]*grad[k][j]*s[j+1][k]/s[0][k]*2.0/3.0*mu;
              R.var(r[4],D) += mult * dt * J/24.0 * s[i+1][k]/s[0][k] *
                                 grad[3][i]*grad[k][j]*s[j+1][k]/s[0][k]*2.0/3.0*mu;
            }

        // add conduction contribution for conservation of energy to right hand side
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t k=0; k<4; ++k) {
            R.var(r[4],A) += mult * dt * J/24.0 * grad[k][j] * T[k] * kc;
            R.var(r[4],B) += mult * dt * J/24.0 * grad[k][j] * T[k] * kc;
            R.var(r[4],C) += mult * dt * J/24.0 * grad[k][j] * T[k] * kc;
            R.var(r[4],D) += mult * dt * J/24.0 * grad[k][j] * T[k] * kc;
          }

      }
    }

    //! \brief Advance unknowns according to the Euler equations
    //! \param[in,out] unk Array of unknowns
    //! \param[in] dt Time step size
    //! \param[in] t Physical time
    //! \author J. Bakosi
    void advance( tk::MeshNodes& unk, tk::real dt, tk::real t ) {
      IGNORE(unk);
      IGNORE(dt);
      IGNORE(t);
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    std::vector< std::string > names() {
      std::vector< std::string > n( m_ncomp );
      // ...
      return n;
    }

    //! Return field output going to file
    //! \param[in] t Physical time
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] U Solution vector at recent time step stage
    //! \return Vector of vectors to be output to file
    //! \details Note that U is overwritten
    std::vector< std::vector< tk::real > >
    output( tk::real t,
            const std::array< std::vector< tk::real >, 3 >& coord,
            tk::MeshNodes& U )
    {
      IGNORE(t);
      IGNORE(coord);
      IGNORE(U);
      std::vector< std::vector< tk::real > > out;
      // ...
      return out;
   }

  private:
    const ncomp_t m_ncomp;              //!< Number of components
    const ncomp_t m_offset;             //!< Offset PDE operates from
};

} // inciter::

#endif // CompNS_h
