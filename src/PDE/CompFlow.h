// *****************************************************************************
/*!
  \file      src/PDE/CompFlow.h
  \author    J. Bakosi
  \date      Thu 01 Sep 2016 08:20:55 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Governing equations describing compressible single-phase flow
  \details   This file implements the time integration of the equations
     governing compressible fluid flow.
*/
// *****************************************************************************
#ifndef CompFlow_h
#define CompFlow_h

#include <algorithm>

#include "Macro.h"
#include "Keywords.h"
#include "CompFlowPhysics.h"
#include "CompFlowProblem.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! \brief CompFlow used polymorphically with tk::PDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/CompFlowPhysics.h
//!   - Problem - problem configuration, see PDE/CompFlowProblem.h
//! \note The default physics is Euler, set in inciter::deck::check_compflow()
template< class Physics, class Problem >
class CompFlow {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! \brief Constructor
    //! \author J. Bakosi
    explicit CompFlow( ncomp_t ) : m_offset( 0 ) {}

    //! Initalize the compressible flow equations, prepare for time integration
    //! \param[in,out] unk Array of unknowns
    //! \param[in] coord Mesh node coordinates
    //! \author J. Bakosi
    void initialize( const std::array< std::vector< tk::real >, 3 >& coord,
                     tk::MeshNodes& unk,
                     tk::real ) const
    {
      //! Set initial conditions using problem configuration policy
      Problem::init( g_inputdeck, coord, unk, 0, m_offset );
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
              tk::MeshNodes& lhso ) const
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
      for (ncomp_t c=0; c<5; ++c) {
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

        for (ncomp_t c=0; c<5; ++c) {
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
    //! \param[in,out] R Right-hand side vector computed
    //! \author J. Bakosi
    void rhs( tk::real mult,
              tk::real dt,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const tk::MeshNodes& U,
              tk::MeshNodes& R ) const
    {
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      Assert( R.nunk() == coord[0].size(), "Number of unknowns in right-hand "
              "side vector incorrect" );

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // zero right hand side for all components
      for (ncomp_t c=0; c<5; ++c) R.fill( c, m_offset, 0.0 );

      for (std::size_t e=0; e<inpoel.size()/4; ++e) {
        const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                               inpoel[e*4+2], inpoel[e*4+3] }};
        // compute element Jacobi determinant
        const std::array< tk::real, 3 >
          ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
          ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
          da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
        const auto J = tk::triple( ba, ca, da );

        // construct tetrahedron element-level matrices

        // consistent mass, nnode*nnode [4][4]
        std::array< std::array< tk::real, 4 >, 4 > mass;
        mass[0][0] = mass[1][1] = mass[2][2] = mass[3][3] = J/60.0;  // diagonal
        mass[0][1] = mass[0][2] = mass[0][3] =                   // off-diagonal
        mass[1][0] = mass[1][2] = mass[1][3] =
        mass[2][0] = mass[2][1] = mass[2][3] =
        mass[3][0] = mass[3][1] = mass[3][2] = J/120.0;

        // shape function derivatives, nnode*ndim [4][3]
        std::array< std::array< tk::real, 3 >, 4 > grad;
        grad[1] = tk::crossdiv( ca, da, J );
        grad[2] = tk::crossdiv( da, ba, J );
        grad[3] = tk::crossdiv( ba, ca, J );
        for (std::size_t i=0; i<3; ++i)
          grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

        // access solution at element nodes at recent time step stage
        std::array< std::array< tk::real, 4 >, 5 > s;
        for (ncomp_t c=0; c<5; ++c) s[c] = U.extract( c, m_offset, N );
        // access pointer to right hand side at component and offset
        std::array< const tk::real*, 5 > r;
        for (ncomp_t c=0; c<5; ++c) r[c] = R.cptr( c, m_offset );

        // ratio of specific heats
        tk::real g =
          g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0];

        // add source to rhs for all equations
        Problem::sourceRhs( coord, 0, mult, dt, J, N, grad, mass, r, s, R,
                            const_cast<tk::MeshNodes&>(U) );

        // compute pressure
        std::array< tk::real, 4 > p;
        for (std::size_t i=0; i<4; ++i)
          p[i] = (g-1.0)*(s[4][i] - (s[1][i]*s[1][i] +
                                     s[2][i]*s[2][i] +
                                     s[3][i]*s[3][i])/2.0/s[0][i]);

        tk::real c = mult * dt * J/24.0;
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<4; ++j)
            for (std::size_t k=0; k<4; ++k) {
              // advection contribution to mass rhs
              R.var(r[0],N[j]) -= c * grad[k][i] * s[i+1][k];

              // advection contribution to momentum rhs
              for (std::size_t l=0; l<3; ++l)
                R.var(r[i+1],N[j]) -= c * grad[k][l] *
                                      s[i+1][k] * s[l+1][k]/s[0][k];

              // pressure gradient contribution to momentum rhs
              R.var(r[i+1],N[j]) -= c * grad[k][i] * p[k];

              // advection and pressure gradient contribution to energy rhs
              R.var(r[4],N[j]) -= c * grad[k][i] *
                                  (s[4][k] + p[k]) * s[i+1][k]/s[0][k];
            }

        // add viscous stress contribution to momentum and energy rhs
        Physics::viscousRhs( mult, dt, J, N, grad, s, r, R );

        // add heat conduction contribution to energy rhs
        Physics::conductRhs( mult, dt, J, N, grad, s, r, R );
      }
    }

    //! Extract the velocity field at cell nodes
    //! \param[in] U Solution vector at recent time step stage
    //! \param[in] N Element node indices    
    //! \return Array of the four values of the three velocity coordinates
    std::vector< std::array< tk::real, 4 > >
    velocity( const tk::MeshNodes& U,
              const std::array< std::vector< tk::real >, 3 >&,
              const std::array< std::size_t, 4 >& N ) const
    {
      std::vector< std::array< tk::real, 4 > > v( 3 );
      v[0] = U.extract( 1, m_offset, N );
      v[1] = U.extract( 2, m_offset, N );
      v[2] = U.extract( 3, m_offset, N );
      auto r = U.extract( 0, m_offset, N );
      std::transform( r.begin(), r.end(), v[0].begin(), v[0].begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      std::transform( r.begin(), r.end(), v[1].begin(), v[1].begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      std::transform( r.begin(), r.end(), v[2].begin(), v[2].begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      return v;
    }

    //! \brief Query if a Dirichlet boundary condition has been set on a side
    //!   set for any component in the PDE system
    //! \param[in] sideset Side set ID
    //! \return True if a Dirichlet boundary condition has been set on the side
    //!   set for any component in the PDE system.
    bool anydirbc( int sideset ) const {
      const auto& bc =
        g_inputdeck.get< tag::param, tag::compflow, tag::bc_dirichlet >();
      for (const auto& s : bc)
        if (static_cast<int>(std::round(s[0])) == sideset)
          return true;
      return false;
    }

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
    //! \param[in] sideset Side set ID
    //! \return Vector of pairs of bool and BC value for all components
    std::vector< std::pair< bool, tk::real > > dirbc( int sideset ) const {
      const auto& bc =
        g_inputdeck.get< tag::param, tag::compflow, tag::bc_dirichlet >();
      std::vector< std::pair< bool, tk::real > > b( 5, { false, 0.0 } );
      for (const auto& s : bc) {
        Assert( s.size() == 3, "Side set vector size incorrect" );
        if (static_cast<int>(std::round(s[0])) == sideset) {
          b[ static_cast<std::size_t>(std::round(s[1]))-1 ] = { true, s[2] };
        }
      }
      return b;
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    std::vector< std::string > names() const { return Problem::names(); }

    //! Return field output going to file
    //! \param[in] U Solution vector at recent time step stage
    //! \return Vector of vectors to be output to file
    std::vector< std::vector< tk::real > >
    output( tk::real t,
            const std::array< std::vector< tk::real >, 3 >& coord,
            const tk::MeshNodes& U ) const
    { return Problem::output( 0, m_offset, t, coord, U ); }

  private:
    const ncomp_t m_offset;             //!< Offset PDE operates from
};

} // inciter::

#endif // CompFlow_h
