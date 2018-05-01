// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/CGCompFlow.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Compressible single-material flow using continuous Galerkin
  \details   This file implements the physics operators governing compressible
    single-material flow using continuous Galerkin discretization.
*/
// *****************************************************************************
#ifndef CGCompFlow_h
#define CGCompFlow_h

#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

#include "Macro.h"
#include "Exception.h"
#include "Vector.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace cg {

//! \brief CompFlow used polymorphically with tk::CGPDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/CompFlow/Physics.h
//!   - Problem - problem configuration, see PDE/CompFlow/Problems.h
//! \note The default physics is Euler, set in inciter::deck::check_compflow()
template< class Physics, class Problem >
class CompFlow {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! \brief Constructor
    explicit CompFlow( ncomp_t ) : m_offset( 0 ) {}

    //! Initalize the compressible flow equations, prepare for time integration
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    void initialize( const std::array< std::vector< tk::real >, 3 >& coord,
                     tk::Fields& unk,
                     tk::real t ) const
    {
      Assert( coord[0].size() == unk.nunk(), "Size mismatch" );
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
      // set initial and boundary conditions using problem policy
      for (ncomp_t i=0; i<coord[0].size(); ++i) {
        const auto s = Problem::solution( 0, x[i], y[i], z[i], t );
        unk(i,0,m_offset) = s[0]; // rho
        unk(i,1,m_offset) = s[1]; // rho * u
        unk(i,2,m_offset) = s[2]; // rho * v
        unk(i,3,m_offset) = s[3]; // rho * w
        unk(i,4,m_offset) = s[4]; // rho * e, e: total = kinetic + internal
      }
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
    void lhs( const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const std::pair< std::vector< std::size_t >,
                               std::vector< std::size_t > >& psup,
              tk::Fields& lhsd,
              tk::Fields& lhso ) const
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
        Assert( J > 0, "Element Jacobian non-positive" );

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
    //! \param[in] t Physical time
    //! \param[in] deltat Size of time step
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real t,
              tk::real deltat,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const tk::Fields& U,
              tk::Fields& Ue,
              tk::Fields& R ) const
    {
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      Assert( R.nunk() == coord[0].size() && R.nprop() == 5,
              "Number of unknowns and/or number of components in right-hand "
              "side vector incorrect" );
      Assert( U.nprop() == 5,
              "Number of components in solution vector must be 5" );
      Assert( R.nprop() == 5, "Number of components in rhs must be 5" );

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // ratio of specific heats
      auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0];

      // 1st stage: update element values from node values (gather-add)
      for (std::size_t e=0; e<inpoel.size()/4; ++e) {

        // access node IDs
        const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                               inpoel[e*4+2], inpoel[e*4+3] }};
        // compute element Jacobi determinant
        const std::array< tk::real, 3 >
          ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
          ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
          da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
        const auto J = tk::triple( ba, ca, da );        // J = 6V
        Assert( J > 0, "Element Jacobian non-positive" );

        // shape function derivatives, nnode*ndim [4][3]
        std::array< std::array< tk::real, 3 >, 4 > grad;
        grad[1] = tk::crossdiv( ca, da, J );
        grad[2] = tk::crossdiv( da, ba, J );
        grad[3] = tk::crossdiv( ba, ca, J );
        for (std::size_t i=0; i<3; ++i)
          grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

        // access solution at element nodes
        std::array< std::array< tk::real, 4 >, 5 > u;
        for (ncomp_t c=0; c<5; ++c) u[c] = U.extract( c, m_offset, N );
        // access solution at elements
        std::array< const tk::real*, 5 > ue;
        for (ncomp_t c=0; c<5; ++c) ue[c] = Ue.cptr( c, m_offset );

        // pressure
        std::array< tk::real, 4 > p;
        for (std::size_t a=0; a<4; ++a)
          p[a] = (g-1.0)*(u[4][a] - (u[1][a]*u[1][a] +
                                     u[2][a]*u[2][a] +
                                     u[3][a]*u[3][a])/2.0/u[0][a]);

        // sum nodal averages to element
        for (ncomp_t c=0; c<5; ++c) {
          Ue.var(ue[c],e) = 0.0;
          for (std::size_t a=0; a<4; ++a)
            Ue.var(ue[c],e) += u[c][a]/4.0;
        }

        // sum flux contributions to element
        tk::real d = deltat/2.0;
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t a=0; a<4; ++a) {
            // mass: advection
            Ue.var(ue[0],e) -= d * grad[a][j] * u[j+1][a];
            // momentum: advection
            for (std::size_t i=0; i<3; ++i)
              Ue.var(ue[i+1],e) -= d * grad[a][j] * u[j+1][a]*u[i+1][a]/u[0][a];
            // momentum: pressure
            Ue.var(ue[j+1],e) -= d * grad[a][j] * p[a];
            // energy: advection and pressure
            Ue.var(ue[4],e) -= d * grad[a][j] *
                              (u[4][a] + p[a]) * u[j+1][a]/u[0][a];
          }

        // add (optional) source to all equations
        std::array< std::array< tk::real, 5 >, 4 > s{{
          Problem::src( 0, x[N[0]], y[N[0]], z[N[0]], t ),
          Problem::src( 0, x[N[1]], y[N[1]], z[N[1]], t ),
          Problem::src( 0, x[N[2]], y[N[2]], z[N[2]], t ),
          Problem::src( 0, x[N[3]], y[N[3]], z[N[3]], t ) }};
        for (std::size_t c=0; c<5; ++c)
          for (std::size_t a=0; a<4; ++a)
            Ue.var(ue[c],e) += d/4.0 * s[a][c];

      }


      // zero right hand side for all components
      for (ncomp_t c=0; c<5; ++c) R.fill( c, m_offset, 0.0 );

      // 2nd stage: form rhs from element values (scatter-add)
      for (std::size_t e=0; e<inpoel.size()/4; ++e) {

        // access node IDs
        const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                               inpoel[e*4+2], inpoel[e*4+3] }};
        // compute element Jacobi determinant
        const std::array< tk::real, 3 >
          ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
          ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
          da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
        const auto J = tk::triple( ba, ca, da );        // J = 6V
        Assert( J > 0, "Element Jacobian non-positive" );

        // shape function derivatives, nnode*ndim [4][3]
        std::array< std::array< tk::real, 3 >, 4 > grad;
        grad[1] = tk::crossdiv( ca, da, J );
        grad[2] = tk::crossdiv( da, ba, J );
        grad[3] = tk::crossdiv( ba, ca, J );
        for (std::size_t i=0; i<3; ++i)
          grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

        // access solution at elements
        std::array< tk::real, 5 > ue;
        for (ncomp_t c=0; c<5; ++c) ue[c] = Ue( e, c, m_offset );
        // access pointer to right hand side at component and offset
        std::array< const tk::real*, 5 > r;
        for (ncomp_t c=0; c<5; ++c) r[c] = R.cptr( c, m_offset );

        // pressure
        auto p = (g-1.0)*(ue[4] -
                   (ue[1]*ue[1] + ue[2]*ue[2] + ue[3]*ue[3])/2.0/ue[0]);

        // scatter-add flux contributions to rhs at nodes
        tk::real d = deltat * J/6.0;
        for (std::size_t j=0; j<3; ++j)
          for (std::size_t a=0; a<4; ++a) {
            // mass: advection
            R.var(r[0],N[a]) += d * grad[a][j] * ue[j+1];
            // momentum: advection
            for (std::size_t i=0; i<3; ++i)
              R.var(r[i+1],N[a]) += d * grad[a][j] * ue[j+1]*ue[i+1]/ue[0];
            // momentum: pressure
            R.var(r[j+1],N[a]) += d * grad[a][j] * p;
            // energy: advection and pressure
            R.var(r[4],N[a]) += d * grad[a][j] * (ue[4] + p) * ue[j+1]/ue[0];
          }

        // add (optional) source to all equations
        auto xc = (x[N[0]] + x[N[1]] + x[N[2]] + x[N[3]]) / 4.0;
        auto yc = (y[N[0]] + y[N[1]] + y[N[2]] + y[N[3]]) / 4.0;
        auto zc = (z[N[0]] + z[N[1]] + z[N[2]] + z[N[3]]) / 4.0;
        auto s = Problem::src( 0, xc, yc, zc, t+deltat/2 );
        for (std::size_t c=0; c<5; ++c)
          for (std::size_t a=0; a<4; ++a)
            R.var(r[c],N[a]) += d/4.0 * s[c];

      }
//         // add viscous stress contribution to momentum and energy rhs
//         Physics::viscousRhs( deltat, J, N, grad, u, r, R );
//         // add heat conduction contribution to energy rhs
//         Physics::conductRhs( deltat, J, N, grad, u, r, R );
    }

    //! Compute the minimum time step size
    //! \param[in] U Solution vector at recent time step
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \return Minimum time step size
    tk::real dt( const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< std::size_t >& inpoel,
                 const tk::Fields& U ) const
    {
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
      // ratio of specific heats
      auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0];
      // compute the minimum dt across all elements we own
      tk::real mindt = std::numeric_limits< tk::real >::max();
      for (std::size_t e=0; e<inpoel.size()/4; ++e) {
        const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                               inpoel[e*4+2], inpoel[e*4+3] }};
        // compute cubic root of element volume as the characteristic length
        const std::array< tk::real, 3 >
          ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
          ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
          da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
        const auto L = std::cbrt( tk::triple( ba, ca, da ) / 6.0 );
        // access solution at element nodes at recent time step
        std::array< std::array< tk::real, 4 >, 5 > u;
        for (ncomp_t c=0; c<5; ++c) u[c] = U.extract( c, m_offset, N );
        // compute the maximum length of the characteristic velocity (fluid
        // velocity + sound velocity) across the four element nodes
        tk::real maxvel = 0.0;
        for (std::size_t j=0; j<4; ++j) {
          auto& r  = u[0][j];    // rho
          auto& ru = u[1][j];    // rho * u
          auto& rv = u[2][j];    // rho * v
          auto& rw = u[3][j];    // rho * w
          auto& re = u[4][j];    // rho * e
          auto p = (g-1.0)*(re - (ru*ru + rv*rv + rw*rw)/2.0/r); // pressure
          if (p < 0) p = 0.0;
          auto c = std::sqrt(g*p/r);     // sound speed
          auto v = std::sqrt((ru*ru + rv*rv + rw*rw)/r/r) + c; // char. velocity
          if (v > maxvel) maxvel = v;
        }
        // compute element dt for the Euler equations
        auto euler_dt = L / maxvel;
        // compute element dt based on the viscous force
        auto viscous_dt = Physics::viscous_dt( L, u );
        // compute element dt based on thermal diffusion
        auto conduct_dt = Physics::conduct_dt( L, u );
        // compute minimum element dt
        auto elemdt = std::min( euler_dt, std::min( viscous_dt, conduct_dt ) );
        // find minimum dt across all elements
        if (elemdt < mindt) mindt = elemdt;
      }
      return mindt;
    }

    //! Extract the velocity field at cell nodes
    //! \param[in] U Solution vector at recent time step
    //! \param[in] N Element node indices    
    //! \return Array of the four values of the velocity field
    std::array< std::array< tk::real, 4 >, 3 >
    velocity( const tk::Fields& U,
              const std::array< std::vector< tk::real >, 3 >&,
              const std::array< std::size_t, 4 >& N ) const
    {
      std::array< std::array< tk::real, 4 >, 3 > v;
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

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    void side( std::unordered_set< int >& conf ) const
    { Problem::side( conf ); }

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
    //! \param[in] t Physical time
    //! \param[in] deltat Time step size
    //! \param[in] ss Pair of side set ID and node IDs on the side set
    //! \param[in] coord Mesh node coordinates
    //! \return Vector of pairs of bool and boundary condition value associated
    //!   to mesh node IDs at which Dirichlet boundary conditions are set. Note
    //!   that instead of the actual boundary condition value, we return the
    //!   increment between t+dt and t, since that is what the solution requires
    //!   as we solve for the soution increments and not the solution itself.
    std::unordered_map< std::size_t, std::vector< std::pair<bool,tk::real> > >
    dirbc( tk::real t,
           tk::real deltat,
           const std::pair< const int, std::vector< std::size_t > >& ss,
           const std::array< std::vector< tk::real >, 3 >& coord ) const
    {
      using tag::param; using tag::compflow; using tag::bcdir;
      using NodeBC = std::vector< std::pair< bool, tk::real > >;
      std::unordered_map< std::size_t, NodeBC > bc;
      const auto& ubc = g_inputdeck.get< param, compflow, bcdir >();
      if (!ubc.empty()) {
        Assert( ubc.size() > 0, "Indexing out of Dirichlet BC eq-vector" );
        const auto& x = coord[0];
        const auto& y = coord[1];
        const auto& z = coord[2];
        for (const auto& b : ubc[0])
          if (std::stoi(b) == ss.first)
            for (auto n : ss.second) {
              Assert( x.size() > n, "Indexing out of coordinate array" );
              auto s = Problem::solinc( 0, x[n], y[n], z[n], t, deltat );
              bc[n] = {{ {true,s[0]}, {true,s[1]}, {true,s[2]}, {true,s[3]},
                         {true,s[4]} }};
            }
      }
      return bc;
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    std::vector< std::string > fieldNames() const
    { return Problem::fieldNames(); }

    //! Return field output going to file
    //! \param[in] t Physical time
    //! \param[in] V Total mesh volume
    //! \param[in] coord Mesh node coordinates
    //! \param[in] v Nodal mesh volumes
    //! \param[in,out] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    std::vector< std::vector< tk::real > >
    fieldOutput( tk::real t,
                 tk::real V,
                 const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< tk::real >& v,
                 tk::Fields& U ) const
    { return Problem::fieldOutput( 0, m_offset, t, V, v, coord, U ); }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const
    { return Problem::names(); }

  private:
    const ncomp_t m_offset;             //!< Offset PDE operates from
};

} // cg::

} // inciter::

#endif // CGCompFlow_h
