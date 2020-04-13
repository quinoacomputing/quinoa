// *****************************************************************************
/*!
  \file      src/PDE/Transport/CGTransport.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Scalar transport using continous Galerkin discretization
  \details   This file implements the physics operators governing transported
     scalars using continuous Galerkin discretization.
*/
// *****************************************************************************
#ifndef CGTransport_h
#define CGTransport_h

#include <vector>
#include <array>
#include <limits>
#include <cmath>
#include <unordered_set>
#include <unordered_map>

#include "Exception.hpp"
#include "Vector.hpp"
#include "DerivedData.hpp"
#include "Around.hpp"
#include "Reconstruction.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "CGPDE.hpp"
#include "History.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace cg {

//! \brief Transport equation used polymorphically with tk::CGPDE
//! \details The template argument(s) specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/Transport/Physics/CG.h
//!   - Problem - problem configuration, see PDE/Transport/Problem.h
//! \note The default physics is CGAdvection, set in
//!    inciter::deck::check_transport()
template< class Physics, class Problem >
class Transport {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor
    //! \param[in] c Equation system index (among multiple systems configured)
    explicit Transport( ncomp_t c ) :
      m_physics( Physics() ),
      m_problem( Problem() ),
      m_system( c ),
      m_ncomp(
        g_inputdeck.get< tag::component >().get< tag::transport >().at(c) ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< tag::transport >(c) )
    {
      m_problem.errchk( m_system, m_ncomp );
    }

    //! Initalize the transport equations using problem policy
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    void initialize( const std::array< std::vector< tk::real >, 3 >& coord,
                     tk::Fields& unk,
                     tk::real t,
                     std::vector< std::size_t >& ) const
    {
      Assert( coord[0].size() == unk.nunk(), "Size mismatch" );
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
      for (ncomp_t i=0; i<x.size(); ++i) {
        int inbox = 0;
        const auto s =
          Problem::solution( m_system, m_ncomp, x[i], y[i], z[i], t, inbox );
        for (ncomp_t c=0; c<m_ncomp; ++c)
          unk( i, c, m_offset ) = s[c];
      }
    }

    //! Set initial condition in user-defined box IC nodes (no-op for Transport)
    void box( tk::real, const std::vector< std::size_t >&, tk::Fields& )
    const {}

    //! Return analytic solution (if defined by Problem) at xi, yi, zi, t
    //! \param[in] xi X-coordinate
    //! \param[in] yi Y-coordinate
    //! \param[in] zi Z-coordinate
    //! \param[in] t Physical time
    //! \return Vector of analytic solution at given location and time
    std::vector< tk::real >
    analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    {
      int inbox = 0;
      auto s = Problem::solution( m_system, m_ncomp, xi, yi, zi, t, inbox );
      return std::vector< tk::real >( begin(s), end(s) );
    }

    //! Compute nodal gradients of primitive variables for ALECG
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] bndel List of elements contributing to chare-boundary nodes
    //! \param[in] gid Local->global node id map
    //! \param[in] bid Local chare-boundary node ids (value) associated to
    //!    global node ids (key)
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] G Nodal gradients of primitive variables
    void grad( const std::array< std::vector< tk::real >, 3 >& coord,
               const std::vector< std::size_t >& inpoel,
               const std::vector< std::size_t >& bndel,
               const std::vector< std::size_t >& gid,
               const std::unordered_map< std::size_t, std::size_t >& bid,
               const tk::Fields& U,
               tk::Fields& G ) const
    {
      chbgrad( m_ncomp, m_offset, coord, inpoel, bndel, gid, bid, {},
               U, egrad, G );
    }

    //! Compute right hand side for ALECG
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] triinpoel Boundary triangle face connecitivity
    //! \param[in] bid Local chare-boundary node ids (value) associated to
    //!    global node ids (key)
    //! \param[in] lid Global->local node ids
    //! \param[in] dfn Dual-face normals
    //! \param[in] psup Points surrounding points
    //! \param[in] bnorm Face normals in boundary points
    //! \param[in] vol Nodal volumes
    //! \param[in] G Nodal gradients in chare-boundary nodes
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    void rhs(
      tk::real,
      const std::array< std::vector< tk::real >, 3 >&  coord,
      const std::vector< std::size_t >& inpoel,
      const std::vector< std::size_t >& triinpoel,
      const std::unordered_map< std::size_t, std::size_t >& bid,
      const std::unordered_map< std::size_t, std::size_t >& lid,
      const std::vector< tk::real >& dfn,
      const std::pair< std::vector< std::size_t >,
                       std::vector< std::size_t > >& psup,
      const std::unordered_map< std::size_t, std::array< tk::real, 4 > >& bnorm,
      const std::vector< tk::real >& vol,
      const tk::Fields& G,
      const tk::Fields& U,
      tk::Fields& R) const
    {
      Assert( G.nprop() == m_ncomp*3,
              "Number of components in gradient vector incorrect" );
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      Assert( R.nunk() == coord[0].size(),
              "Number of unknowns and/or number of components in right-hand "
              "side vector incorrect" );

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // zero right hand side for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) R.fill( c, m_offset, 0.0 );

      // access pointer to right hand side at component and offset
      std::vector< const tk::real* > r( m_ncomp );
      for (ncomp_t c=0; c<m_ncomp; ++c) r[c] = R.cptr( c, m_offset );

      // compute/assemble gradients in points
      auto Grad = nodegrad( m_ncomp, m_offset, coord, inpoel, lid, bid,
                            vol, {}, U, G, egrad );

      // compute derived data structures
      auto esued = tk::genEsued( inpoel, 4, tk::genEsup( inpoel, 4 ) );

      // domain-edge integral
      for (std::size_t p=0,k=0; p<U.nunk(); ++p) {  // for each point p
        for (auto q : tk::Around(psup,p)) {
          // access dual-face normals for edge p-q
          std::array< tk::real, 3 > n{ dfn[k+0], dfn[k+1], dfn[k+2] };
          k += 6;
          // compute primitive variables at edge-end points (for Transport,
          // these are the same as the conserved variables)
          std::vector< tk::real > uL( m_ncomp, 0.0 );
          std::vector< tk::real > uR( m_ncomp, 0.0 );
          for (std::size_t c=0; c<m_ncomp; ++c) {
            uL[c] = U(p,c,m_offset);
            uR[c] = U(q,c,m_offset);
          }
          // compute MUSCL reconstruction in edge-end points
          tk::muscl( {p,q}, coord, Grad, uL, uR );
          // evaluate prescribed velocity
          auto v =
            Problem::prescribedVelocity( m_system, m_ncomp, x[p], y[p], z[p] );
          // sum donain-edge contributions
          for (auto e : tk::cref_find(esued,{p,q})) {
            const auto [ N, grad, u, J ] =
              egrad( m_ncomp, m_offset, e, coord, inpoel, {}, U );
            auto J48 = J/48.0;
            for (const auto& [a,b] : tk::lpoed) {
              auto s = tk::orient( {N[a],N[b]}, {p,q} );
              for (std::size_t j=0; j<3; ++j) {
                for (std::size_t c=0; c<m_ncomp; ++c) {
                  R.var(r[c],p) -= J48 * s * (grad[a][j] - grad[b][j])
                                   * v[c][j]*(uL[c] + uR[c])
                    - J48 * std::abs(s * (grad[a][j] - grad[b][j]))
                          * std::abs(tk::dot(v[c],n)) * (uR[c] - uL[c]);
                }
              }
            }
          }
        }
      }

      // boundary-edge integrals
      for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
        // access node IDs
        const std::array< std::size_t, 3 >
          N{ tk::cref_find( lid, triinpoel[e*3+0] ),
             tk::cref_find( lid, triinpoel[e*3+1] ),
             tk::cref_find( lid, triinpoel[e*3+2] ) };
        // apply symmetry BCs
        if (bnorm.find(N[0]) != end(bnorm)) continue;
        // node coordinates
        std::array< tk::real, 3 > xp{ x[N[0]], x[N[1]], x[N[2]] },
                                  yp{ y[N[0]], y[N[1]], y[N[2]] },
                                  zp{ z[N[0]], z[N[1]], z[N[2]] };
        // compute face area
        auto A = tk::area( xp, yp, zp );
        auto A24 = A / 24.0;
        auto A6 = A / 6.0;
        // compute face normal
        auto n = tk::normal( xp, yp, zp );
        // access solution at element nodes
        std::vector< std::array< tk::real, 3 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, m_offset, N );
        // evaluate prescribed velocity
        auto v =
          Problem::prescribedVelocity( m_system, m_ncomp, xp[0], yp[0], zp[0] );
        // sum boundary integral contributions
        for (std::size_t c=0; c<m_ncomp; ++c) {
          auto vdotn = tk::dot( v[c], n );
          for (const auto& [a,b] : tk::lpoet) {
            auto Bab = A24 * vdotn * (u[c][a] + u[c][b]);
            R.var(r[c],N[a]) -= Bab + A6 * vdotn * u[c][a];
            R.var(r[c],N[b]) -= Bab;
          }
        }
      }
    }

    //! Compute right hand side for DiagCG (CG+FCT)
    //! \param[in] deltat Size of time step
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] Ue Element-centered solution vector at intermediate step
    //!    (used here internally as a scratch array)
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real,
              tk::real deltat,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const tk::Fields& U,
              tk::Fields& Ue,
              tk::Fields& R ) const
    {
      using tag::transport;
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      Assert( R.nunk() == coord[0].size(),
              "Number of unknowns in right-hand side vector incorrect" );

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // 1st stage: update element values from node values (gather-add)
      for (std::size_t e=0; e<inpoel.size()/4; ++e) {
        // access node IDs
        const std::array< std::size_t, 4 >
          N{{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] }};
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
        std::vector< std::array< tk::real, 4 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, m_offset, N );
        // access solution at elements
        std::vector< const tk::real* > ue( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) ue[c] = Ue.cptr( c, m_offset );

        // get prescribed velocity
        const std::array< std::vector<std::array<tk::real,3>>, 4 > vel{{
          Problem::prescribedVelocity( m_system, m_ncomp,
                                       x[N[0]], y[N[0]], z[N[0]] ),
          Problem::prescribedVelocity( m_system, m_ncomp,
                                       x[N[1]], y[N[1]], z[N[1]] ),
          Problem::prescribedVelocity( m_system, m_ncomp,
                                       x[N[2]], y[N[2]], z[N[2]] ),
          Problem::prescribedVelocity( m_system, m_ncomp,
                                       x[N[3]], y[N[3]], z[N[3]] ) }};

        // sum flux (advection) contributions to element
        auto d = deltat/2.0;
        for (std::size_t c=0; c<m_ncomp; ++c)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t a=0; a<4; ++a)
              Ue.var(ue[c],e) -= d * grad[a][j] * vel[a][c][j]*u[c][a];
      }


      // 2nd stage: form rhs from element values (scatter-add)
      for (std::size_t e=0; e<inpoel.size()/4; ++e) {
        // access node IDs
        const std::array< std::size_t, 4 >
          N{{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] }};
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
        std::vector< tk::real > ue( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) ue[c] = Ue( e, c, m_offset );
        // access pointer to right hand side at component and offset
        std::vector< const tk::real* > r( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) r[c] = R.cptr( c, m_offset );
        // access solution at nodes of element
        std::vector< std::array< tk::real, 4 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, m_offset, N );

        // get prescribed velocity
        auto xc = (x[N[0]] + x[N[1]] + x[N[2]] + x[N[3]]) / 4.0;
        auto yc = (y[N[0]] + y[N[1]] + y[N[2]] + y[N[3]]) / 4.0;
        auto zc = (z[N[0]] + z[N[1]] + z[N[2]] + z[N[3]]) / 4.0;
        const auto vel =
          Problem::prescribedVelocity( m_system, m_ncomp, xc, yc, zc );

        // scatter-add flux contributions to rhs at nodes
        tk::real d = deltat * J/6.0;
        for (std::size_t c=0; c<m_ncomp; ++c)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t a=0; a<4; ++a)
              R.var(r[c],N[a]) += d * grad[a][j] * vel[c][j]*ue[c];

        // add (optional) diffusion contribution to right hand side
        m_physics.diffusionRhs(m_system, m_ncomp, deltat, J, grad, N, u, r, R);
      }
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
      using tag::transport;
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
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
        std::vector< std::array< tk::real, 4 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, m_offset, N );
        // get velocity for problem
        const std::array< std::vector<std::array<tk::real,3>>, 4 > vel{{
          Problem::prescribedVelocity( m_system, m_ncomp,
                                       x[N[0]], y[N[0]], z[N[0]] ),
          Problem::prescribedVelocity( m_system, m_ncomp,
                                       x[N[1]], y[N[1]], z[N[1]] ),
          Problem::prescribedVelocity( m_system, m_ncomp,
                                       x[N[2]], y[N[2]], z[N[2]] ),
          Problem::prescribedVelocity( m_system, m_ncomp,
                                       x[N[3]], y[N[3]], z[N[3]] ) }};
        // compute the maximum length of the characteristic velocity (advection
        // velocity) across the four element nodes
        tk::real maxvel = 0.0;
        for (ncomp_t c=0; c<m_ncomp; ++c)
          for (std::size_t i=0; i<4; ++i) {
            auto v = tk::length( vel[i][c] );
            if (v > maxvel) maxvel = v;
          }
        // compute element dt for the advection
        auto advection_dt = L / maxvel;
        // compute element dt based on diffusion
        auto diffusion_dt = m_physics.diffusion_dt( m_system, m_ncomp, L, u );
        // compute minimum element dt
        auto elemdt = std::min( advection_dt, diffusion_dt );
        // find minimum dt across all elements
        if (elemdt < mindt) mindt = elemdt;
      }
      return mindt;
    }

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
    //! \param[in] t Physical time
    //! \param[in] deltat Time step size
    //! \param[in] ss Pair of side set ID and list of node IDs on the side set
    //! \param[in] coord Mesh node coordinates
    //! \return Vector of pairs of bool and boundary condition value associated
    //!   to mesh node IDs at which Dirichlet boundary conditions are set. Note
    //!   that instead of the actual boundary condition value, we return the
    //!   increment between t+dt and t, since that is what the solution requires
    //!   as we solve for the soution increments and not the solution itself.
    std::map< std::size_t, std::vector< std::pair<bool,tk::real> > >
    dirbc( tk::real t,
           tk::real deltat,
           const std::pair< const int, std::vector< std::size_t > >& ss,
           const std::array< std::vector< tk::real >, 3 >& coord ) const
    {
      using tag::param; using tag::transport; using tag::bcdir;
      using NodeBC = std::vector< std::pair< bool, tk::real > >;
      std::map< std::size_t, NodeBC > bc;
      const auto& ubc = g_inputdeck.get< param, transport, tag::bc, bcdir >();
      if (!ubc.empty()) {
        Assert( ubc.size() > m_system, "Indexing out of Dirichlet BC eq-vector" );
        const auto& x = coord[0];
        const auto& y = coord[1];
        const auto& z = coord[2];
        for (const auto& b : ubc[m_system])
          if (std::stoi(b) == ss.first)
            for (auto n : ss.second) {
              Assert( x.size() > n, "Indexing out of coordinate array" );
              const auto s = solinc( m_system, m_ncomp, x[n], y[n], z[n],
                                     t, deltat, Problem::solution );
              auto& nbc = bc[n] = NodeBC( m_ncomp );
              for (ncomp_t c=0; c<m_ncomp; ++c)
                nbc[c] = { true, s[c] };
            }
      }
      return bc;
    }

    //! Set symmetry boundary conditions at nodes
    void
    symbc( tk::Fields&,
           const std::unordered_map<std::size_t,std::array<tk::real,4>>& )
    const {}

    //! Query nodes at which symmetry boundary conditions are set
    //! \param[in] bface Boundary-faces mapped to side set ids
    //! \param[in] triinpoel Boundary-face connectivity
    //! \param[in,out] nodes Node ids at which symmetry BCs are set
    void
    symbcnodes( const std::map< int, std::vector< std::size_t > >& bface,
                const std::vector< std::size_t >& triinpoel,
                std::unordered_set< std::size_t >& nodes ) const
    {
      using tag::param; using tag::transport; using tag::bcsym;
      const auto& bc = g_inputdeck.get< param, transport, tag::bc, bcsym >();
      if (!bc.empty() && bc.size() > m_system) {
        const auto& ss = bc[ m_system ];// side sets with sym bcs specified
        for (const auto& s : ss) {
          auto k = bface.find( std::stoi(s) );
          if (k != end(bface)) {
            for (auto f : k->second) {  // face ids on symbc side set
              nodes.insert( triinpoel[f*3+0] );
              nodes.insert( triinpoel[f*3+1] );
              nodes.insert( triinpoel[f*3+2] );
            }
          }
        }
      }
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    //! \details This functions should be written in conjunction with
    //!   fieldOutput(), which provides the vector of fields to be output
    std::vector< std::string > fieldNames() const {
      std::vector< std::string > n;
      const auto& depvar =
        g_inputdeck.get< tag::param, tag::transport, tag::depvar >().
        at(m_system);
      // will output numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_numerical" );
      // will output analytic solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_analytic" );
      // will output error for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_error" );
      return n;
    }

    //! Return surface field names to be output to file
    //! \return Vector of strings labelling surface fields output in file
    //! \details This functions should be written in conjunction with
    //!   surfOutput(), which provides the vector of surface fields to be output
    std::vector< std::string > surfNames() const {
      std::vector< std::string > n;     // punt for now
      return n;
    }

    //! Return surface field output going to file
    std::vector< std::vector< tk::real > >
    surfOutput( const std::map< int, std::vector< std::size_t > >&,
                tk::Fields& ) const
    {
      std::vector< std::vector< tk::real > > s; // punt for now
      return s;
    }

    //! Return time history field names to be output to file
    //! \return Vector of strings labelling time history fields output in file
    std::vector< std::string > histNames() const {
      std::vector< std::string > s; // punt for now
      return s;
    }

    //! Return time history field output evaluated at time history points
    std::vector< std::vector< tk::real > >
    histOutput( const std::vector< HistData >&,
                const std::vector< std::size_t >&,
                const tk::Fields& ) const
    {
      std::vector< std::vector< tk::real > > s; // punt for now
      return s;
    }

    //! Return field output going to file
    //! \param[in] t Physical time
    //! \param[in] V Total mesh volume
    //! \param[in] coord Mesh node coordinates
    //! \param[in] v Nodal volumes
    //! \param[in,out] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    //! \details This functions should be written in conjunction with names(),
    //!   which provides the vector of field names
    //! \note U is overwritten
    std::vector< std::vector< tk::real > >
    fieldOutput( tk::real t,
                 tk::real V,
                 const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< tk::real >& v,
                 tk::Fields& U ) const
    {
      std::vector< std::vector< tk::real > > out;
      // will output numerical solution for all components
      auto E = U;
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( U.extract( c, m_offset ) );
      // evaluate analytic solution at time t
      std::vector< std::size_t > inbox;
      initialize( coord, U, t, inbox );
      // will output analytic solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( U.extract( c, m_offset ) );
      // will output error for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) {
        auto u = U.extract( c, m_offset );
        auto e = E.extract( c, m_offset );
        Assert( u.size() == e.size(), "Size mismatch" );
        Assert( u.size() == v.size(), "Size mismatch" );
        for (std::size_t i=0; i<u.size(); ++i)
          e[i] = std::pow( e[i] - u[i], 2.0 ) * v[i] / V;
        out.push_back( e );
      }
      return out;
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const {
      std::vector< std::string > n;
      const auto& depvar =
        g_inputdeck.get< tag::param, tag::transport, tag::depvar >().
        at(m_system);
      // construct the name of the numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) );
      return n;
    }

  private:
    const Physics m_physics;            //!< Physics policy
    const Problem m_problem;            //!< Problem policy
    const ncomp_t m_system;             //!< Equation system index
    const ncomp_t m_ncomp;              //!< Number of components in this PDE
    const ncomp_t m_offset;             //!< Offset this PDE operates from

    //! Compute element contribution to nodal gradient
    //! \param[in] e Element whose contribution to compute
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] U Solution vector at recent time step
    //! \return Tuple of element contribution
    //! \note The function signature must follow tk::ElemGradFn
    static tk::ElemGradFn::result_type
    egrad( ncomp_t ncomp,
           ncomp_t offset,
           std::size_t e,
           const std::array< std::vector< tk::real >, 3 >& coord,
           const std::vector< std::size_t >& inpoel,
           const std::tuple< std::vector< tk::real >,
                             std::vector< tk::real > >&,
           const tk::Fields& U )
    {
      // access node cooordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
      // access node IDs
      const std::array< std::size_t, 4 >
        N{{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] }};
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
      std::vector< std::array< tk::real, 4 > > u( ncomp );
      for (ncomp_t c=0; c<ncomp; ++c) u[c] = U.extract( c, offset, N );
      // return data needed to scatter add element contribution to gradient
      return std::tuple< std::array< std::size_t, 4 >,
                         std::array< std::array< tk::real, 3 >, 4 >,
                         std::vector< std::array< tk::real, 4 > >,
                         tk::real >
               { std::move(N), std::move(grad), std::move(u), std::move(J) };
    }
};

} // cg::
} // inciter::

#endif // Transport_h
