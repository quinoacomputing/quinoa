// *****************************************************************************
/*!
  \file      src/PDE/Transport/CGTransport.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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
    using ncomp_t = tk::ncomp_t;
    using real = tk::real;
    using eq = tag::transport;

    static constexpr real muscl_eps = 1.0e-9;
    static constexpr real muscl_const = 1.0/3.0;
    static constexpr real muscl_m1 = 1.0 - muscl_const;
    static constexpr real muscl_p1 = 1.0 + muscl_const;

  public:
    //! Constructor
    explicit Transport() :
      m_physics( Physics() ),
      m_problem( Problem() ),
      m_ncomp(g_inputdeck.get< tag::ncomp >())
    {
      m_problem.errchk( m_ncomp );
    }

    //! Determine nodes that lie inside the user-defined IC box
    void
    IcBoxNodes( const tk::UnsMesh::Coords&,
                const std::vector< std::size_t >&,
                const std::unordered_map< std::size_t, std::set< std::size_t > >&,
                std::vector< std::unordered_set< std::size_t > >&,
                std::unordered_map< std::size_t, std::set< std::size_t > >&,
                std::size_t& ) const {}

    //! Initalize the transport equations using problem policy
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    void
    initialize( const std::array< std::vector< real >, 3 >& coord,
                tk::Fields& unk,
                real t,
                real,
                const std::vector< std::unordered_set< std::size_t > >&,
                const std::vector< tk::real >&,
                const std::unordered_map< std::size_t, std::set< std::size_t > >&
              ) const
    {
      Assert( coord[0].size() == unk.nunk(), "Size mismatch" );
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
      for (ncomp_t i=0; i<x.size(); ++i) {
        auto s = Problem::initialize( m_ncomp, m_mat_blk, x[i], y[i],
                                      z[i], t );
        for (ncomp_t c=0; c<m_ncomp; ++c)
          unk( i, c ) = s[c];
      }
    }

    //! Query a velocity
    //! \note Since this function does not touch its output argument, that
    //!   means this system does not define a "velocity".
    void velocity( const tk::Fields&, tk::UnsMesh::Coords& ) const {}

    //! Query the sound speed
    //! \note Since this function does not touch its output argument, that
    //!   means this system does not define a "sound speed".
    void soundspeed( const tk::Fields&, std::vector< tk::real >& ) const {}

    //! Return analytic solution (if defined by Problem) at xi, yi, zi, t
    //! \param[in] xi X-coordinate
    //! \param[in] yi Y-coordinate
    //! \param[in] zi Z-coordinate
    //! \param[in] t Physical time
    //! \return Vector of analytic solution at given location and time
    std::vector< real >
    analyticSolution( real xi, real yi, real zi, real t ) const
    { return Problem::analyticSolution( m_ncomp, m_mat_blk, xi, yi, zi, t ); }

    //! Return analytic solution for conserved variables
    //! \param[in] xi X-coordinate at which to evaluate the analytic solution
    //! \param[in] yi Y-coordinate at which to evaluate the analytic solution
    //! \param[in] zi Z-coordinate at which to evaluate the analytic solution
    //! \param[in] t Physical time at which to evaluate the analytic solution
    //! \return Vector of analytic solution at given location and time
    std::vector< tk::real >
    solution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return Problem::initialize( m_ncomp, m_mat_blk, xi, yi, zi, t ); }

    //! Compute nodal gradients of primitive variables for ALECG
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] bndel List of elements contributing to chare-boundary nodes
    //! \param[in] gid Local->global node id map
    //! \param[in] bid Local chare-boundary node ids (value) associated to
    //!    global node ids (key)
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] G Nodal gradients of primitive variables
    void chBndGrad( const std::array< std::vector< real >, 3 >& coord,
      const std::vector< std::size_t >& inpoel,
      const std::vector< std::size_t >& bndel,
      const std::vector< std::size_t >& gid,
      const std::unordered_map< std::size_t, std::size_t >& bid,
      const tk::Fields& U,
      tk::Fields& G ) const
    {
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );

      // compute gradients of primitive variables in points
      G.fill( 0.0 );

      // access node cooordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      for (auto e : bndel) {  // elements contributing to chare boundary nodes
        // access node IDs
        std::size_t N[4] =
          { inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
        // compute element Jacobi determinant, J = 6V
        real bax = x[N[1]]-x[N[0]];
        real bay = y[N[1]]-y[N[0]];
        real baz = z[N[1]]-z[N[0]];
        real cax = x[N[2]]-x[N[0]];
        real cay = y[N[2]]-y[N[0]];
        real caz = z[N[2]]-z[N[0]];
        real dax = x[N[3]]-x[N[0]];
        real day = y[N[3]]-y[N[0]];
        real daz = z[N[3]]-z[N[0]];
        auto J = tk::triple( bax, bay, baz, cax, cay, caz, dax, day, daz );
        auto J24 = J/24.0;
        // shape function derivatives, nnode*ndim [4][3]
        real g[4][3];
        tk::crossdiv( cax, cay, caz, dax, day, daz, J,
                      g[1][0], g[1][1], g[1][2] );
        tk::crossdiv( dax, day, daz, bax, bay, baz, J,
                      g[2][0], g[2][1], g[2][2] );
        tk::crossdiv( bax, bay, baz, cax, cay, caz, J,
                      g[3][0], g[3][1], g[3][2] );
        for (std::size_t i=0; i<3; ++i)
          g[0][i] = -g[1][i] - g[2][i] - g[3][i];
        // scatter-add gradient contributions to boundary nodes
        for (std::size_t a=0; a<4; ++a) {
          auto i = bid.find( gid[N[a]] );
          if (i != end(bid))
            for (std::size_t c=0; c<m_ncomp; ++c)
              for (std::size_t b=0; b<4; ++b)
                for (std::size_t j=0; j<3; ++j)
                  G(i->second,c*3+j) += J24 * g[b][j] * U(N[b],c);
        }
      }
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
    //! \param[in] esup Elements surrounding points
    //! \param[in] symbctri Vector with 1 at symmetry BC nodes
    //! \param[in] vol Nodal volumes
    //! \param[in] edgeid Local node id pair -> edge id map
    //! \param[in] G Nodal gradients in chare-boundary nodes
    //! \param[in] U Solution vector at recent time step
    //! \param[in] W Mesh velocity
    //! \param[in,out] R Right-hand side vector computed
    void rhs(
      real,
      const std::array< std::vector< real >, 3 >&  coord,
      const std::vector< std::size_t >& inpoel,
      const std::vector< std::size_t >& triinpoel,
      const std::vector< std::size_t >&,
      const std::unordered_map< std::size_t, std::size_t >& bid,
      const std::unordered_map< std::size_t, std::size_t >& lid,
      const std::vector< real >& dfn,
      const std::pair< std::vector< std::size_t >,
                       std::vector< std::size_t > >& psup,
      const std::pair< std::vector< std::size_t >,
                       std::vector< std::size_t > >& esup,
      const std::vector< int >& symbctri,
      const std::vector< real >& vol,
      const std::vector< std::size_t >&,
      const std::vector< std::size_t >& edgeid,
      const std::vector< std::unordered_set< std::size_t > >&,
      const tk::Fields& G,
      const tk::Fields& U,
      [[maybe_unused]] const tk::Fields& W,
      const std::vector< tk::real >&,
      real,
      tk::Fields& R ) const
    {
      Assert( G.nprop() == m_ncomp*3,
              "Number of components in gradient vector incorrect" );
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      Assert( R.nunk() == coord[0].size(),
              "Number of unknowns and/or number of components in right-hand "
              "side vector incorrect" );

      // compute/assemble gradients in points
      auto Grad = nodegrad( coord, inpoel, lid, bid, vol, esup, U, G );

      // zero right hand side for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) R.fill( c, 0.0 );

      // compute domain-edge integral
      domainint( coord, inpoel, edgeid, psup, dfn, U, Grad, R );

      // compute boundary integrals
      bndint( coord, triinpoel, symbctri, U, R );
    }

    //! Compute boundary pressure integrals (force) (no-op for transport)
    void bndPressureInt(
      const std::array< std::vector< real >, 3 >&,
      const std::vector< std::size_t >&,
      const std::vector< int >&,
      const tk::Fields&,
      const std::array< tk::real, 3 >&,
      std::vector< real >& ) const
    { }


    //! Compute the minimum time step size (for unsteady time stepping)
    //! \param[in] U Solution vector at recent time step
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] t Physical time
    //! \return Minimum time step size
    real dt( const std::array< std::vector< real >, 3 >& coord,
             const std::vector< std::size_t >& inpoel,
             tk::real t,
             tk::real,
             const tk::Fields& U,
             const std::vector< tk::real >&,
             const std::vector< tk::real >& ) const
    {
      using tag::transport;
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
      // compute the minimum dt across all elements we own
      auto mindt = std::numeric_limits< tk::real >::max();
      auto eps = std::numeric_limits< tk::real >::epsilon();
      auto large = std::numeric_limits< tk::real >::max();
      for (std::size_t e=0; e<inpoel.size()/4; ++e) {
        const std::array< std::size_t, 4 >
          N{{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] }};
        // compute cubic root of element volume as the characteristic length
        const std::array< real, 3 >
          ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
          ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
          da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
        const auto L = std::cbrt( tk::triple( ba, ca, da ) / 6.0 );
        // access solution at element nodes at recent time step
        std::vector< std::array< real, 4 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, N );
        // get velocity for problem
        const std::array< std::vector<std::array<real,3>>, 4 > vel{{
          Problem::prescribedVelocity( m_ncomp,
                                       x[N[0]], y[N[0]], z[N[0]], t ),
          Problem::prescribedVelocity( m_ncomp,
                                       x[N[1]], y[N[1]], z[N[1]], t ),
          Problem::prescribedVelocity( m_ncomp,
                                       x[N[2]], y[N[2]], z[N[2]], t ),
          Problem::prescribedVelocity( m_ncomp,
                                       x[N[3]], y[N[3]], z[N[3]], t ) }};
        // compute the maximum length of the characteristic velocity (advection
        // velocity) across the four element nodes
        real maxvel = 0.0;
        for (ncomp_t c=0; c<m_ncomp; ++c)
          for (std::size_t i=0; i<4; ++i) {
            auto v = tk::length( vel[i][c] );
            if (v > maxvel) maxvel = v;
          }
        // compute element dt for the advection
        auto advection_dt = std::abs(maxvel) > eps ? L / maxvel : large;
        // compute element dt based on diffusion
        auto diffusion_dt = m_physics.diffusion_dt( m_ncomp, L, u );
        // compute minimum element dt
        auto elemdt = std::min( advection_dt, diffusion_dt );
        // find minimum dt across all elements
        if (elemdt < mindt) mindt = elemdt;
      }
      return mindt * g_inputdeck.get< tag::cfl >();
    }

    //! Compute a time step size for each mesh node (for steady time stepping)
    void dt( uint64_t,
             const std::vector< tk::real >&,
             const tk::Fields&,
             std::vector< tk::real >& ) const {}

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
    //! \param[in] t Physical time
    //! \param[in] deltat Time step size
    //! \param[in] tp Physical time for each mesh node
    //! \param[in] dtp Time step size for each mesh node
    //! \param[in] ss Pair of side set ID and list of node IDs on the side set
    //! \param[in] coord Mesh node coordinates
    //! \param[in] increment If true, evaluate the solution increment between
    //!   t and t+dt for Dirichlet BCs. If false, evlauate the solution instead.
    //! \return Vector of pairs of bool and boundary condition value associated
    //!   to mesh node IDs at which Dirichlet boundary conditions are set. Note
    //!   that if increment is true, instead of the actual boundary condition
    //!   value, we return the increment between t+deltat and t, since,
    //!   depending on client code and solver, that may be what the solution
    //!   requires.
    std::map< std::size_t, std::vector< std::pair<bool,real> > >
    dirbc( real t,
           real deltat,
           const std::vector< tk::real >& tp,
           const std::vector< tk::real >& dtp,
           const std::pair< const int, std::vector< std::size_t > >& ss,
           const std::array< std::vector< real >, 3 >& coord,
           bool increment ) const
    {
      using NodeBC = std::vector< std::pair< bool, real > >;
      std::map< std::size_t, NodeBC > bc;
      const auto& ubc = g_inputdeck.get< tag::bc >()[0].get< tag::dirichlet >();
      const auto steady = g_inputdeck.get< tag::steady_state >();
      if (!ubc.empty()) {
        Assert( ubc.size() > 0, "Indexing out of Dirichlet BC eq-vector" );
        const auto& x = coord[0];
        const auto& y = coord[1];
        const auto& z = coord[2];
        for (const auto& b : ubc)
          if (static_cast<int>(b) == ss.first)
            for (auto n : ss.second) {
              Assert( x.size() > n, "Indexing out of coordinate array" );
              if (steady) { t = tp[n]; deltat = dtp[n]; }
              const auto s = increment ?
                solinc( m_ncomp, m_mat_blk, x[n], y[n], z[n],
                        t, deltat, Problem::initialize ) :
                Problem::initialize( m_ncomp, m_mat_blk, x[n], y[n],
                                     z[n], t+deltat );
              auto& nbc = bc[n] = NodeBC( m_ncomp );
              for (ncomp_t c=0; c<m_ncomp; ++c)
                nbc[c] = { true, s[c] };
            }
      }
      return bc;
    }

    //! Set symmetry boundary conditions at nodes
    void
    symbc(
      tk::Fields&,
      tk::Fields&,
      const std::array< std::vector< real >, 3 >&,
      const std::unordered_map< int,
              std::unordered_map< std::size_t,
                std::array< real, 4 > > >&,
      const std::unordered_set< std::size_t >& ) const {}

    //! Set farfield boundary conditions at nodes
    void farfieldbc(
      tk::Fields&,
      const std::array< std::vector< real >, 3 >&,
      const std::unordered_map< int,
              std::unordered_map< std::size_t,
                std::array< real, 4 > > >&,
      const std::unordered_set< std::size_t >& ) const {}

    //! Apply user defined time dependent BCs (no-op for transport)
    void
    timedepbc( tk::real,
      tk::Fields&,
      const std::vector< std::unordered_set< std::size_t > >&,
      const std::vector< tk::Table<5> >& ) const {}

    //! Return a map that associates user-specified strings to functions
    //! \return Map that associates user-specified strings to functions that
    //!  compute relevant quantities to be output to file
    std::map< std::string, tk::GetVarFn > OutVarFn() const {
      std::map< std::string, tk::GetVarFn > OutFnMap;
      OutFnMap["material_indicator"] = transport::matIndicatorOutVar;

      return OutFnMap;
    }

    //! Return analytic field names to be output to file
    //! \return Vector of strings labelling analytic fields output in file
    std::vector< std::string > analyticFieldNames() const {
      std::vector< std::string > n;
      auto depvar = g_inputdeck.get< tag::depvar >()[0];
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_analytic" );
      return n;
    }

    //! Return surface field names to be output to file
    //! \return Vector of strings labelling surface fields output in file
    //! \details This functions should be written in conjunction with
    //!   surfOutput(), which provides the vector of surface fields to be output
    std::vector< std::string > surfNames() const { return {}; }

    //! Return nodal surface field output going to file
    std::vector< std::vector< real > >
    surfOutput( const std::map< int, std::vector< std::size_t > >&,
                const tk::Fields& ) const { return {}; }

    //! Return elemental surface field output (on triangle faces) going to file
    std::vector< std::vector< real > >
    elemSurfOutput( const std::map< int, std::vector< std::size_t > >&,
      const std::vector< std::size_t >&,
      const tk::Fields& ) const { return {}; }

    //! Return time history field names to be output to file
    //! \return Vector of strings labelling time history fields output in file
    std::vector< std::string > histNames() const { return {}; }

    //! Return time history field output evaluated at time history points
    std::vector< std::vector< real > >
    histOutput( const std::vector< HistData >&,
                const std::vector< std::size_t >&,
                const tk::Fields& ) const { return {}; }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const {
      std::vector< std::string > n;
      const auto& depvar =
        g_inputdeck.get< tag::depvar >().at(0);
      // construct the name of the numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) );
      return n;
    }

  private:
    const Physics m_physics;            //!< Physics policy
    const Problem m_problem;            //!< Problem policy
    const ncomp_t m_ncomp;              //!< Number of components in this PDE
    //! EOS material block
    const std::vector< EOS > m_mat_blk;

    //! \brief Compute/assemble nodal gradients of primitive variables for
    //!   ALECG in all points
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] lid Global->local node ids
    //! \param[in] bid Local chare-boundary node ids (value) associated to
    //!    global node ids (key)
    //! \param[in] vol Nodal volumes
    //! \param[in] esup Elements surrounding points
    //! \param[in] U Solution vector at recent time step
    //! \param[in] G Nodal gradients of primitive variables in chare-boundary nodes
    //! \return Gradients of primitive variables in all mesh points
    tk::Fields
    nodegrad( const std::array< std::vector< real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const std::unordered_map< std::size_t, std::size_t >& lid,
              const std::unordered_map< std::size_t, std::size_t >& bid,
              const std::vector< real >& vol,
              const std::pair< std::vector< std::size_t >,
                               std::vector< std::size_t > >& esup,
              const tk::Fields& U,
              const tk::Fields& G ) const
    {
      // allocate storage for nodal gradients of primitive variables
      tk::Fields Grad( U.nunk(), m_ncomp*3 );
      Grad.fill( 0.0 );

      // access node cooordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // compute gradients of primitive variables in points
      auto npoin = U.nunk();
      #pragma omp simd
      for (std::size_t p=0; p<npoin; ++p)
        for (auto e : tk::Around(esup,p)) {
          // access node IDs
          std::size_t N[4] =
            { inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
          // compute element Jacobi determinant, J = 6V
          real bax = x[N[1]]-x[N[0]];
          real bay = y[N[1]]-y[N[0]];
          real baz = z[N[1]]-z[N[0]];
          real cax = x[N[2]]-x[N[0]];
          real cay = y[N[2]]-y[N[0]];
          real caz = z[N[2]]-z[N[0]];
          real dax = x[N[3]]-x[N[0]];
          real day = y[N[3]]-y[N[0]];
          real daz = z[N[3]]-z[N[0]];
          auto J = tk::triple( bax, bay, baz, cax, cay, caz, dax, day, daz );
          auto J24 = J/24.0;
          // shape function derivatives, nnode*ndim [4][3]
          real g[4][3];
          tk::crossdiv( cax, cay, caz, dax, day, daz, J,
                        g[1][0], g[1][1], g[1][2] );
          tk::crossdiv( dax, day, daz, bax, bay, baz, J,
                        g[2][0], g[2][1], g[2][2] );
          tk::crossdiv( bax, bay, baz, cax, cay, caz, J,
                        g[3][0], g[3][1], g[3][2] );
          for (std::size_t i=0; i<3; ++i)
            g[0][i] = -g[1][i] - g[2][i] - g[3][i];
          // scatter-add gradient contributions to boundary nodes
          for (std::size_t c=0; c<m_ncomp; ++c)
            for (std::size_t b=0; b<4; ++b)
              for (std::size_t i=0; i<3; ++i)
                Grad(p,c*3+i) += J24 * g[b][i] * U(N[b],c);
        }

      // put in nodal gradients of chare-boundary points
      for (const auto& [g,b] : bid) {
        auto i = tk::cref_find( lid, g );
        for (ncomp_t c=0; c<Grad.nprop(); ++c)
          Grad(i,c) = G(b,c);
      }

      // divide weak result in gradients by nodal volume
      for (std::size_t p=0; p<npoin; ++p)
        for (std::size_t c=0; c<m_ncomp*3; ++c)
          Grad(p,c) /= vol[p];

      return Grad;
    }

    //! \brief Compute MUSCL reconstruction in edge-end points using a MUSCL
    //!    procedure with van Leer limiting
    //! \param[in] p Left node id of edge-end
    //! \param[in] q Right node id of edge-end
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] G Gradient of all unknowns in mesh points
    //! \param[in,out] uL Primitive variables at left edge-end point
    //! \param[in,out] uR Primitive variables at right edge-end point
    void
    muscl( std::size_t p,
           std::size_t q,
           const tk::UnsMesh::Coords& coord,
           const tk::Fields& G,
           std::vector< real >& uL,
           std::vector< real >& uR ) const
    {
      Assert( uL.size() == m_ncomp && uR.size() == m_ncomp, "Size mismatch" );
      Assert( G.nprop()/3 == m_ncomp, "Size mismatch" );

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // edge vector
      std::array< real, 3 > vw{ x[q]-x[p], y[q]-y[p], z[q]-z[p] };

      std::vector< real >
        delta1( m_ncomp, 0.0 ), delta2( m_ncomp, 0.0 ), delta3( m_ncomp, 0.0 );

      // MUSCL reconstruction of edge-end-point primitive variables
      for (std::size_t c=0; c<m_ncomp; ++c) {
        // gradients
        std::array< real, 3 >
          g1{ G(p,c*3+0), G(p,c*3+1), G(p,c*3+2) },
          g2{ G(q,c*3+0), G(q,c*3+1), G(q,c*3+2) };

        delta2[c] = uR[c] - uL[c];
        delta1[c] = 2.0 * tk::dot(g1,vw) - delta2[c];
        delta3[c] = 2.0 * tk::dot(g2,vw) - delta2[c];

        // form limiters
        auto rL = (delta2[c] + muscl_eps) / (delta1[c] + muscl_eps);
        auto rR = (delta2[c] + muscl_eps) / (delta3[c] + muscl_eps);
        auto rLinv = (delta1[c] + muscl_eps) / (delta2[c] + muscl_eps);
        auto rRinv = (delta3[c] + muscl_eps) / (delta2[c] + muscl_eps);

        auto phiL = (std::abs(rL) + rL) / (std::abs(rL) + 1.0);
        auto phiR = (std::abs(rR) + rR) / (std::abs(rR) + 1.0);
        auto phi_L_inv = (std::abs(rLinv) + rLinv) / (std::abs(rLinv) + 1.0);
        auto phi_R_inv = (std::abs(rRinv) + rRinv) / (std::abs(rRinv) + 1.0);

        // update unknowns with reconstructed unknowns
        uL[c] += 0.25*(delta1[c]*muscl_m1*phiL + delta2[c]*muscl_p1*phi_L_inv);
        uR[c] -= 0.25*(delta3[c]*muscl_m1*phiR + delta2[c]*muscl_p1*phi_R_inv);
      }
    }

    //! Compute domain-edge integral for ALECG
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] edgeid Local node id pair -> edge id map
    //! \param[in] psup Points surrounding points
    //! \param[in] dfn Dual-face normals
    //! \param[in] U Solution vector at recent time step
    //! \param[in] G Nodal gradients
    //! \param[in,out] R Right-hand side vector computed
    void domainint( const std::array< std::vector< real >, 3 >& coord,
                    const std::vector< std::size_t >& inpoel,
                    const std::vector< std::size_t >& edgeid,
                    const std::pair< std::vector< std::size_t >,
                                     std::vector< std::size_t > >& psup,
                    const std::vector< real >& dfn,
                    const tk::Fields& U,
                    const tk::Fields& G,
                    tk::Fields& R ) const
    {
      // access node cooordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // compute derived data structures
      auto esued = tk::genEsued( inpoel, 4, tk::genEsup( inpoel, 4 ) );

      // access pointer to right hand side at component
      std::vector< const real* > r( m_ncomp );
      for (ncomp_t c=0; c<m_ncomp; ++c) r[c] = R.cptr( c );

      // domain-edge integral
      for (std::size_t p=0,k=0; p<U.nunk(); ++p) {
        for (auto q : tk::Around(psup,p)) {
          // access dual-face normals for edge p-q
          auto ed = edgeid[k++];
          std::array< tk::real, 3 > n{ dfn[ed*6+0], dfn[ed*6+1], dfn[ed*6+2] };

          std::vector< tk::real > uL( m_ncomp, 0.0 );
          std::vector< tk::real > uR( m_ncomp, 0.0 );
          for (std::size_t c=0; c<m_ncomp; ++c) {
            uL[c] = U(p,c);
            uR[c] = U(q,c);
          }
          // compute MUSCL reconstruction in edge-end points
          muscl( p, q, coord, G, uL, uR );

          // evaluate prescribed velocity
          auto v =
            Problem::prescribedVelocity( m_ncomp, x[p], y[p], z[p], 0.0 );
          // sum donain-edge contributions
          for (auto e : tk::cref_find(esued,{p,q})) {
            const std::array< std::size_t, 4 >
              N{{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] }};
            // compute element Jacobi determinant
            const std::array< tk::real, 3 >
              ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
              ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
              da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
            const auto J = tk::triple( ba, ca, da );        // J = 6V
            // shape function derivatives, nnode*ndim [4][3]
            std::array< std::array< tk::real, 3 >, 4 > grad;
            grad[1] = tk::crossdiv( ca, da, J );
            grad[2] = tk::crossdiv( da, ba, J );
            grad[3] = tk::crossdiv( ba, ca, J );
            for (std::size_t i=0; i<3; ++i)
              grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];
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
    }

    //! Compute boundary integrals for ALECG
    //! \param[in] coord Mesh node coordinates
    //! \param[in] triinpoel Boundary triangle face connecitivity with local ids
    //! \param[in] symbctri Vector with 1 at symmetry BC boundary triangles
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    void bndint( const std::array< std::vector< real >, 3 >& coord,
                 const std::vector< std::size_t >& triinpoel,
                 const std::vector< int >& symbctri,
                 const tk::Fields& U,
                 tk::Fields& R ) const
    {
      // access node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // boundary integrals: compute fluxes in edges
      std::vector< real > bflux( triinpoel.size() * m_ncomp * 2 );

      for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
        // access node IDs
        std::array< std::size_t, 3 >
          N{ triinpoel[e*3+0], triinpoel[e*3+1], triinpoel[e*3+2] };
        // apply symmetry BCs
        if (symbctri[e]) continue;
        // node coordinates
        std::array< tk::real, 3 > xp{ x[N[0]], x[N[1]], x[N[2]] },
                                  yp{ y[N[0]], y[N[1]], y[N[2]] },
                                  zp{ z[N[0]], z[N[1]], z[N[2]] };
        // access solution at element nodes
        std::vector< std::array< real, 3 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, N );
        // evaluate prescribed velocity
        auto v =
          Problem::prescribedVelocity( m_ncomp, xp[0], yp[0], zp[0], 0.0 );
        // compute face area
        auto A6 = tk::area( x[N[0]], x[N[1]], x[N[2]],
                            y[N[0]], y[N[1]], y[N[2]],
                            z[N[0]], z[N[1]], z[N[2]] ) / 6.0;
        auto A24 = A6/4.0;
        // compute face normal
        auto n = tk::normal( xp, yp, zp );
        // store flux in boundary elements
        for (std::size_t c=0; c<m_ncomp; ++c) {
          auto eb = (e*m_ncomp+c)*6;
          auto vdotn = tk::dot( v[c], n );
          auto Bab = A24 * vdotn * (u[c][0] + u[c][1]);
          bflux[eb+0] = Bab + A6 * vdotn * u[c][0];
          bflux[eb+1] = Bab;
          Bab = A24 * vdotn * (u[c][1] + u[c][2]);
          bflux[eb+2] = Bab + A6 * vdotn * u[c][1];
          bflux[eb+3] = Bab;
          Bab = A24 * vdotn * (u[c][2] + u[c][0]);
          bflux[eb+4] = Bab + A6 * vdotn * u[c][2];
          bflux[eb+5] = Bab;
        }
      }

      // access pointer to right hand side at component
      std::vector< const real* > r( m_ncomp );
      for (ncomp_t c=0; c<m_ncomp; ++c) r[c] = R.cptr( c );

      // boundary integrals: sum flux contributions to points
      for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
        std::size_t N[3] =
          { triinpoel[e*3+0], triinpoel[e*3+1], triinpoel[e*3+2] };
        for (std::size_t c=0; c<m_ncomp; ++c) {
          auto eb = (e*m_ncomp+c)*6;
          R.var(r[c],N[0]) -= bflux[eb+0] + bflux[eb+5];
          R.var(r[c],N[1]) -= bflux[eb+1] + bflux[eb+2];
          R.var(r[c],N[2]) -= bflux[eb+3] + bflux[eb+4];
        }
      }

      tk::destroy(bflux);
    }
};

} // cg::
} // inciter::

#endif // Transport_h
