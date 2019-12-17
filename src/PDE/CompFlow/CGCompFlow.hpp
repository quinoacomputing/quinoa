// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/CGCompFlow.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include "DerivedData.hpp"
#include "Exception.hpp"
#include "Vector.hpp"
#include "EoS/EoS.hpp"
#include "Mesh/Around.hpp"
#include "Reconstruction.hpp"
#include "Integrate/Riemann/HLLC.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

static constexpr tk::real muscl_eps = 1.0e-9;
static constexpr tk::real muscl_const = 1.0/3.0;
static constexpr tk::real muscl_m1 = 1.0 - muscl_const;
static constexpr tk::real muscl_p1 = 1.0 + muscl_const;

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
    //! \param[in] c Equation system index (among multiple systems configured)
    explicit CompFlow( ncomp_t c ) :
      m_physics(),
      m_problem(),
      m_system( c ),
      m_ncomp(
        g_inputdeck.get< tag::component >().get< tag::compflow >().at(c) ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< tag::compflow >(c) )
    {
       Assert( m_ncomp == 5, "Number of CompFlow PDE components must be 5" );
    }

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
        const auto s =
          Problem::solution( m_system, m_ncomp, x[i], y[i], z[i], t );
        unk(i,0,m_offset) = s[0]; // rho
        unk(i,1,m_offset) = s[1]; // rho * u
        unk(i,2,m_offset) = s[2]; // rho * v
        unk(i,3,m_offset) = s[3]; // rho * w
        unk(i,4,m_offset) = s[4]; // rho * e, e: total = kinetic + internal
      }
    }

    //! Return analytic solution (if defined by Problem) at xi, yi, zi, t
    //! \param[in] xi X-coordinate
    //! \param[in] yi Y-coordinate
    //! \param[in] zi Z-coordinate
    //! \param[in] t Physical time
    //! \return Vector of analytic solution at given location and time
    std::vector< tk::real >
    analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    {
      auto s = Problem::solution( m_system, m_ncomp, xi, yi, zi, t );
      return std::vector< tk::real >( begin(s), end(s) );
    }

    //! Compute nodal gradients of primitive variables for ALECG
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
//    //! \param[in] bndel List of elements contributing to chare-boundary nodes
//    //! \param[in] gid Local->global node id map
//    //! \param[in] bid Local chare-boundary node ids (value) associated to
//    //!    global node ids (key)
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
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      Assert( G.nprop() == m_ncomp*3,
              "Number of components in gradient vector incorrect" );

      // compute gradients of primitive variables in points
      G.fill( 0.0 );

      for (auto e : bndel) {  // for elements contributing to chare boundary
        const auto [ N, grad, u, J ] = egrad( e, coord, inpoel, U );
        auto J24 = J/24.0;
        for (std::size_t a=0; a<4; ++a) {
          auto i = bid.find( gid[N[a]] );
          if (i != end(bid))    // only contribute to chare-boundary node
            for (std::size_t b=0; b<4; ++b)
              for (std::size_t j=0; j<3; ++j)
                for (std::size_t c=0; c<m_ncomp; ++c)
                  G(i->second,c*3+j,0) += J24 * grad[b][j] * u[c][b];
        }
      }
    }

    //! Compute right hand side for ALECG
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] esued Elements surrounding edges
    //! \param[in] triinpoel Boundary triangle face connecitivity
    //! \param[in] gid Local->global node id map
    //! \param[in] bid Local chare-boundary node ids (value) associated to
    //!    global node ids (key)
    //! \param[in] lid Global->local node ids
    //! \param[in] dfnorm Dual-face normals associated to edges
    //! \param[in] bnorm Face normals in boundary points
    //! \param[in] vol Nodal volumes
    //! \param[in] G Nodal gradients
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real /* t */,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const std::unordered_map< tk::UnsMesh::Edge,
                      std::vector< std::size_t >, tk::UnsMesh::Hash<2>,
                      tk::UnsMesh::Eq<2> >& esued,
              const std::pair< std::vector< std::size_t >,
                               std::vector< std::size_t > >& psup,
              const std::vector< std::size_t >& triinpoel,
              const std::vector< std::size_t >& gid,
              const std::unordered_map< std::size_t, std::size_t >& bid,
              const std::unordered_map< std::size_t, std::size_t >& lid,
              const std::unordered_map< tk::UnsMesh::Edge,
                      std::array< tk::real, 3 >,
                      tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> >& dfnorm,
              const std::unordered_map< std::size_t,
                      std::array< tk::real, 4 > >& bnorm,
              const std::vector< tk::real >& vol,
              const tk::Fields& G,
              const tk::Fields& U,
              tk::Fields& R ) const
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
      for (ncomp_t c=0; c<5; ++c) R.fill( c, m_offset, 0.0 );

      // access pointer to right hand side at component and offset
      std::array< const tk::real*, 5 > r;
      for (ncomp_t c=0; c<5; ++c) r[c] = R.cptr( c, m_offset );
      
      // allocate storage for nodal gradients of primitive variables
      tk::Fields Grad( U.nunk(), m_ncomp*3 );
      Grad.fill( 0.0 );

      // copy in nodal gradients of chare-boundary points
      for (const auto& [g,b] : bid) {
        auto i = tk::cref_find( lid, g );
        for (ncomp_t c=0; c<Grad.nprop(); ++c)
          Grad(i,c,0) = G(b,c,0);
      }
      
      //// for verification only, will go away once correct
      //tk::Fields V( U.nunk(), 3 );
      //V.fill( 0.0 );

      // compute gradients of primitive variables in internal points
      for (std::size_t e=0; e<inpoel.size()/4; ++e) {
        const auto [ N, grad, u, J ] = egrad( e, coord, inpoel, U );
        auto J24 = J/24.0;
        for (std::size_t a=0; a<4; ++a) {
          auto i = bid.find( gid[N[a]] );
          if (i == end(bid))    // only contribute to internal nodes
            for (std::size_t b=0; b<4; ++b)
              for (std::size_t j=0; j<3; ++j)
                for (std::size_t c=0; c<m_ncomp; ++c)
                  Grad(N[a],c*3+j,0) += J24 * grad[b][j] * u[c][b];
        }
      }

      // divide weak result in gradients by nodal volume
      for (std::size_t p=0; p<Grad.nunk(); ++p)
        for (std::size_t c=0; c<Grad.nprop(); ++c)
           Grad(p,c,0) /= vol[p];

      // domain-edge integral
      for (std::size_t p=0; p<U.nunk(); ++p) {  // for each point p
        for (auto q : tk::Around(psup,p)) {     // for each edge p-q
          // access and orient dual-face normals for edge p-q
          auto n = tk::cref_find( dfnorm, {gid[p],gid[q]} );
          if (gid[p] > gid[q]) { n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2]; }
          // Access primitive variables at edge-end points
          std::array< std::vector< tk::real >, 2 >
            ru{ std::vector<tk::real>(5,0.0), std::vector<tk::real>(5,0.0) };
#if 0
          for (std::size_t c=0; c<5; ++c) {
            ru[0][c] = U(p, c, m_offset);
            ru[1][c] = U(q, c, m_offset);
          }
#else
          // density
          ru[0][0] = U(p, 0, m_offset);
          ru[1][0] = U(q, 0, m_offset);
          // divide out density
          for (std::size_t c=1; c<5; ++c) {
            ru[0][c] =  U(p, c, m_offset)/ ru[0][0];
            ru[1][c] =  U(q, c, m_offset)/ ru[1][0];
          }
          // convert to internal energy
          for (std::size_t d=0; d<3; ++d) {
            ru[0][4] -= 0.5*ru[0][1+d]*ru[0][1+d];
            ru[1][4] -= 0.5*ru[1][1+d]*ru[1][1+d];
          }
          // compute MUSCL reconstruction in edge-end points
          tk::muscl( {p,q}, coord, Grad, ru );
          // convert back to conserved
          // convert to internal energy
          for (std::size_t d=0; d<3; ++d) {
            ru[0][4] += 0.5*ru[0][1+d]*ru[0][1+d];
            ru[1][4] += 0.5*ru[1][1+d]*ru[1][1+d];
          }
          // multiply density
          for (std::size_t c=1; c<5; ++c) {
            ru[0][c] *= ru[0][0];
            ru[1][c] *= ru[1][0];
         }
#endif
            
          // edge vector = outward face normal of the dual mesh face
          std::array< tk::real, 3 > fn{ 0., 0., 0. };

          // compute domain integral
          for (auto e : tk::cref_find(esued,{p,q})) {
            // access node IDs
            const std::array< std::size_t, 4 >
                N{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
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
            // sum flux contributions to nodes
            auto J24 = J/24.0;
            for (const auto& [a,b] : tk::lpoed) {
              auto s = tk::orient( {N[a],N[b]}, {p,q} );
              for (std::size_t j=0; j<3; ++j) {
                fn[j] += J24 * s * (grad[a][j] - grad[b][j]);
                //V(p,j,0) -= 2.0*J48 * s * (grad[a][j] - grad[b][j]);
              }
            }
          }
          tk::real A = 0;
          for (std::size_t i=0; i<3; ++i) A += fn[i] * fn[i];
          A = std::sqrt(A);
          for (std::size_t i=0; i<3; ++i) fn[i] /= A;
              
          // Compute Riemann flux using edge-end point states
          //auto f = HLLC::flux( fn, ru, {{0.0,0.0,0.0}} );
          auto f = rusanov_flux( fn, ru );
          for (std::size_t c=0; c<m_ncomp; ++c) R.var(r[c],p) -= A*f[c];
                
        }
      }
      
      //// Optional source
      //for (std::size_t p=0; p<U.nunk(); ++p) {
      //  auto s = Problem::src( m_system, m_ncomp, x[p], y[p], z[p], t );
      //  for (std::size_t c=0; c<5; ++c) R.var(r[c],p) += vol[p] * s[c];
      //}

      // Boundary integrals
      for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
        // access node IDs
        const std::array< std::size_t, 3 >
          N{ triinpoel[e*3+0], triinpoel[e*3+1], triinpoel[e*3+2] };
        // node coordinates
        std::array< tk::real, 3 > xp{ x[N[0]], x[N[1]], x[N[2]] },
                                  yp{ y[N[0]], y[N[1]], y[N[2]] },
                                  zp{ z[N[0]], z[N[1]], z[N[2]] };
        // compute face area
        auto A = tk::area( xp, yp, zp );
        auto A24 = A / 24.0;
        auto A6 = A / 6.;
        // compute face normal
        auto n = tk::normal( xp, yp, zp );
        // access solution at element nodes
        std::vector< std::array< tk::real, 3 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, m_offset, N );
        // compute fluxes
        std::array< std::array< tk::real, 3 >, 5 > f;
        if ( bnorm.find(N[0]) == bnorm.end() )
          flux( n, u, f );
        else {
          wall_flux(n, u, f);
        }

        // sum boundary integral contributions to boundary nodes
        for (const auto& [a,b] : tk::lpoet) {
          for (std::size_t c=0; c<m_ncomp; ++c) {
            auto Bab = A24 * (f[c][a] + f[c][b]);
            R.var(r[c],N[a]) -= Bab + A6 * f[c][a];
            R.var(r[c],N[b]) -= Bab;
          }
          //for (std::size_t j=0; j<3; ++j ) {
          //  V(N[a],j,0) -= 2.0 * A24 * n[j];
          //  V(N[b],j,0) -= 2.0 * A24 * n[j];
          //  V(N[a],j,0) -= A6 * n[j];
          //}

        }
      }
      
      // test 2*sum_{vw in v} D_i^{vw} + 2*sum_{vw in v} B_i^{vw} + B_i^v = 0
      // for boundary points (this only makes sense in serial)
      // for (std::size_t p=0; p<coord[0].size(); ++p)
      //     for (std::size_t j=0; j<3; ++j)
      //       if (std::abs(V(p,j,0)) > 1.0e-15)
      //         std::cout << 'b' << std::endl;

    }

    static void flux(
        const std::array< tk::real, 3 >& fn,
        const std::vector< std::array< tk::real, 3 > >& u,
        std::array< std::array< tk::real, 3 >, 5 >& f)
    {

      for (std::size_t i=0; i<3; ++i) {
        auto dinv = 1 / u[0][i];
        auto p = eos_pressure< tag::compflow >(
            0, u[0][i], u[1][i]*dinv, u[2][i]*dinv, u[3][i]*dinv, u[4][i] );
      
        tk::real div = 0;
        for (std::size_t d=0; d<3; ++d) div += fn[d] * u[1+d][i];
        div *= dinv;

        f[0][i] = u[0][i] * div;
        for (std::size_t d=0; d<3; ++d) f[1+d][i] = u[1+d][i]*div + p*fn[d];
        f[4][i] = (u[4][i] + p) * div;
      }

    }
    static void wall_flux(
        const std::array< tk::real, 3 >& fn,
        const std::vector< std::array< tk::real, 3 > >& u,
        std::array< std::array< tk::real, 3 >, 5 >& f)
    {

      for (std::size_t i=0; i<3; ++i) {
        auto dinv = 1 / u[0][i];
        auto p = eos_pressure< tag::compflow >(
            0, u[0][i], u[1][i]*dinv, u[2][i]*dinv, u[3][i]*dinv, u[4][i] );
      
        f[0][i] = 0;
        for (std::size_t d=0; d<3; ++d) f[1+d][i] = p*fn[d];
        f[4][i] = 0;
      }

    }
  
    static auto rusanov_flux( const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u )
    {
      std::vector< tk::real > flx( u[0].size(), 0 );

      const auto & UL = u[0];
      const auto & UR = u[1];

      // Primitive variables
      auto rhol = UL[0];
      auto rhor = UR[0];

      auto ul = UL[1]/rhol;
      auto vl = UL[2]/rhol;
      auto wl = UL[3]/rhol;

      auto ur = UR[1]/rhor;
      auto vr = UR[2]/rhor;
      auto wr = UR[3]/rhor;

      auto pl = eos_pressure< tag::compflow >( 0, rhol, ul, vl, wl, UL[4] );
      auto pr = eos_pressure< tag::compflow >( 0, rhor, ur, vr, wr, UR[4] );

      auto al = eos_soundspeed< tag::compflow >( 0, rhol, pl );
      auto ar = eos_soundspeed< tag::compflow >( 0, rhor, pr );

      // Face-normal velocities
      tk::real vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];
      tk::real vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];


      // Numerical fluxes
      flx[0]  = 0.5 * ( UL[0] * vnl );
      flx[1]  = 0.5 * ( UL[1] * vnl + pl*fn[0] );
      flx[2]  = 0.5 * ( UL[2] * vnl + pl*fn[1] );
      flx[3]  = 0.5 * ( UL[3] * vnl + pl*fn[2] );
      flx[4]  = 0.5 * ( ( UL[4] + pl ) * vnl );
      
      flx[0] += 0.5 * ( UR[0] * vnr );
      flx[1] += 0.5 * ( UR[1] * vnr + pr*fn[0] );
      flx[2] += 0.5 * ( UR[2] * vnr + pr*fn[1] );
      flx[3] += 0.5 * ( UR[3] * vnr + pr*fn[2] );
      flx[4] += 0.5 * ( ( UR[4] + pr ) * vnr );
      
      auto sl = std::abs(vnl) + al;
      auto sr = std::abs(vnr) + ar;
      auto smax = std::max( sl, sr );

      flx[0] -= 0.5 * smax * ( UR[0] - UL[0] ); 
      flx[1] -= 0.5 * smax * ( UR[1] - UL[1] ); 
      flx[2] -= 0.5 * smax * ( UR[2] - UL[2] ); 
      flx[3] -= 0.5 * smax * ( UR[3] - UL[3] ); 
      flx[4] -= 0.5 * smax * ( UR[4] - UL[4] ); 

      return flx;
    }


    //! Compute right hand side for DiagCG (CG-FCT)
    //! \param[in] t Physical time
    //! \param[in] deltat Size of time step
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] Ue Element-centered solution vector at intermediate step
    //!    (used here internally as a scratch array)
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
      Assert( R.nunk() == coord[0].size(),
              "Number of unknowns and/or number of components in right-hand "
              "side vector incorrect" );

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

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
          p[a] = eos_pressure< tag::compflow >
                   ( m_system, u[0][a], u[1][a]/u[0][a], u[2][a]/u[0][a],
                     u[3][a]/u[0][a], u[4][a] );

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
        std::array< std::vector< tk::real >, 4 > s{{
          Problem::src( m_system, m_ncomp, x[N[0]], y[N[0]], z[N[0]], t ),
          Problem::src( m_system, m_ncomp, x[N[1]], y[N[1]], z[N[1]], t ),
          Problem::src( m_system, m_ncomp, x[N[2]], y[N[2]], z[N[2]], t ),
          Problem::src( m_system, m_ncomp, x[N[3]], y[N[3]], z[N[3]], t ) }};
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
        auto p = eos_pressure< tag::compflow >
                   ( m_system, ue[0], ue[1]/ue[0], ue[2]/ue[0], ue[3]/ue[0],
                     ue[4] );

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
        auto s = Problem::src( m_system, m_ncomp, xc, yc, zc, t+deltat/2 );
        for (std::size_t c=0; c<5; ++c)
          for (std::size_t a=0; a<4; ++a)
            R.var(r[c],N[a]) += d/4.0 * s[c];

      }
//         // add viscous stress contribution to momentum and energy rhs
//         m_physics.viscousRhs( deltat, J, N, grad, u, r, R );
//         // add heat conduction contribution to energy rhs
//         m_physics.conductRhs( deltat, J, N, grad, u, r, R );
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
      auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0][0];
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
          auto p = eos_pressure< tag::compflow >
                     ( m_system, r, ru/r, rv/r, rw/r, re );
          if (p < 0) p = 0.0;
          auto c = eos_soundspeed< tag::compflow >( m_system, r, p );
          auto v = std::sqrt((ru*ru + rv*rv + rw*rw)/r/r) + c; // char. velocity
          if (v > maxvel) maxvel = v;
        }
        // compute element dt for the Euler equations
        auto euler_dt = L / maxvel;
        // compute element dt based on the viscous force
        auto viscous_dt = m_physics.viscous_dt( L, u );
        // compute element dt based on thermal diffusion
        auto conduct_dt = m_physics.conduct_dt( L, g, u );
        // compute minimum element dt
        auto elemdt = std::min( euler_dt, std::min( viscous_dt, conduct_dt ) );
        // find minimum dt across all elements
        if (elemdt < mindt) mindt = elemdt;
      }
      return mindt;
    }

    //! Extract the velocity field at cell nodes. Currently unused.
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
    { m_problem.side( conf ); }

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
    //! \param[in] t Physical time
    //! \param[in] deltat Time step size
    //! \param[in] ss Pair of side set ID and (local) node IDs on the side set
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
      using tag::param; using tag::compflow; using tag::bcdir;
      using NodeBC = std::vector< std::pair< bool, tk::real > >;
      std::map< std::size_t, NodeBC > bc;
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
              auto s = m_problem.solinc( m_system, m_ncomp, x[n], y[n], z[n],
                                         t, deltat );
              bc[n] = {{ {true,s[0]}, {true,s[1]}, {true,s[2]}, {true,s[3]},
                         {true,s[4]} }};
            }
      }
      return bc;
    }

    //! Set symmetry boundary conditions at nodes
    //! \param[in] U Solution vector at recent time step
    //! \param[in] bnorm Face normals in boundary points: key local node id,
    //!    value: unit normal
    void
    symbc( tk::Fields& U,
           const std::unordered_map<std::size_t,std::array<tk::real,4>>& bnorm )
    const {
      for (const auto& [ i, nr ] : bnorm ) {
        std::array< tk::real, 3 >
          n{ nr[0], nr[1], nr[2] },
          v{ U(i,1,m_offset), U(i,2,m_offset), U(i,3,m_offset) };
        auto v_dot_n = tk::dot( v, n );
        U(i,1,m_offset) -= v_dot_n * n[0];
        U(i,2,m_offset) -= v_dot_n * n[1];
        U(i,3,m_offset) -= v_dot_n * n[2];
      }
    }

    //! Query nodes at which symmetry boundary conditions are set
    //! \param[in] bface Boundary-faces mapped to side set ids
    //! \param[in] triinpoel Boundary-face connectivity
    //! \param[in,out] nodes Node ids at which symmetry BCs are set
    void
    symbcnodes( const std::map< int, std::vector< std::size_t > >& bface,
                const std::vector< std::size_t >& triinpoel,
                std::unordered_set< std::size_t >& nodes ) const
    {
      using tag::param; using tag::compflow; using tag::bcsym;
      const auto& bc = g_inputdeck.get< param, compflow, bcsym >();
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
    std::vector< std::string > fieldNames() const
    { return m_problem.fieldNames( m_ncomp ); }

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
    {
      return
        m_problem.fieldOutput( m_system, m_ncomp, m_offset, t, V, v, coord, U );
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const
    { return m_problem.names( m_ncomp ); }

  private:
    const Physics m_physics;            //!< Physics policy
    const Problem m_problem;            //!< Problem policy
    const ncomp_t m_system;             //!< Equation system index
    const ncomp_t m_ncomp;              //!< Number of components in this PDE
    const ncomp_t m_offset;             //!< Offset PDE operates from
    
    //! Compute element contribution to nodal gradient
    //! \param[in] e Element whose contribution to compute
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] U Solution vector at recent time step
    //! \return Tuple of element contribution
    auto egrad( std::size_t e,
                const std::array< std::vector< tk::real >, 3 >& coord,
                const std::vector< std::size_t >& inpoel,
                const tk::Fields& U ) const
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
      std::vector< std::array< tk::real, 4 > > u( m_ncomp );
      for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, m_offset, N );
      // divide out density
      for (std::size_t c=1; c<5; ++c)
        for (std::size_t j=0; j<4; ++j )
          u[c][j] /= u[0][j];
      // convert to internal energy
      for (std::size_t d=0; d<3; ++d)
        for (std::size_t j=0; j<4; ++j )
          u[4][j] -= 0.5*u[1+d][j]*u[1+d][j];
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

#endif // CGCompFlow_h
