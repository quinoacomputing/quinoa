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
#include "Macro.hpp"
#include "Exception.hpp"
#include "Vector.hpp"
#include "EoS/EoS.hpp"
#include "Mesh/Around.hpp"

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

    //! Compute right hand side for ALECG
    //! \param[in] t Physical time
    //! \param[in] deltat Size of time step
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] psup Points surrounding points
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real t,
              tk::real deltat,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const std::pair< std::vector< std::size_t >,
                               std::vector< std::size_t > >& psup,
              const tk::Fields& U,
              tk::Fields& R ) const
    {
      std::cout << "CGCompflow::rhs" << std::endl;

      constexpr int num_vars = 7;
      
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      //------------------------------------------------------------------------
      // Update solution quntities
      //------------------------------------------------------------------------

      // get the number of nodes in the mesh
      const auto num_nodes = x.size();
      std::vector<tk::real> nodal_pressure(num_nodes);
      std::vector<tk::real> nodal_soundspeed(num_nodes);

      for ( std::size_t n=0; n<num_nodes; ++n ) {
        
        // access solution
        auto rho  = U( n, 0, m_offset ); 
        auto inv_rho = 1 / rho;
        auto velx = U( n, 1, m_offset ) * inv_rho; 
        auto vely = U( n, 2, m_offset ) * inv_rho; 
        auto velz = U( n, 3, m_offset ) * inv_rho; 
        auto ener = U( n, 4, m_offset );
          
        // get updated pressure and sound speed
        nodal_pressure[n] = eos_pressure< tag::compflow >
          ( m_system, rho, velx, vely, velz, ener );
        nodal_soundspeed[n] = eos_soundspeed< tag::compflow >( m_system, rho, nodal_pressure[n] );

      }

      //------------------------------------------------------------------------
      // Compute gradients
      //------------------------------------------------------------------------

      std::vector< tk::real > dudx_ave;
      dudx_ave.reserve( 3 * num_vars * num_nodes );

      for ( std::size_t n=0; n<num_nodes; ++n ) {

        // get current solution
        auto x0 = std::array<tk::real, 3>{ x[n], y[n], z[n] };

        std::array< tk::real, num_vars > u0;
        for (ncomp_t c=0; c<5; ++c) u0[c] = U( n, c, m_offset );
        u0[5] = nodal_pressure[n];
        u0[6] = nodal_soundspeed[n];

        auto umax = u0;
        auto umin = u0;

        // coefficients
        tk::real dxdx{0}, dxdy{0}, dxdz{0}, dydy{0}, dydz{0}, dzdz{0};

        std::array< tk::real, num_vars > dudx = {0, 0, 0, 0, 0, 0, 0};
        std::array< tk::real, num_vars > dudy = {0, 0, 0, 0, 0, 0, 0};
        std::array< tk::real, num_vars > dudz = {0, 0, 0, 0, 0, 0, 0};

        for (auto q : tk::Around(psup,n) ) {

          // get neighbor solution
          std::array< tk::real, num_vars > u1;
          for (ncomp_t c=0; c<5; ++c) u1[c] = U( q, c, m_offset );
          u1[5] = nodal_pressure[n];
          u1[6] = nodal_soundspeed[n];
          
          // keep track of max and min
          for ( std::size_t i=0; i<num_vars; ++i ) {
            umax[i] = std::max( umax[i], u1[i] );
            umin[i] = std::min( umin[i], u1[i] );
          }

          // deltas
          tk::real dx = x[q] - x0[0];
          tk::real dy = y[q] - x0[0];
          tk::real dz = z[q] - x0[0];

          // compute contributions
          dxdx += dx*dx;
          dxdy += dx*dy;
          dxdz += dx*dz;
          dydy += dy*dy;
          dydz += dy*dz;
          dzdz += dz*dz;

          for ( std::size_t i=0; i<num_vars; ++i ){
            auto du = u1[i] - u0[i];
            dudx[i] += du * dx;
            dudy[i] += du * dy;
            dudz[i] += du * dz;
          }

        } // neighbor

        // First compute minors
        auto min1 = dydy*dzdz-dydz*dydz;
        auto min2 = dxdy*dzdz-dydz*dxdz;
        auto min3 = dxdy*dydz-dydy*dxdz;
        // Now determinants
        auto denom = 1 / (dxdx*min1-dxdy*min2+dxdz*min3);

        for ( std::size_t i=0; i<num_vars; ++i ) {
          // component 1
          dudx_ave.emplace_back( (
              dudx[i] * min1
            - dxdy*(dudy[i]*dzdz-dydz*dudz[i])
            + dxdz*(dudy[i]*dydz-dydy*dudz[i])
          ) * denom );
          // component 2
          dudx_ave.emplace_back( (
              dxdx*(dudy[i]*dzdz-dydz*dudz[i])
            - dudx[i] * min2
            + dxdz*(dudz[i]*dxdy-dxdz*dudy[i])
          ) * denom );
          // component 3
          dudx_ave.emplace_back( (
              dxdx*(dudz[i]*dydy-dydz*dudy[i])
            - dxdy*(dudz[i]*dxdy-dxdz*dudy[i])
            + dudx[i] * min3
          ) * denom );
        }// vars


      } // points

 
      //------------------------------------------------------------------------
      // Compute Limited values
      //------------------------------------------------------------------------
      
      constexpr tk::real muscl_eps = 1.e-09;
      constexpr tk::real muscl_const = 0.3333333;
      constexpr auto muscl_con_m1 = muscl_const - 1;
      constexpr auto muscl_con_p1 = muscl_const + 1;
      
      // edges to node connectivity
      auto esup = tk::genEsup( inpoel, 4 );
      auto edge_points = tk::genInpoed(inpoel,4,esup);
      auto num_edge = edge_points.size() / 2;

      std::vector< tk::real > edge_unknowns( 2 * num_edge * num_vars );

      for ( std::size_t e=0; e<num_edge; ++e ) {
        auto p1 = edge_points[2*e];
        auto p2 = edge_points[2*e + 1];

        // edge vector, doubled since a factor of two is needed
        std::array< tk::real, 3 > l{ 2. * ( x[p2] - x[p1] ),
                                     2. * ( y[p2] - y[p1] ),
                                     2. * ( z[p2] - z[p1] ) };

        for ( std::size_t ivar = 0; ivar<num_vars; ++ivar ) {

          tk::real u1, u2;
          switch (ivar) {
            case (5):
              u1 = nodal_pressure[p1];
              u2 = nodal_pressure[p2];
              break;
            case (6):
              u1 = nodal_soundspeed[p1];
              u2 = nodal_soundspeed[p2];
              break;
            default:
              u1 = U( p1, ivar, m_offset );
              u2 = U( p2, ivar, m_offset );
          }

          // form delta_2 (u_j - u_i)
          auto delta_2 = u2 - u1;

          // form delta_1 (u_i - u_i-1) and delta_3 (u_j+1 - u_j)
          tk::real dot{0};
          for ( std::size_t i=0; i<3; ++i )
            dot += dudx_ave[p1*3*num_vars + ivar*3 + i] * l[i];
          auto delta_1 = dot - delta_2;
          
          dot = 0;
          for ( std::size_t i=0; i<3; ++i )
            dot += dudx_ave[p2*3*num_vars + ivar*3 + i] * l[i];
          auto delta_3 = dot - delta_2;

          // form limiters
          auto r_L = (delta_2 + muscl_eps) / (delta_1 + muscl_eps);
          auto r_R = (delta_2 + muscl_eps) / (delta_3 + muscl_eps);
          auto r_L_inv = (delta_1 + muscl_eps) / (delta_2 + muscl_eps);
          auto r_R_inv = (delta_3 + muscl_eps) / (delta_2 + muscl_eps);
          auto phi_L = (std::abs(r_L) + r_L) / (std::abs(r_L) + 1.);
          auto phi_R = (std::abs(r_R) + r_R) / (std::abs(r_R) + 1.);
          auto phi_L_inv = (std::abs(r_L_inv) + r_L_inv) / (std::abs(r_L_inv) + 1.);
          auto phi_R_inv = (std::abs(r_R_inv) + r_R_inv) / (std::abs(r_R_inv) + 1.);

          // final form of higher-order unknown
          edge_unknowns[e*2*num_vars + ivar] = u1 +
            0.25*(delta_1*muscl_con_m1*phi_L + delta_2*muscl_con_p1*phi_L_inv);

          edge_unknowns[e*2*num_vars + num_vars + ivar] = u2 -
            0.25*(delta_3*muscl_con_m1*phi_R + delta_2*muscl_con_p1*phi_R_inv);
          
        } // var

      }
      
      //------------------------------------------------------------------------
      // Compute some element geometry parameters
      //------------------------------------------------------------------------
      
      auto num_elem = inpoel.size()/4;
      std::vector< std::array< std::array< tk::real, 3 >, 4 > > elem_grad(num_elem);
      std::vector< tk::real > elem_volume(num_elem);

      for ( std::size_t e=0; e<num_elem; ++e ) {
        
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
        elem_volume[e] = J/6.;
        
        // shape function derivatives, nnode*ndim [4][3]
        auto & grad = elem_grad[e];
        grad[1] = tk::crossdiv( ca, da, J );
        grad[2] = tk::crossdiv( da, ba, J );
        grad[3] = tk::crossdiv( ba, ca, J );
        for (std::size_t i=0; i<3; ++i)
          grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];
      }

      //------------------------------------------------------------------------
      // Compute some edge geometry parameters
      //------------------------------------------------------------------------
     
      // elements surrounding edges
      auto esued = tk::genEsued( inpoel, 4, esup );

      //for ( std::size_t ed=0; ed<num_edge; ++ed ) {
      //  for ( auto el : tk::Around(esued, ed) ) {
      //  }
      //}

      //------------------------------------------------------------------------
      // Compute RHS
      //------------------------------------------------------------------------
      
      //for ( std::size_t e=0; e<num_edge; ++e ) {
      //  auto p1 = edge_points[2*e];
      //  auto p2 = edge_points[2*e + 1];
      //  auto ul = &edge_unknowns[e*2*num_vars];
      //  auto ur = &edge_unknowns[e*2*num_vars + num_vars];
      //  //auto f = flux( ul, ur, );
      //  
      //}

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
};

} // cg::

} // inciter::

#endif // CGCompFlow_h
