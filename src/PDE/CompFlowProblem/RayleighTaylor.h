// *****************************************************************************
/*!
  \file      src/PDE/CompFlowProblem/RayleighTaylor.h
  \author    F. Gonzalez
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem configurations for the compressible flow equations
  \details   This file defines a policy classe for the compressible flow
    equations, defined in PDE/CompFlow.h. See PDE/CompFlow.h for general
    requirements on flow equations problem policy classes.
*/
// *****************************************************************************
#ifndef CompFlowProblemRayleighTaylor_h
#define CompFlowProblemRayleighTaylor_h

#include <string>
#include <unordered_set>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! CompFlow system of PDEs problem: Rayleigh-Taylor
//! \see Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
//!   equations with relevance to Inertial Confinement Fusion", Journal of
//!   Computational Physics 267 (2014) 196-209.
class CompFlowProblemRayleighTaylor {
  public:

    //! Set initial conditions
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    static void init( const std::array< std::vector< tk::real >, 3 >& coord,
                      const std::vector< std::size_t >&,
                      tk::Fields& unk,
                      tk::ctr::ncomp_type e,
                      tk::ctr::ncomp_type offset,
                      tk::real )
    {
      Assert( coord[0].size() == unk.nunk(), "Size mismatch" );
      // manufactured solution parameters
      const auto& a =
        g_inputdeck.get< tag::param, tag::compflow, tag::alpha >()[e];
      const auto& bx =
        g_inputdeck.get< tag::param, tag::compflow, tag::betax >()[e];
      const auto& by =
        g_inputdeck.get< tag::param, tag::compflow, tag::betay >()[e];
      const auto& bz =
        g_inputdeck.get< tag::param, tag::compflow, tag::betaz >()[e];
      const auto& p0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::p0 >()[e];
      const auto& r0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::r0 >()[e];
      // ratio of specific heats
      tk::real g =
        g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[e];

      // set initial and boundary conditions
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
      for (ncomp_t i=0; i<x.size(); ++i) {
        auto& r  = unk(i,0,offset); // rho
        auto& ru = unk(i,1,offset); // rho * u
        auto& rv = unk(i,2,offset); // rho * v
        auto& rw = unk(i,3,offset); // rho * w
        auto& re = unk(i,4,offset); // rho * e
        // pressure field
        tk::real p = p0 + a*(bx*x[i]*x[i] + by*y[i]*y[i] + bz*z[i]*z[i]);
        r = r0-(bx*x[i]*x[i] + by*y[i]*y[i] + bz*z[i]*z[i]);
        ru = r*(z[i]*std::sin(M_PI*x[i]));
        rv = r*(z[i]*std::cos(M_PI*y[i]));
        rw = r*(-M_PI*z[i]*z[i]*(std::cos(M_PI*x[i])-std::sin(M_PI*y[i]))/2.0);
        re = p/(g-1) + (ru*ru + rv*rv + rw*rw)/2.0;
      }
    }

    //! Add source term to rhs for Rayleigh-Taylor manufactured solution
    //! \param[in] t Physical time
    //! \param[in] coord Mesh node coordinates
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] dt Size of time step
    //! \param[in] N Element node indices
    //! \param[in] mass Element mass matrix, nnode*nnode [4][4]
    //! \param[in] r Pointers to right hand side at component and offset
    //! \param[in,out] R Right-hand side vector contributing to
    static void
    sourceRhs( tk::real t,
               const std::array< std::vector< tk::real >, 3 >& coord,
               tk::ctr::ncomp_type e,
               tk::real dt,
               const std::array< std::size_t, 4 >& N,
               const std::array< std::array< tk::real, 4 >, 4 >& mass,
               const std::array< const tk::real*, 5 >& r,
               tk::Fields& R )
    {
      // manufactured solution parameters
      const auto& a =
        g_inputdeck.get< tag::param, tag::compflow, tag::alpha >()[e];
      const auto& bx =
        g_inputdeck.get< tag::param, tag::compflow, tag::betax >()[e];
      const auto& by =
        g_inputdeck.get< tag::param, tag::compflow, tag::betay >()[e];
      const auto& bz =
        g_inputdeck.get< tag::param, tag::compflow, tag::betaz >()[e];
      const auto& kappa =
        g_inputdeck.get< tag::param, tag::compflow, tag::kappa >()[e];
      const auto& p0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::p0 >()[e];
      const auto& r0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::r0 >()[e];
      // ratio of specific heats
      tk::real g =
        g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[e];

      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // temporal component of velocity field
      tk::real f = std::cos(kappa*M_PI*t);

      // derivative of temporal component of velocity field
      tk::real df = -kappa*M_PI*std::sin(kappa*M_PI*t);

      // spatial component of velocity field
      std::array< std::array< tk::real, 4 >, 3 >
        gx{{ {{ z[N[0]]*std::sin(M_PI*x[N[0]]),
                z[N[1]]*std::sin(M_PI*x[N[1]]),
                z[N[2]]*std::sin(M_PI*x[N[2]]),
                z[N[3]]*std::sin(M_PI*x[N[3]]) }},
             {{ z[N[0]]*std::cos(M_PI*y[N[0]]),
                z[N[1]]*std::cos(M_PI*y[N[1]]),
                z[N[2]]*std::cos(M_PI*y[N[2]]),
                z[N[3]]*std::cos(M_PI*y[N[3]]) }},
             {{ -M_PI*z[N[0]]*z[N[0]]*
                  (std::cos(M_PI*x[N[0]])-std::sin(M_PI*y[N[0]]))/2.0,
                -M_PI*z[N[1]]*z[N[1]]*
                  (std::cos(M_PI*x[N[1]])-std::sin(M_PI*y[N[1]]))/2.0,
                -M_PI*z[N[2]]*z[N[2]]*
                  (std::cos(M_PI*x[N[2]])-std::sin(M_PI*y[N[2]]))/2.0,
                -M_PI*z[N[3]]*z[N[3]]*
                  (std::cos(M_PI*x[N[3]])-std::sin(M_PI*y[N[3]]))/2.0 }} }};

      // gradients of spatial components of velocity field
      std::array< std::array< std::array< tk::real, 4 >, 3 >, 3 >
        dg{{ {{ {{ M_PI*z[N[0]]*std::cos(M_PI*x[N[0]]),
                   M_PI*z[N[1]]*std::cos(M_PI*x[N[1]]),
                   M_PI*z[N[2]]*std::cos(M_PI*x[N[2]]),
                   M_PI*z[N[3]]*std::cos(M_PI*x[N[3]]) }}, // dg1/dx1
                {{ 0.0, 0.0, 0.0, 0.0 }}, // dg1/dx2
                {{ std::sin(M_PI*x[N[0]]),
                   std::sin(M_PI*x[N[1]]),
                   std::sin(M_PI*x[N[2]]),
                   std::sin(M_PI*x[N[3]]) }} }}, // dg1/dx3
             {{ {{ 0.0, 0.0, 0.0, 0.0 }}, // dg2/dx1
                {{ -M_PI*z[N[0]]*std::sin(M_PI*y[N[0]]),
                   -M_PI*z[N[1]]*std::sin(M_PI*y[N[1]]),
                   -M_PI*z[N[2]]*std::sin(M_PI*y[N[2]]),
                   -M_PI*z[N[3]]*std::sin(M_PI*y[N[3]]) }}, // dg2/dx2
                {{ std::cos(M_PI*y[N[0]]),
                   std::cos(M_PI*y[N[1]]),
                   std::cos(M_PI*y[N[2]]),
                   std::cos(M_PI*y[N[3]]) }} }}, // dg2/dx3
             {{ {{ M_PI*M_PI*z[N[0]]*z[N[0]]*std::sin(M_PI*x[N[0]])/2.0,
                   M_PI*M_PI*z[N[1]]*z[N[1]]*std::sin(M_PI*x[N[1]])/2.0,
                   M_PI*M_PI*z[N[2]]*z[N[2]]*std::sin(M_PI*x[N[2]])/2.0,
                   M_PI*M_PI*z[N[3]]*z[N[3]]*std::sin(M_PI*x[N[3]])/2.0 }},
                {{ M_PI*M_PI*z[N[0]]*z[N[0]]*std::cos(M_PI*y[N[0]])/2.0,
                   M_PI*M_PI*z[N[1]]*z[N[1]]*std::cos(M_PI*y[N[1]])/2.0,
                   M_PI*M_PI*z[N[2]]*z[N[2]]*std::cos(M_PI*y[N[2]])/2.0,
                   M_PI*M_PI*z[N[3]]*z[N[3]]*std::cos(M_PI*y[N[3]])/2.0 }},
                {{ -M_PI*z[N[0]]*
                    (std::cos(M_PI*x[N[0]])-std::sin(M_PI*y[N[0]])),
                   -M_PI*z[N[1]]*
                    (std::cos(M_PI*x[N[1]])-std::sin(M_PI*y[N[1]])),
                   -M_PI*z[N[2]]*
                    (std::cos(M_PI*x[N[2]])-std::sin(M_PI*y[N[2]])),
                   -M_PI*z[N[3]]*
                    (std::cos(M_PI*x[N[3]])-std::sin(M_PI*y[N[3]])) }} }} }};

      // density field (NOTE: need to change variable name)
      std::array< tk::real, 4 > rho{{
        r0 - (bx*x[N[0]]*x[N[0]] + by*y[N[0]]*y[N[0]] + bz*z[N[0]]*z[N[0]]),
        r0 - (bx*x[N[1]]*x[N[1]] + by*y[N[1]]*y[N[1]] + bz*z[N[1]]*z[N[1]]),
        r0 - (bx*x[N[2]]*x[N[2]] + by*y[N[2]]*y[N[2]] + bz*z[N[2]]*z[N[2]]),
        r0 - (bx*x[N[3]]*x[N[3]] + by*y[N[3]]*y[N[3]] + bz*z[N[3]]*z[N[3]]) }};

      // density gradient
      std::array< std::array< tk::real, 4 >, 3 >
        dr{{ {{ -2.0*bx*x[N[0]],
                -2.0*bx*x[N[1]],
                -2.0*bx*x[N[2]],
                -2.0*bx*x[N[3]] }},
             {{ -2.0*by*y[N[0]],
                -2.0*by*y[N[1]],
                -2.0*by*y[N[2]],
                -2.0*by*y[N[3]] }},
             {{ -2.0*bz*z[N[0]],
                -2.0*bz*z[N[1]],
                -2.0*bz*z[N[2]],
                -2.0*bz*z[N[3]] }} }};

      // pressure field
      std::array< tk::real, 4 > p{{
        p0 + a*(bx*x[N[0]]*x[N[0]] + by*y[N[0]]*y[N[0]] + bz*z[N[0]]*z[N[0]]),
        p0 + a*(bx*x[N[1]]*x[N[1]] + by*y[N[1]]*y[N[1]] + bz*z[N[1]]*z[N[1]]),
        p0 + a*(bx*x[N[2]]*x[N[2]] + by*y[N[2]]*y[N[2]] + bz*z[N[2]]*z[N[2]]),
        p0 + a*(bx*x[N[3]]*x[N[3]] + by*y[N[3]]*y[N[3]] + bz*z[N[3]]*z[N[3]])}};

      // pressure gradient
      std::array< std::array< tk::real, 4 >, 3 >
        dp{{ {{ 2.0*a*bx*x[N[0]],
                2.0*a*bx*x[N[1]],
                2.0*a*bx*x[N[2]],
                2.0*a*bx*x[N[3]] }},
             {{ 2.0*a*by*y[N[0]],
                2.0*a*by*y[N[1]],
                2.0*a*by*y[N[2]],
                2.0*a*by*y[N[3]] }},
             {{ 2.0*a*bz*z[N[0]],
                2.0*a*bz*z[N[1]],
                2.0*a*bz*z[N[2]],
                2.0*a*bz*z[N[3]] }} }};

      // density source
      std::array< tk::real, 4 > Sr{{
        f*(gx[0][0]*dr[0][0] + gx[1][0]*dr[1][0] + gx[2][0]*dr[2][0]),
        f*(gx[0][1]*dr[0][1] + gx[1][1]*dr[1][1] + gx[2][1]*dr[2][1]),
        f*(gx[0][2]*dr[0][2] + gx[1][2]*dr[1][2] + gx[2][2]*dr[2][2]),
        f*(gx[0][3]*dr[0][3] + gx[1][3]*dr[1][3] + gx[2][3]*dr[2][3]) }};

      // specific total energy
      std::array< tk::real, 4 > E{{
        p[0]/(rho[0]*(g-1.0)) + 0.5*f*f*(gx[0][0]*gx[0][0] + gx[1][0]*gx[1][0] +
          gx[2][0]*gx[2][0]),
        p[1]/(rho[1]*(g-1.0)) + 0.5*f*f*(gx[0][1]*gx[0][1] + gx[1][1]*gx[1][1] +
          gx[2][1]*gx[2][1]),
        p[2]/(rho[2]*(g-1.0)) + 0.5*f*f*(gx[0][2]*gx[0][2] + gx[1][2]*gx[1][2] +
          gx[2][2]*gx[2][2]),
        p[3]/(rho[3]*(g-1.0)) + 0.5*f*f*(gx[0][3]*gx[0][3] + gx[1][3]*gx[1][3] +
          gx[2][3]*gx[2][3]) }};

      // temporal derivatives of specific total energy
      std::array< tk::real, 4 > dEt{{
        f*df*(gx[0][0]*gx[0][0] + gx[1][0]*gx[1][0] + gx[2][0]*gx[2][0]),
        f*df*(gx[0][1]*gx[0][1] + gx[1][1]*gx[1][1] + gx[2][1]*gx[2][1]),
        f*df*(gx[0][2]*gx[0][2] + gx[1][2]*gx[1][2] + gx[2][2]*gx[2][2]),
        f*df*(gx[0][3]*gx[0][3] + gx[1][3]*gx[1][3] + gx[2][3]*gx[2][3]) }};

      // spatial derivatives of specific total energy
      std::array< std::array< tk::real, 4 >, 3 >
        dE{{ {{ dp[0][0]/(rho[0]*(g-1)) - (p[0]*dr[0][0])/(rho[0]*rho[0]*(g-1))
                  + f*f*(gx[0][0]*dg[0][0][0] + gx[2][0]*dg[2][0][0]),
                dp[0][1]/(rho[1]*(g-1)) - (p[1]*dr[0][1])/(rho[1]*rho[1]*(g-1))
                  + f*f*(gx[0][1]*dg[0][0][1] + gx[2][1]*dg[2][0][1]),
                dp[0][2]/(rho[2]*(g-1)) - (p[2]*dr[0][2])/(rho[2]*rho[2]*(g-1))
                  + f*f*(gx[0][2]*dg[0][0][2] + gx[2][2]*dg[2][0][2]),
                dp[0][3]/(rho[3]*(g-1)) - (p[3]*dr[0][3])/(rho[3]*rho[3]*(g-1))
                  + f*f*(gx[0][3]*dg[0][0][3] + gx[2][3]*dg[2][0][3]) }},
             {{ dp[1][0]/(rho[0]*(g-1)) - (p[0]*dr[1][0])/(rho[0]*rho[0]*(g-1))
                  + f*f*(gx[1][0]*dg[1][1][0] + gx[2][0]*dg[2][1][0]),
                dp[1][1]/(rho[1]*(g-1)) - (p[1]*dr[1][1])/(rho[1]*rho[1]*(g-1))
                  + f*f*(gx[1][1]*dg[1][1][1] + gx[2][1]*dg[2][1][1]),
                dp[1][2]/(rho[2]*(g-1)) - (p[2]*dr[1][2])/(rho[2]*rho[2]*(g-1))
                  + f*f*(gx[1][2]*dg[1][1][2] + gx[2][2]*dg[2][1][2]),
                dp[1][3]/(rho[3]*(g-1)) - (p[3]*dr[1][3])/(rho[3]*rho[3]*(g-1))
                  + f*f*(gx[1][3]*dg[1][1][3] + gx[2][3]*dg[2][1][3]) }},
             {{ dp[2][0]/(rho[0]*(g-1)) - (p[0]*dr[2][0])/(rho[0]*rho[0]*(g-1))
                  + f*f*(gx[0][0]*dg[0][2][0] + gx[1][0]*dg[1][2][0] +
                  gx[2][0]*dg[2][2][0]),
                dp[2][1]/(rho[1]*(g-1)) - (p[1]*dr[2][1])/(rho[1]*rho[1]*(g-1))
                  + f*f*(gx[0][1]*dg[0][2][1] + gx[1][1]*dg[1][2][1] +
                  gx[2][1]*dg[2][2][1]),
                dp[2][2]/(rho[2]*(g-1)) - (p[2]*dr[2][2])/(rho[2]*rho[2]*(g-1))
                  + f*f*(gx[0][2]*dg[0][2][2] + gx[1][2]*dg[1][2][2] +
                  gx[2][2]*dg[2][2][2]),
                dp[2][3]/(rho[3]*(g-1)) - (p[3]*dr[2][3])/(rho[3]*rho[3]*(g-1))
                  + f*f*(gx[0][3]*dg[0][2][3] + gx[1][3]*dg[1][2][3] +
                  gx[2][3]*dg[2][2][3]) }} }};

      // energy source
      std::array< tk::real, 4 > Se{{
        rho[0]*dEt[0] + E[0]*Sr[0] +
          rho[0]*f*(gx[0][0]*dE[0][0] + gx[1][0]*dE[1][0] + gx[2][0]*dE[2][0]) +
          f*(gx[0][0]*dp[0][0] + gx[1][0]*dp[1][0] + gx[2][0]*dp[2][0]),
        rho[1]*dEt[1] + E[1]*Sr[1] +
          rho[1]*f*(gx[0][1]*dE[0][1] + gx[1][1]*dE[1][1] + gx[2][1]*dE[2][1]) +
          f*(gx[0][1]*dp[0][1] + gx[1][1]*dp[1][1] + gx[2][1]*dp[2][1]),
        rho[2]*dEt[2] + E[2]*Sr[2] +
          rho[2]*f*(gx[0][2]*dE[0][2] + gx[1][2]*dE[1][2] + gx[2][2]*dE[2][2]) +
          f*(gx[0][2]*dp[0][2] + gx[1][2]*dp[1][2] + gx[2][2]*dp[2][2]),
        rho[3]*dEt[3] + E[3]*Sr[3] +
          rho[3]*f*(gx[0][3]*dE[0][3] + gx[1][3]*dE[1][3] + gx[2][3]*dE[2][3]) +
          f*(gx[0][3]*dp[0][3] + gx[1][3]*dp[1][3] + gx[2][3]*dp[2][3]) }};

      // momentum source
      std::array< std::array< tk::real, 4 >, 3 >
        Sm{{ {{ rho[0]*gx[0][0]*df + f*gx[0][0]*Sr[0] +
                  rho[0]*f*f*(gx[0][0]*dg[0][0][0] + gx[2][0]*dg[0][2][0]),
                rho[1]*gx[0][1]*df + f*gx[0][1]*Sr[1] +
                  rho[1]*f*f*(gx[0][1]*dg[0][0][1] + gx[2][1]*dg[0][2][1]),
                rho[2]*gx[0][2]*df + f*gx[0][2]*Sr[2] +
                  rho[2]*f*f*(gx[0][2]*dg[0][0][2] + gx[2][2]*dg[0][2][2]),
                rho[3]*gx[0][3]*df + f*gx[0][3]*Sr[3] +
                  rho[3]*f*f*(gx[0][3]*dg[0][0][3] + gx[2][3]*dg[0][2][3]) }},
             {{ rho[0]*gx[1][0]*df + f*gx[1][0]*Sr[0] +
                  rho[0]*f*f*(gx[1][0]*dg[1][1][0] + gx[2][0]*dg[1][2][0] ),
                rho[1]*gx[1][1]*df + f*gx[1][1]*Sr[1] +
                  rho[1]*f*f*(gx[1][1]*dg[1][1][1] + gx[2][1]*dg[1][2][1] ),
                rho[2]*gx[1][2]*df + f*gx[1][2]*Sr[2] +
                  rho[2]*f*f*(gx[1][2]*dg[1][1][2] + gx[2][2]*dg[1][2][2] ),
                rho[3]*gx[1][3]*df + f*gx[1][3]*Sr[3] +
                  rho[3]*f*f*(gx[1][3]*dg[1][1][3] + gx[2][3]*dg[1][2][3] ) }},
             {{ rho[0]*gx[2][0]*df + f*gx[2][0]*Sr[0] +
                  rho[0]*f*f*(gx[0][0]*dg[2][0][0] + gx[1][0]*dg[2][1][0] +
                  gx[2][0]*dg[2][2][0]),
                rho[1]*gx[2][1]*df + f*gx[2][1]*Sr[1] +
                  rho[1]*f*f*(gx[0][1]*dg[2][0][1] + gx[1][1]*dg[2][1][1] +
                  gx[2][1]*dg[2][2][1]),
                rho[2]*gx[2][2]*df + f*gx[2][2]*Sr[2] +
                  rho[2]*f*f*(gx[0][2]*dg[2][0][2] + gx[1][2]*dg[2][1][2] +
                  gx[2][2]*dg[2][2][2]),
                rho[3]*gx[2][3]*df + f*gx[2][3]*Sr[3] +
                  rho[3]*f*f*(gx[0][3]*dg[2][0][3] + gx[1][3]*dg[2][1][3] +
                  gx[2][3]*dg[2][2][3]) }} }};

      // add momentum and energy source at element nodes
      for (std::size_t j=0; j<4; ++j)
        for (std::size_t k=0; k<4; ++k) {
          // source contribution to mass rhs
          R.var(r[0],N[j]) += dt * mass[j][k] * Sr[k];
          // source contribution to momentum rhs
          for (std::size_t l=0; l<3; ++l)
            R.var(r[l+1],N[j]) += dt * mass[j][k] * Sm[l][k];
          // source contribution to enerhy rhs
          R.var(r[4],N[j]) += dt * mass[j][k] * Se[k];
        }
    }

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    static void side( std::unordered_set< int >& conf ) {
      using tag::param; using tag::compflow; using tag::bcdir;
      for (const auto& s : g_inputdeck.get< param, compflow, bcdir >())
        for (const auto& i : s)
          conf.insert( std::stoi(i) );
    }

//     //! \brief Query Dirichlet boundary condition value on a given side set for
//     //!    all components in this PDE system
//     //! \param[in] sideset Side set ID
//     //! \return Vector of pairs of bool and BC value for all components
//     static std::vector< std::pair< bool, tk::real > > dirbc( int sideset ) {
//       using tag::param; using tag::compflow; using tag::bcdir;
//       std::vector< std::pair< bool, tk::real > > bc( 5, { false, 0.0 } );
//       for (const auto& s : g_inputdeck.get< param, compflow, bcdir >())
//         for (const auto& i : s)
//           if (std::stoi(i) == sideset)
//             for (auto& b : bc)
//                b = { true, 0.0 };
//       return bc;
//     }

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
    //! \return Vector of pairs of bool and boundary condition value associated
    //!   to mesh node IDs at which Dirichlet boundary conditions are set. Note
    //!   that instead of the actual boundary condition value, we return the
    //!   increment between t+dt and t, since that is what the solution requires
    //!   as we solve for the soution increments and not the solution itself.
    static std::unordered_map< std::size_t,
                               std::vector< std::pair< bool, tk::real > > >
    dirbc( tk::ctr::ncomp_type,
           tk::real,
           tk::real,
           const std::pair< const int, std::vector< std::size_t > >&,
           const std::array< std::vector< tk::real >, 3 >& )
    {
      using tag::param; using tag::compflow; using tag::bcdir;
      using NodeBC = std::vector< std::pair< bool, tk::real > >;
      std::unordered_map< std::size_t, NodeBC > bc;
      // TODO: include functionality from above dirbc() commented
      return bc;
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    static std::vector< std::string > fieldNames() {
      std::vector< std::string > n;
      n.push_back( "density numerical" );
      n.push_back( "density analytical" );
      n.push_back( "x-velocity numerical" );
      n.push_back( "x-velocity analytical" );
      n.push_back( "y-velocity numerical" );
      n.push_back( "y-velocity analytical" );
      n.push_back( "z-velocity numerical" );
      n.push_back( "z-velocity analytical" );
      n.push_back( "specific total energy numerical" );
      n.push_back( "specific total energy analytical" );
      return n;
    }

    //! Return field output going to file
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] t Physical time
    //! \param[in] coord Mesh node coordinates
    //! \param[in] U Solution vector at recent time step stage
    //! \return Vector of vectors to be output to file
    static std::vector< std::vector< tk::real > >
    fieldOutput( tk::ctr::ncomp_type e,
                 tk::ctr::ncomp_type offset,
                 tk::real t,
                 tk::real,
                 const std::vector< tk::real >&,
                 const std::array< std::vector< tk::real >, 3 >& coord,
                 const tk::Fields& U )
    {
      // manufactured solution parameters
      const auto& a =
        g_inputdeck.get< tag::param, tag::compflow, tag::alpha >()[e];
      const auto& bx =
        g_inputdeck.get< tag::param, tag::compflow, tag::betax >()[e];
      const auto& by =
        g_inputdeck.get< tag::param, tag::compflow, tag::betay >()[e];
      const auto& bz =
        g_inputdeck.get< tag::param, tag::compflow, tag::betaz >()[e];
      const auto& kappa =
        g_inputdeck.get< tag::param, tag::compflow, tag::kappa >()[e];
      const auto& p0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::p0 >()[e];
      const auto& r0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::r0 >()[e];
      // ratio of specific heats
      tk::real g =
        g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[e];

      std::vector< std::vector< tk::real > > out;
      const auto r  = U.extract( 0, offset );
      const auto ru = U.extract( 1, offset );
      const auto rv = U.extract( 2, offset );
      const auto rw = U.extract( 3, offset );
      const auto re = U.extract( 4, offset );

      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      std::vector< tk::real > rho = r;
      out.push_back( rho );
      for (std::size_t i=0; i<rho.size(); ++i) {
        rho[i] = r0-(bx*x[i]*x[i] + by*y[i]*y[i] + bz*z[i]*z[i]);
      }
      out.push_back( rho);

      std::vector< tk::real > u = ru;
      std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( u );
      for (std::size_t i=0; i<u.size(); ++i) {
        u[i] = rho[i]*std::cos(kappa*M_PI*t)*(z[i]*std::sin(M_PI*x[i]));
      }
      out.push_back( u );

      std::vector< tk::real > v = rv;
      std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( v );
      for (std::size_t i=0; i<v.size(); ++i) {
        v[i] = rho[i]*std::cos(kappa*M_PI*t)*(z[i]*std::cos(M_PI*y[i]));
      }
      out.push_back( v );

      std::vector< tk::real > w = rw;
      std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( w );
      for (std::size_t i=0; i<w.size(); ++i) {
        w[i] = rho[i]*std::cos(kappa*M_PI*t)*(-M_PI*z[i]*z[i]*
               (std::cos(M_PI*x[i])-std::sin(M_PI*y[i]))/2.0);
      }
      out.push_back( w );

      std::vector< tk::real > E = re;
      std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( E );
      for (std::size_t i=0; i<E.size(); ++i) {
        tk::real p = p0 + a*(bx*x[i]*x[i] + by*y[i]*y[i] + bz*z[i]*z[i]);
        E[i] = p/(g-1) + (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0;
      }
      out.push_back( E );

      return out;
   }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    static std::vector< std::string > names()
    { return { "r", "ru", "rv", "rw", "re" }; }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::RAYLEIGH_TAYLOR; }
};

} // inciter::

#endif // CompFlowProblemRayleighTaylor_h
