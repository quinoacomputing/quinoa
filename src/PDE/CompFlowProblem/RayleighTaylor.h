// *****************************************************************************
/*!
  \file      src/PDE/CompFlowProblem/RayleighTaylor.h
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

  private:
    //! Evaluate analytical solution at (x,y,z,t) for all components
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \param[in] t Time where to evaluate the solution
    //! \return Values of all components evaluated at (x,y,z,t)
    static std::array< tk::real, 5 >
    solution( tk::ctr::ncomp_type e,
              tk::real x, tk::real y, tk::real z, tk::real t )
   {
      using tag::param; using tag::compflow; using std::sin; using std::cos;
      // manufactured solution parameters
      const auto& a = g_inputdeck.get< param, compflow, tag::alpha >()[e];
      const auto& bx = g_inputdeck.get< param, compflow, tag::betax >()[e];
      const auto& by = g_inputdeck.get< param, compflow, tag::betay >()[e];
      const auto& bz = g_inputdeck.get< param, compflow, tag::betaz >()[e];
      const auto& p0 = g_inputdeck.get< param, compflow, tag::p0 >()[e];
      const auto& r0 = g_inputdeck.get< param, compflow, tag::r0 >()[e];
      const auto& k = g_inputdeck.get< param, compflow, tag::kappa >()[e];
      // ratio of specific heats
      const tk::real g = g_inputdeck.get< param, compflow, tag::gamma >()[e];
      // spatial component of density and pressure fields
      const tk::real gx = bx*x*x + by*y*y + bz*z*z;
      // density
      const tk::real r = r0 - gx;
      // pressure
      const tk::real p = p0 + a*gx;
      // velocity
      const tk::real ft = cos(k*M_PI*t);
      const tk::real u = ft * z * sin(M_PI*x);
      const tk::real v = ft * z * cos(M_PI*y);
      const tk::real w = ft*(-0.5*M_PI*z*z*(cos(M_PI*x)-sin(M_PI*y)));
      // total specific energy
      const tk::real rE = p/(g-1.0) + 0.5*r*(u*u + v*v + w*w);
      return {{ r, r*u, r*v, r*w, rE }};
    }

    //! \brief Evaluate the increment from t to t+dt of the analytical solution
    //!   at (x,y,z) for all components
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \param[in] t Time where to evaluate the solution increment starting from
    //! \param[in] dt Time increment at which evaluate the solution increment to
    //! \return Increment in values of all components evaluated at (x,y,z,t+dt)
    static std::array< tk::real, 5 >
    solinc( tk::ctr::ncomp_type e,
            tk::real x, tk::real y, tk::real z, tk::real t, tk::real dt )
    {
      auto st1 = solution( e, x, y, z, t );
      auto st2 = solution( e, x, y, z, t+dt );
      std::transform( begin(st1), end(st1), begin(st2), begin(st2),
                      []( tk::real s, tk::real& d ){ return d -= s; } );
      return st2;
    }

  public:

    //! Set initial conditions
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] t Physical time
    static void init( const std::array< std::vector< tk::real >, 3 >& coord,
                      const std::vector< std::size_t >&,
                      tk::Fields& unk,
                      tk::ctr::ncomp_type e,
                      tk::ctr::ncomp_type offset,
                      tk::real t )
    {
      Assert( coord[0].size() == unk.nunk(), "Size mismatch" );
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
      // set initial and boundary conditions
      for (ncomp_t i=0; i<coord[0].size(); ++i) {
        const auto s = solution( e, x[i], y[i], z[i], t );
        unk(i,0,offset) = s[0]; // rho
        unk(i,1,offset) = s[1]; // rho * u
        unk(i,2,offset) = s[2]; // rho * v
        unk(i,3,offset) = s[3]; // rho * w
        unk(i,4,offset) = s[4]; // rho * e, e: total = kinetic + internal energy
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
      using tag::param; using tag::compflow; using std::sin; using std::cos;
      // manufactured solution parameters
      const auto a = g_inputdeck.get< param, compflow, tag::alpha >()[e];
      const auto bx = g_inputdeck.get< param, compflow, tag::betax >()[e];
      const auto by = g_inputdeck.get< param, compflow, tag::betay >()[e];
      const auto bz = g_inputdeck.get< param, compflow, tag::betaz >()[e];
      const auto k = g_inputdeck.get< param, compflow, tag::kappa >()[e];
      // ratio of specific heats
      const tk::real g = g_inputdeck.get< param, compflow, tag::gamma >()[e];

      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // density gradient
      const std::array< std::array< tk::real, 4 >, 3 >
        drdx{{  {{ -2.0*bx*x[N[0]],
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

      // pressure gradient
      const std::array< std::array< tk::real, 4 >, 3 >
        dpdx{{ {{ 2.0*a*bx*x[N[0]],
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

      // solution
      const std::array< std::array< tk::real, 5 >, 4 >
        s{{ solution( e, x[N[0]], y[N[0]], z[N[0]], t ),
            solution( e, x[N[1]], y[N[1]], z[N[1]], t ),
            solution( e, x[N[2]], y[N[2]], z[N[2]], t ),
            solution( e, x[N[3]], y[N[3]], z[N[3]], t ) }};

      // velocity
      const std::array< tk::real, 4 > u{{ s[0][1]/s[0][0],
                                          s[1][1]/s[1][0],
                                          s[2][1]/s[2][0],
                                          s[3][1]/s[3][0] }};
      const std::array< tk::real, 4 > v{{ s[0][2]/s[0][0],
                                          s[1][2]/s[1][0],
                                          s[2][2]/s[2][0],
                                          s[3][2]/s[3][0] }};
      const std::array< tk::real, 4 > w{{ s[0][3]/s[0][0],
                                          s[1][3]/s[1][0],
                                          s[2][3]/s[2][0],
                                          s[3][3]/s[3][0] }};

      // pressure
      const std::array< tk::real, 4 > p{{
        (g-1.0)*(s[0][4] - 0.5*s[0][0]*(u[0]*u[0] + v[0]*v[0] + w[0]*w[0])),
        (g-1.0)*(s[1][4] - 0.5*s[1][0]*(u[1]*u[1] + v[1]*v[1] + w[1]*w[1])),
        (g-1.0)*(s[2][4] - 0.5*s[2][0]*(u[2]*u[2] + v[2]*v[2] + w[2]*w[2])),
        (g-1.0)*(s[3][4] - 0.5*s[3][0]*(u[3]*u[3] + v[3]*v[3] + w[3]*w[3])) }};

      // velocity time derivatives
      const std::array< tk::real, 4 >
        dudt{{ -k*M_PI*sin(k*M_PI*t)*z[N[0]]*sin(M_PI*x[N[0]]),
               -k*M_PI*sin(k*M_PI*t)*z[N[1]]*sin(M_PI*x[N[1]]),
               -k*M_PI*sin(k*M_PI*t)*z[N[2]]*sin(M_PI*x[N[2]]),
               -k*M_PI*sin(k*M_PI*t)*z[N[3]]*sin(M_PI*x[N[3]]) }};
      const std::array< tk::real, 4 >
        dvdt{{ -k*M_PI*sin(k*M_PI*t)*z[N[0]]*cos(M_PI*y[N[0]]),
               -k*M_PI*sin(k*M_PI*t)*z[N[1]]*cos(M_PI*y[N[1]]),
               -k*M_PI*sin(k*M_PI*t)*z[N[2]]*cos(M_PI*y[N[2]]),
               -k*M_PI*sin(k*M_PI*t)*z[N[3]]*cos(M_PI*y[N[3]]) }};
      const std::array< tk::real, 4 >
        dwdt{{ k*M_PI*sin(k*M_PI*t)*0.5*M_PI*z[N[0]]*z[N[0]]
                *(cos(M_PI*x[N[0]]) - sin(M_PI*y[N[0]])),
               k*M_PI*sin(k*M_PI*t)*0.5*M_PI*z[N[1]]*z[N[1]]
                *(cos(M_PI*x[N[1]]) - sin(M_PI*y[N[1]])),
               k*M_PI*sin(k*M_PI*t)*0.5*M_PI*z[N[2]]*z[N[2]]
                *(cos(M_PI*x[N[2]]) - sin(M_PI*y[N[2]])),
               k*M_PI*sin(k*M_PI*t)*0.5*M_PI*z[N[3]]*z[N[3]]
                *(cos(M_PI*x[N[3]]) - sin(M_PI*y[N[3]])) }};

      // velocity spatial derivatives
      const tk::real ft = cos(k*M_PI*t);
      const std::array< std::array< tk::real, 4 >, 3 >
        dudx{{ {{ M_PI*ft*z[N[0]]*cos(M_PI*x[N[0]]),
                  M_PI*ft*z[N[1]]*cos(M_PI*x[N[1]]),
                  M_PI*ft*z[N[2]]*cos(M_PI*x[N[2]]),
                  M_PI*ft*z[N[3]]*cos(M_PI*x[N[3]]) }},
               {{ 0.0, 0.0, 0.0, 0.0 }},
               {{ ft*sin(M_PI*x[N[0]]),
                  ft*sin(M_PI*x[N[1]]),
                  ft*sin(M_PI*x[N[2]]),
                  ft*sin(M_PI*x[N[3]]) }} }};
      const std::array< std::array< tk::real, 4 >, 3 >
        dvdx{{ {{ 0.0, 0.0, 0.0, 0.0 }},
               {{ -M_PI*ft*z[N[0]]*sin(M_PI*y[N[0]]),
                  -M_PI*ft*z[N[1]]*sin(M_PI*y[N[1]]),
                  -M_PI*ft*z[N[2]]*sin(M_PI*y[N[2]]),
                  -M_PI*ft*z[N[3]]*sin(M_PI*y[N[3]]) }},
               {{ ft*cos(M_PI*y[N[0]]),
                  ft*cos(M_PI*y[N[1]]),
                  ft*cos(M_PI*y[N[2]]),
                  ft*cos(M_PI*y[N[3]]) }} }};
      const std::array< std::array< tk::real, 4 >, 3 >
        dwdx{{ {{ M_PI*ft*0.5*M_PI*z[N[0]]*z[N[0]]*sin(M_PI*x[N[0]]),
                  M_PI*ft*0.5*M_PI*z[N[1]]*z[N[1]]*sin(M_PI*x[N[1]]),
                  M_PI*ft*0.5*M_PI*z[N[2]]*z[N[2]]*sin(M_PI*x[N[2]]),
                  M_PI*ft*0.5*M_PI*z[N[3]]*z[N[3]]*sin(M_PI*x[N[3]]) }},
               {{ M_PI*ft*0.5*M_PI*z[N[0]]*z[N[0]]*cos(M_PI*y[N[0]]),
                  M_PI*ft*0.5*M_PI*z[N[1]]*z[N[1]]*cos(M_PI*y[N[1]]),
                  M_PI*ft*0.5*M_PI*z[N[2]]*z[N[2]]*cos(M_PI*y[N[2]]),
                  M_PI*ft*0.5*M_PI*z[N[3]]*z[N[3]]*cos(M_PI*y[N[3]]) }},
               {{ -ft*M_PI*z[N[0]]*(cos(M_PI*x[N[0]]) - sin(M_PI*y[N[0]])),
                  -ft*M_PI*z[N[1]]*(cos(M_PI*x[N[1]]) - sin(M_PI*y[N[1]])),
                  -ft*M_PI*z[N[2]]*(cos(M_PI*x[N[2]]) - sin(M_PI*y[N[2]])),
                  -ft*M_PI*z[N[3]]*(cos(M_PI*x[N[3]]) - sin(M_PI*y[N[3]])) }}}};

      // energy time derivative
      const std::array< tk::real, 4 >
        dedt{{ u[0]*dudt[0] + v[0]*dvdt[0] + w[0]*dwdt[0],
               u[1]*dudt[1] + v[1]*dvdt[1] + w[1]*dwdt[1],
               u[2]*dudt[2] + v[2]*dvdt[2] + w[2]*dwdt[2],
               u[3]*dudt[3] + v[3]*dvdt[3] + w[3]*dwdt[3] }};

      // energy spatial derivatives
      const std::array< std::array< tk::real, 4 >, 3 >
        dedx{{ {{ dpdx[0][0]/s[0][0]/(g-1.0) -
                    p[0]/(g-1.0)/s[0][0]/s[0][0]*drdx[0][0] +
                    u[0]*dudx[0][0] + v[0]*dvdx[0][0] + w[0]*dwdx[0][0],
                  dpdx[0][1]/s[1][0]/(g-1.0) -
                    p[1]/(g-1.0)/s[1][0]/s[1][0]*drdx[0][1] +
                    u[1]*dudx[0][1] + v[1]*dvdx[0][1] + w[1]*dwdx[0][1],
                  dpdx[0][2]/s[2][0]/(g-1.0) -
                    p[2]/(g-1.0)/s[2][0]/s[2][0]*drdx[0][2] +
                    u[2]*dudx[0][2] + v[2]*dvdx[0][2] + w[2]*dwdx[0][2],
                  dpdx[0][3]/s[3][0]/(g-1.0) -
                    p[3]/(g-1.0)/s[3][0]/s[3][0]*drdx[0][3] +
                    u[3]*dudx[0][3] + v[3]*dvdx[0][3] + w[3]*dwdx[0][3] }},
               {{ dpdx[1][0]/s[0][0]/(g-1.0) -
                    p[0]/(g-1.0)/s[0][0]/s[0][0]*drdx[1][0] +
                    u[0]*dudx[1][0] + v[0]*dvdx[1][0] + w[0]*dwdx[1][0],
                  dpdx[1][1]/s[1][0]/(g-1.0) -
                    p[1]/(g-1.0)/s[1][0]/s[1][0]*drdx[1][1] +
                    u[1]*dudx[1][1] + v[1]*dvdx[1][1] + w[1]*dwdx[1][1],
                  dpdx[1][2]/s[2][0]/(g-1.0) -
                    p[2]/(g-1.0)/s[2][0]/s[2][0]*drdx[1][2] +
                    u[2]*dudx[1][2] + v[2]*dvdx[1][2] + w[2]*dwdx[1][2],
                  dpdx[1][3]/s[3][0]/(g-1.0) -
                    p[3]/(g-1.0)/s[3][0]/s[3][0]*drdx[1][3] +
                    u[3]*dudx[1][3] + v[3]*dvdx[1][3] + w[3]*dwdx[1][3] }},
               {{ dpdx[2][0]/s[0][0]/(g-1.0) -
                    p[0]/(g-1.0)/s[0][0]/s[0][0]*drdx[2][0] +
                    u[0]*dudx[2][0] + v[0]*dvdx[2][0] + w[0]*dwdx[2][0],
                  dpdx[2][1]/s[1][0]/(g-1.0) -
                    p[1]/(g-1.0)/s[1][0]/s[1][0]*drdx[2][1] +
                    u[1]*dudx[2][1] + v[1]*dvdx[2][1] + w[1]*dwdx[2][1],
                  dpdx[2][2]/s[2][0]/(g-1.0) -
                    p[2]/(g-1.0)/s[2][0]/s[2][0]*drdx[2][2] +
                    u[2]*dudx[2][2] + v[2]*dvdx[2][2] + w[2]*dwdx[2][2],
                  dpdx[2][3]/s[3][0]/(g-1.0) -
                    p[3]/(g-1.0)/s[3][0]/s[3][0]*drdx[2][3] +
                    u[3]*dudx[2][3] + v[3]*dvdx[2][3] + w[3]*dwdx[2][3] }} }};

      // density source
      const std::array< tk::real, 4 >
        Sr{{ u[0]*drdx[0][0] + v[0]*drdx[1][0] + w[0]*drdx[2][0],
             u[1]*drdx[0][1] + v[1]*drdx[1][1] + w[1]*drdx[2][1],
             u[2]*drdx[0][2] + v[2]*drdx[1][2] + w[2]*drdx[2][2],
             u[3]*drdx[0][3] + v[3]*drdx[1][3] + w[3]*drdx[2][3] }};

      // momentum source
      const std::array< std::array< tk::real, 4 >, 3 >
        Sm{{ {{ s[0][0]*dudt[0] + u[0]*Sr[0] + dpdx[0][0] +
                  s[0][1]*dudx[0][0] + s[0][2]*dudx[1][0] + s[0][3]*dudx[2][0],
                s[1][0]*dudt[1] + u[1]*Sr[1] + dpdx[0][1] +
                  s[1][1]*dudx[0][1] + s[1][2]*dudx[1][1] + s[1][3]*dudx[2][1],
                s[2][0]*dudt[2] + u[2]*Sr[2] + dpdx[0][2] +
                  s[2][1]*dudx[0][2] + s[2][2]*dudx[1][2] + s[2][3]*dudx[2][2],
                s[3][0]*dudt[3] + u[3]*Sr[3] + dpdx[0][3] +
                  s[3][1]*dudx[0][3] + s[3][2]*dudx[1][3] + s[3][3]*dudx[3][3],
             }},
             {{ s[0][0]*dvdt[0] + v[0]*Sr[0] + dpdx[1][0] +
                  s[0][1]*dvdx[0][0] + s[0][2]*dvdx[1][0] + s[0][3]*dvdx[2][0],
                s[1][0]*dvdt[1] + v[1]*Sr[1] + dpdx[1][1] +
                  s[1][1]*dvdx[0][1] + s[1][2]*dvdx[1][1] + s[1][3]*dvdx[2][1],
                s[2][0]*dvdt[2] + v[2]*Sr[2] + dpdx[1][2] +
                  s[2][1]*dvdx[0][2] + s[2][2]*dvdx[1][2] + s[2][3]*dvdx[2][2],
                s[3][0]*dvdt[3] + v[3]*Sr[3] + dpdx[1][3] +
                  s[3][1]*dvdx[0][3] + s[3][2]*dvdx[1][3] + s[3][3]*dvdx[3][3]
             }},
             {{ s[0][0]*dwdt[0] + w[0]*Sr[0] + dpdx[2][0] +
                  s[0][1]*dwdx[0][0] + s[0][2]*dwdx[1][0] + s[0][3]*dwdx[2][0],
                s[1][0]*dwdt[1] + w[1]*Sr[1] + dpdx[2][1] +
                  s[1][1]*dwdx[0][1] + s[1][2]*dwdx[1][1] + s[1][3]*dwdx[2][1],
                s[2][0]*dwdt[2] + w[2]*Sr[2] + dpdx[2][2] +
                  s[2][1]*dwdx[0][2] + s[2][2]*dwdx[1][2] + s[2][3]*dwdx[2][2],
                s[3][0]*dwdt[3] + w[3]*Sr[3] + dpdx[2][3] +
                  s[3][1]*dwdx[0][3] + s[3][2]*dwdx[1][3] + s[3][3]*dwdx[3][3]
             }} }};

      // energy source
      const std::array< tk::real, 4 >
        Se{{ s[0][0]*dedt[0] + s[0][4]/s[0][0]*Sr[0]
               + s[0][1]*dedx[0][0] + s[0][2]*dedx[1][0] + s[0][3]*dedx[2][0]
               + u[0]*dpdx[0][0] + v[0]*dpdx[1][0] + w[0]*dpdx[2][0],
             s[1][0]*dedt[1] + s[1][4]/s[1][0]*Sr[1]
               + s[1][1]*dedx[0][1] + s[1][2]*dedx[1][1] + s[1][3]*dedx[2][1]
               + u[1]*dpdx[0][1] + v[1]*dpdx[1][1] + w[1]*dpdx[2][1],
             s[2][0]*dedt[2] + s[2][4]/s[2][0]*Sr[2]
               + s[2][1]*dedx[0][2] + s[2][2]*dedx[1][2] + s[2][3]*dedx[2][2]
               + u[2]*dpdx[0][2] + v[2]*dpdx[1][2] + w[2]*dpdx[2][2],
             s[3][0]*dedt[3] + s[3][4]/s[3][0]*Sr[3]
               + s[3][1]*dedx[0][3] + s[3][2]*dedx[1][3] + s[3][3]*dedx[3][3]
               + u[3]*dpdx[0][3] + v[3]*dpdx[1][3] + w[3]*dpdx[2][3] }};

      // add density, momentum, and energy source at element nodes
      for (std::size_t alpha=0; alpha<4; ++alpha)
        for (std::size_t beta=0; beta<4; ++beta) {
          // source contribution to mass rhs
          R.var(r[0],N[alpha]) += dt * mass[alpha][beta] * Sr[beta];
          // source contribution to momentum rhs
          for (std::size_t i=0; i<3; ++i)
            R.var(r[i+1],N[alpha]) += dt * mass[alpha][beta] * Sm[i][beta];
          // source contribution to energy rhs
          R.var(r[4],N[alpha]) += dt * mass[alpha][beta] * Se[beta];
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

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] t Physical time
    //! \param[in] deltat Time step size
    //! \param[in] side Pair of side set ID and node IDs on the side set
    //! \param[in] coord Mesh node coordinates
    //! \return Vector of pairs of bool and boundary condition value associated
    //!   to mesh node IDs at which Dirichlet boundary conditions are set. Note
    //!   that instead of the actual boundary condition value, we return the
    //!   increment between t+dt and t, since that is what the solution requires
    //!   as we solve for the soution increments and not the solution itself.
    static std::unordered_map< std::size_t,
                               std::vector< std::pair< bool, tk::real > > >
    dirbc( tk::ctr::ncomp_type e,
           tk::real t,
           tk::real deltat,
           const std::pair< const int, std::vector< std::size_t > >& side,
           const std::array< std::vector< tk::real >, 3 >& coord )
    {
      using tag::param; using tag::compflow; using tag::bcdir;
      using NodeBC = std::vector< std::pair< bool, tk::real > >;
      std::unordered_map< std::size_t, NodeBC > bc;
      const auto& ubc = g_inputdeck.get< param, compflow, bcdir >();
      if (!ubc.empty()) {
        Assert( ubc.size() > e, "Indexing out of Dirichlet BC eq-vector" );
        const auto& x = coord[0];
        const auto& y = coord[1];
        const auto& z = coord[2];
        for (const auto& b : ubc[e])
          if (std::stoi(b) == side.first)
            for (auto n : side.second) {
              Assert( x.size() > n, "Indexing out of coordinate array" );
              auto s = solinc( e, x[n], y[n], z[n], t, deltat );
              bc[n] = {{ {true,s[0]}, {true,s[1]}, {true,s[2]}, {true,s[3]},
                         {true,s[4]} }};
            }
      }
      return bc;
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    static std::vector< std::string > fieldNames() {
      std::vector< std::string > n;
      n.push_back( "density_numerical" );
      n.push_back( "x-velocity_numerical" );
      n.push_back( "y-velocity_numerical" );
      n.push_back( "z-velocity_numerical" );
      n.push_back( "specific_total_energy_numerical" );
      n.push_back( "pressure_numerical" );
      n.push_back( "density_analytical" );
      n.push_back( "x-velocity_analytical" );
      n.push_back( "y-velocity_analytical" );
      n.push_back( "z-velocity_analytical" );
      n.push_back( "specific_total_energy_analytical" );
      n.push_back( "pressure_analytical" );
      n.push_back( "err(rho)" );
      n.push_back( "err(e)" );
      n.push_back( "err(p)" );
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
                 tk::real V,
                 const std::vector< tk::real >& vol,
                 const std::array< std::vector< tk::real >, 3 >& coord,
                 tk::Fields& U )
    {
      // ratio of specific heats
      tk::real g =
        g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[e];

      std::vector< std::vector< tk::real > > out;
      auto r = U.extract( 0, offset );
      auto u = U.extract( 1, offset );
      auto v = U.extract( 2, offset );
      auto w = U.extract( 3, offset );
      auto E = U.extract( 4, offset );

      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      out.push_back( r );
      std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( u );
      std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( v );
      std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( w );
      std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( E );

      auto p = r;
      for (std::size_t i=0; i<r.size(); ++i)
        p[i] = (g-1.0)*r[i]*(E[i] - (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0);
      out.push_back( p );

      auto er = r, ee = r, ep = r;
      for (std::size_t i=0; i<r.size(); ++i) {
        auto s = solution( e, x[i], y[i], z[i], t );
        er[i] = std::pow( r[i] - s[0], 2.0 ) * vol[i] / V;
        ee[i] = std::pow( E[i] - s[4]/s[0], 2.0 ) * vol[i] / V;
        auto ap = (g-1.0)*(s[4] - (s[1]*s[1] + s[2]*s[2] + s[3]*s[3])/2.0/s[0]);
        r[i] = s[0];
        u[i] = s[1]/s[0];
        v[i] = s[2]/s[0];
        w[i] = s[3]/s[0];
        E[i] = s[4]/s[0];
        p[i] = (g-1.0)*r[i]*(E[i] - (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0);
        ep[i] = std::pow( ap - p[i], 2.0 ) * vol[i] / V;
      }

      out.push_back( r );
      out.push_back( u );
      out.push_back( v );
      out.push_back( w );
      out.push_back( E );
      out.push_back( p );

      out.push_back( er );
      out.push_back( ee );
      out.push_back( ep );

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
