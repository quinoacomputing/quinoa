// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/RayleighTaylor.h
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problems.h
    for general requirements on Problem policy classes for CompFlow.
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
      const auto a = g_inputdeck.get< param, compflow, tag::alpha >()[e];
      const auto bx = g_inputdeck.get< param, compflow, tag::betax >()[e];
      const auto by = g_inputdeck.get< param, compflow, tag::betay >()[e];
      const auto bz = g_inputdeck.get< param, compflow, tag::betaz >()[e];
      const auto p0 = g_inputdeck.get< param, compflow, tag::p0 >()[e];
      const auto r0 = g_inputdeck.get< param, compflow, tag::r0 >()[e];
      const auto k = g_inputdeck.get< param, compflow, tag::kappa >()[e];
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
      const tk::real u = ft*z*sin(M_PI*x);
      const tk::real v = ft*z*cos(M_PI*y);
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

    //! Compute and return source term for Rayleigh-Taylor manufactured solution
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \param[in] t Physical time at which to evaluate the source
    //! \return Array of reals containing the source for all components
    static std::array< tk::real, 5 >
    src( tk::ctr::ncomp_type e, tk::real x, tk::real y, tk::real z, tk::real t )
    {
      using tag::param; using tag::compflow; using std::sin; using std::cos;

      // manufactured solution parameters
      auto a = g_inputdeck.get< param, compflow, tag::alpha >()[e];
      auto bx = g_inputdeck.get< param, compflow, tag::betax >()[e];
      auto by = g_inputdeck.get< param, compflow, tag::betay >()[e];
      auto bz = g_inputdeck.get< param, compflow, tag::betaz >()[e];
      auto k = g_inputdeck.get< param, compflow, tag::kappa >()[e];
      auto p0 = g_inputdeck.get< param, compflow, tag::p0 >()[e];
      // ratio of specific heats
      tk::real g = g_inputdeck.get< param, compflow, tag::gamma >()[e];

      // evaluate solution at x,y,z,t
      auto s = solution( e, x, y, z, t );

      // density, velocity, energy, pressure
      auto rho = s[0];
      auto u = s[1]/s[0];
      auto v = s[2]/s[0];
      auto w = s[3]/s[0];
      auto E = s[4]/s[0];
      auto p = p0 + a*(bx*x*x + by*y*y + bz*z*z);

      // spatial gradients
      std::array< tk::real, 3 > drdx{{ -2.0*bx*x, -2.0*by*y, -2.0*bz*z }};
      std::array< tk::real, 3 > dpdx{{ 2.0*a*bx*x, 2.0*a*by*y, 2.0*a*bz*z }};
      tk::real ft = cos(k*M_PI*t);
      std::array< tk::real, 3 > dudx{{ ft*M_PI*z*cos(M_PI*x),
                                       0.0,
                                       ft*sin(M_PI*x) }};
      std::array< tk::real, 3 > dvdx{{ 0.0,
                                       -ft*M_PI*z*sin(M_PI*y),
                                       ft*cos(M_PI*y) }};
      std::array< tk::real, 3 > dwdx{{ ft*M_PI*0.5*M_PI*z*z*sin(M_PI*x),
                                       ft*M_PI*0.5*M_PI*z*z*cos(M_PI*y),
                                      -ft*M_PI*z*(cos(M_PI*x) - sin(M_PI*y)) }};
      std::array< tk::real, 3 > dedx{{
        dpdx[0]/rho/(g-1.0) - p/(g-1.0)/rho/rho*drdx[0]
        + u*dudx[0] + v*dvdx[0] + w*dwdx[0],
        dpdx[1]/rho/(g-1.0) - p/(g-1.0)/rho/rho*drdx[1]
        + u*dudx[1] + v*dvdx[1] + w*dwdx[1],
        dpdx[2]/rho/(g-1.0) - p/(g-1.0)/rho/rho*drdx[2]
        + u*dudx[2] + v*dvdx[2] + w*dwdx[2] }};
      
      // time derivatives
      auto dudt = -k*M_PI*sin(k*M_PI*t)*z*sin(M_PI*x);
      auto dvdt = -k*M_PI*sin(k*M_PI*t)*z*cos(M_PI*y);
      auto dwdt =  k*M_PI*sin(k*M_PI*t)/2*M_PI*z*z*(cos(M_PI*x) - sin(M_PI*y));
      auto dedt = u*dudt + v*dvdt + w*dwdt;

      std::array< tk::real, 5 > r;
      // density source
      r[0] = u*drdx[0] + v*drdx[1] + w*drdx[2];
      // momentum source
      r[1] = rho*dudt+u*r[0]+dpdx[0] + s[1]*dudx[0]+s[2]*dudx[1]+s[3]*dudx[2];
      r[2] = rho*dvdt+v*r[0]+dpdx[1] + s[1]*dvdx[0]+s[2]*dvdx[1]+s[3]*dvdx[2];
      r[3] = rho*dwdt+w*r[0]+dpdx[2] + s[1]*dwdx[0]+s[2]*dwdx[1]+s[3]*dwdx[2];
      // energy source
      r[4] = rho*dedt + E*r[0] + s[1]*dedx[0]+s[2]*dedx[1]+s[3]*dedx[2]
             + u*dpdx[0]+v*dpdx[1]+w*dpdx[2];

      return r;
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
      n.push_back( "err(u)" );
      n.push_back( "err(v)" );
      n.push_back( "err(w)" );
      return n;
    }

    //! Return field output going to file
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] t Physical time
    //! \param[in] coord Mesh node coordinates
    //! \param[in] U Solution vector at recent time step
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

      auto er = r, ee = r, ep = r, eu = r, ev = r, ew = r;
      for (std::size_t i=0; i<r.size(); ++i) {
        auto s = solution( e, x[i], y[i], z[i], t );
        er[i] = std::pow( r[i] - s[0], 2.0 ) * vol[i] / V;
        ee[i] = std::pow( E[i] - s[4]/s[0], 2.0 ) * vol[i] / V;
        eu[i] = std::pow( u[i] - s[1]/s[0], 2.0 ) * vol[i] / V;
        ev[i] = std::pow( v[i] - s[2]/s[0], 2.0 ) * vol[i] / V;
        ew[i] = std::pow( w[i] - s[3]/s[0], 2.0 ) * vol[i] / V;
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
      out.push_back( eu );
      out.push_back( ev );
      out.push_back( ew );

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
