// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/NLEnergyGrowth.h
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problems.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowProblemNLEnergyGrowth_h
#define CompFlowProblemNLEnergyGrowth_h

#include <string>
#include <unordered_set>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! CompFlow system of PDEs problem: nonlinear energy growth (NLEG)
//! \see Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
//!   equations with relevance to Inertial Confinement Fusion", Journal of
//!   Computational Physics 267 (2014) 196-209.
class CompFlowProblemNLEnergyGrowth {

  private:
    //! Compute internal energy parameter
    //! \param[in] bx Parameter betax
    //! \param[in] by Parameter betay
    //! \param[in] bz Parameter betaz
    //! \param[in] x X coordinate to evaluate at
    //! \param[in] y Y coordinate to evaluate at
    //! \param[in] z Z coordinate to evaluate at
    //! \return Internal energy parameter
    static tk::real hx( tk::real bx, tk::real by, tk::real bz,
                        tk::real x, tk::real y, tk::real z )
    { return std::cos(bx*M_PI*x) * std::cos(by*M_PI*y) * std::cos(bz*M_PI*z); }

    //! Compute a power of the internal energy
    //! \param[in] ce Internal energy parameter
    //! \param[in] kappa Internal energy parameter
    //! \param[in] t Physical time
    //! \param[in] h Internal energy parameter
    //! \param[in] p Power
    //! \return Internal energy raised to power p
    static tk::real ec( tk::real ce, tk::real kappa, tk::real t, tk::real h,
                        tk::real p )
    { return std::pow( -3.0*(ce + kappa*h*h*t), p ); }

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
      using tag::param; using tag::compflow;
      // manufactured solution parameters
      const auto ce = g_inputdeck.get< param, compflow, tag::ce >()[e];
      const auto r0 = g_inputdeck.get< param, compflow, tag::r0 >()[e];
      const auto a = g_inputdeck.get< param, compflow, tag::alpha >()[e];
      const auto k = g_inputdeck.get< param, compflow, tag::kappa >()[e];
      const auto bx = g_inputdeck.get< param, compflow, tag::betax >()[e];
      const auto by = g_inputdeck.get< param, compflow, tag::betay >()[e];
      const auto bz = g_inputdeck.get< param, compflow, tag::betaz >()[e];
      // spatial component of density field
      const tk::real gx = 1.0 - x*x - y*y - z*z;
      // internal energy parameter
      const auto h = hx( bx, by, bz, x, y, z );
      // temporal component of the density field
      tk::real ft = std::exp( -a*t );
      // solution at t
      auto r = r0 + ft*gx;
      return {{ r, 0.0, 0.0, 0.0, r*ec(ce,k,t,h,-1.0/3.0) }};
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

    //! Compute and return source term for NLEG manufactured solution
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
      const auto a = g_inputdeck.get< param, compflow, tag::alpha >()[e];
      const auto bx = g_inputdeck.get< param, compflow, tag::betax >()[e];
      const auto by = g_inputdeck.get< param, compflow, tag::betay >()[e];
      const auto bz = g_inputdeck.get< param, compflow, tag::betaz >()[e];
      const auto ce = g_inputdeck.get< param, compflow, tag::ce >()[e];
      const auto kappa = g_inputdeck.get< param, compflow, tag::kappa >()[e];
      const auto r0 = g_inputdeck.get< param, compflow, tag::r0 >()[e];
      // ratio of specific heats
      const auto g = g_inputdeck.get< param, compflow, tag::gamma >()[e];
      // spatial component of density field
      const auto gx = 1.0 - x*x - y*y - z*z;
      // derivative of spatial component of density field
      const std::array< tk::real, 3 > dg{{ -2.0*x, -2.0*y, -2.0*z }};
      // spatial component of energy field
      const auto h = hx( bx, by, bz, x, y, z );
      // derivative of spatial component of energy field
      std::array< tk::real, 3 >
        dh{{ -bx*M_PI*sin(bx*M_PI*x)*cos(by*M_PI*y)*cos(bz*M_PI*z),
             -by*M_PI*cos(bx*M_PI*x)*sin(by*M_PI*y)*cos(bz*M_PI*z),
             -bz*M_PI*cos(bx*M_PI*x)*cos(by*M_PI*y)*sin(bz*M_PI*z) }};
      // temporal function f and its derivative
      const auto ft = std::exp(-a*t);
      const auto dfdt = -a*ft;
      // density and its derivatives
      const auto rho = r0 + ft*gx;
      const std::array< tk::real, 3 > drdx{{ ft*dg[0], ft*dg[1], ft*dg[2] }};
      const auto drdt = gx*dfdt;
      // internal energy and its derivatives
      const auto ie = ec( ce, kappa, t, h, -1.0/3.0 );
      const std::array< tk::real, 3 > dedx{{
        2.0*std::pow(ie,4.0)*kappa*h*dh[0]*t,
        2.0*std::pow(ie,4.0)*kappa*h*dh[1]*t,
        2.0*std::pow(ie,4.0)*kappa*h*dh[2]*t }};
      const auto dedt = kappa*h*h*std::pow(ie,4.0);
      // sources
      std::array< tk::real, 5 > r;
      // density source
      r[0] = drdt;
      // momentum source
      r[1] = (g-1.0)*(rho*dedx[0] + ie*drdx[0]);
      r[2] = (g-1.0)*(rho*dedx[1] + ie*drdx[1]);
      r[3] = (g-1.0)*(rho*dedx[2] + ie*drdx[2]);
      // energy source
      r[4] = rho*dedt + ie*drdt;
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
      auto r  = U.extract( 0, offset );
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

      auto er = r, ee = r;
      for (std::size_t i=0; i<r.size(); ++i) {
        auto s = solution( e, x[i], y[i], z[i], t );
        er[i] = std::pow( r[i] - s[0], 2.0 ) * vol[i] / V;
        ee[i] = std::pow( E[i] - s[4]/s[0], 2.0 ) * vol[i] / V;
        r[i] = s[0];
        u[i] = s[1]/s[0];
        v[i] = s[2]/s[0];
        w[i] = s[3]/s[0];
        E[i] = s[4]/s[0];
        p[i] = (g-1.0)*r[i]*(E[i] - (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0);
      }

      out.push_back( r );
      out.push_back( u );
      out.push_back( v );
      out.push_back( w );
      out.push_back( E );
      out.push_back( p );

      out.push_back( er );
      out.push_back( ee );

      return out;
   }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    static std::vector< std::string > names()
    { return { "r", "ru", "rv", "rw", "re" }; }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::NL_ENERGY_GROWTH; }
};

} // inciter::

#endif // CompFlowProblemNLEnergyGrowth_h
