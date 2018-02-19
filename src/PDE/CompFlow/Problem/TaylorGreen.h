// *****************************************************************************
/*!
  \file      src/PDE/CompFlowProblem/TaylorGreen.h
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problems.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowProblemTaylorGreen_h
#define CompFlowProblemTaylorGreen_h

#include <string>
#include <unordered_set>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! CompFlow system of PDEs problem: Taylor-Green
//! \see G.I. Taylor, A.E. Green, "Mechanism of the Production of Small Eddies
//!   from Large Ones", Proc. R. Soc. Lond. A 1937 158 499-521; DOI:
//!   10.1098/rspa.1937.0036. Published 3 February 1937
//! \see Waltz, et. al, "Verification of a three-dimensional unstructured finite
//!   element method using analytic and manufactured solutions", Computers and
//!   Fluids, 2013, Vol.81, pp.57-67.
class CompFlowProblemTaylorGreen {

  public:

    //! Evaluate analytical solution at (x,y,0) for all components
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \return Values of all components evaluated at (x,y,0)
    static std::array< tk::real, 5 >
    solution( tk::ctr::ncomp_type e,
              tk::real x, tk::real y, tk::real, tk::real )
    {
      using tag::param; using tag::compflow; using std::sin; using std::cos;
      // ratio of specific heats
      const tk::real g = g_inputdeck.get< param, compflow, tag::gamma >()[e];
      // density
      const tk::real r = 1.0;
      // pressure
      const tk::real p = 10.0 + r/4.0*(cos(2.0*M_PI*x) + cos(2.0*M_PI*y));
      // velocity
      const tk::real u =  sin(M_PI*x) * cos(M_PI*y);
      const tk::real v = -cos(M_PI*x) * sin(M_PI*y);
      const tk::real w = 0.0;
      // total specific energy
      const tk::real rE = p/(g-1.0) + 0.5*r*(u*u + v*v + w*w);
      return {{ r, r*u, r*v, r*w, rE }};
    }

    //! \brief Evaluate the increment from t to t+dt of the analytical solution
    //!   at (x,y,z) for all components
    //! \return Increment in values of all components: all zero for this problem
    static std::array< tk::real, 5 >
    solinc( tk::ctr::ncomp_type,
            tk::real, tk::real, tk::real, tk::real, tk::real )
    {
      return {{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
    }

    //! Compute and return source term for Rayleigh-Taylor manufactured solution
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \return Array of reals containing the source for all components
    static std::array< tk::real, 5 >
    src( tk::ctr::ncomp_type, tk::real x, tk::real y, tk::real, tk::real ) {
      using tag::param; using tag::compflow; using std::sin; using std::cos;
      // ratio of specific heats
      std::array< tk::real, 5 > r{{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
      // only energy source
      r[4] = 3.0*M_PI/8.0*( cos(3.0*M_PI*x)*cos(M_PI*y) -
                            cos(3.0*M_PI*y)*cos(M_PI*x) );
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
      n.push_back( "density_analytical" );
      n.push_back( "x-velocity_numerical" );
      n.push_back( "x-velocity_analytical" );
      n.push_back( "err(u)" );
      n.push_back( "y-velocity_numerical" );
      n.push_back( "y-velocity_analytical" );
      n.push_back( "err(v)" );
      n.push_back( "z-velocity_numerical" );
      n.push_back( "z-velocity_analytical" );
      n.push_back( "specific_total_energy_numerical" );
      n.push_back( "specific_total_energy_analytical" );
      n.push_back( "err(E)" );
      n.push_back( "pressure_numerical" );
      n.push_back( "pressure_analytical" );
      return n;
    }

    //! Return field output going to file
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] V Total mesh volume
    //! \param[in] vol Nodal mesh volumes
    //! \param[in] coord Mesh node coordinates
    //! \param[in] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    static std::vector< std::vector< tk::real > >
    fieldOutput( tk::ctr::ncomp_type e,
                 tk::ctr::ncomp_type offset,
                 tk::real,
                 tk::real V,
                 const std::vector< tk::real >& vol,
                 const std::array< std::vector< tk::real >, 3 >& coord,
                 tk::Fields& U )
    {
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

      out.push_back( r );
      out.push_back( std::vector< tk::real >( r.size(), 1.0 ) );

      std::vector< tk::real > u = ru;
      std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( u );
      std::vector< tk::real > ua = ru;
      for (std::size_t i=0; i<ua.size(); ++i)
        ua[i] = std::sin(M_PI*x[i]) * std::cos(M_PI*y[i]);
      out.push_back( ua );

      // error in x-velocity
      auto err = u;
      for (std::size_t i=0; i<u.size(); ++i)
         err[i] = std::pow( ua[i] - u[i], 2.0 ) * vol[i] / V;
       out.push_back( err );

      std::vector< tk::real > v = rv;
      std::vector< tk::real > va = rv;
      std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( v );
      for (std::size_t i=0; i<va.size(); ++i)
        va[i] = -std::cos(M_PI*x[i]) * std::sin(M_PI*y[i]);
      out.push_back( va );

      // error in v-velocity
      for (std::size_t i=0; i<v.size(); ++i)
        err[i] = std::pow( va[i] - v[i], 2.0 ) * vol[i] / V;
      out.push_back( err );

      std::vector< tk::real > w = rw;
      std::vector< tk::real > wa = rw;
      std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( w );
      for (std::size_t i=0; i<wa.size(); ++i)
        wa[i] = 0.0;
      out.push_back( wa );

      std::vector< tk::real > E = re;
      std::vector< tk::real > Ea = re;
      std::vector< tk::real > Pa( r.size(), 0.0 );
      std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( E );
      for (std::size_t i=0; i<Ea.size(); ++i) {
        Pa[i] = 10.0 +
          r[i]/4.0*(std::cos(2.0*M_PI*x[i]) + std::cos(2.0*M_PI*y[i]));
        Ea[i] = Pa[i]/(g-1.0)/r[i] +
                0.5*(ua[i]*ua[i] + va[i]*va[i] + wa[i]*wa[i])/r[i];
      }
      out.push_back( Ea );

      // error in total specific energy
      for (std::size_t i=0; i<v.size(); ++i)
        err[i] = std::pow( Ea[i] - E[i], 2.0 ) * vol[i] / V;
      out.push_back( err );

      std::vector< tk::real > P( r.size(), 0.0 );
      for (std::size_t i=0; i<P.size(); ++i)
        P[i] = (g-1.0)*(E[i] - (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0/r[i]);
      out.push_back( P );
      out.push_back( Pa );

      return out;
   }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    static std::vector< std::string > names()
    { return { "r", "ru", "rv", "rw", "re" }; }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::TAYLOR_GREEN; }
};

} // inciter::

#endif // CompFlowProblemTaylorGreen_h
