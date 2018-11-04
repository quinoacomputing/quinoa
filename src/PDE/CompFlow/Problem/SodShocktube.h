// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/SodShocktube.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Problem configuration for Sod's shock-tube
  \details   This file defines a policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problems.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowProblemSodShocktube_h
#define CompFlowProblemSodShocktube_h

#include <string>
#include <unordered_set>

#include "Types.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! CompFlow system of PDEs problem: Sod shock-tube
//! \see G. A. Sod. A Survey of Several Finite Difference Methods for Systems of
//!   Nonlinear Hyperbolic Conservation Laws. J. Comput. Phys., 27:1–31, 1978.
class CompFlowProblemSodShocktube {

  public:

    //! Evaluate analytical solution at (x,y,0) for all components
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] x X coordinate where to evaluate the solution
//    //! \param[in] t Time at which to evaluate the solution
    //! \return Values of all components evaluated at (x,y,0)
    static std::array< tk::real, 5 >
    solution( tk::ctr::ncomp_type e,
              tk::real x, tk::real, tk::real, tk::real /*t*/ )
    {
      using tag::param; using tag::compflow;
      // ratio of specific heats
      const tk::real g = g_inputdeck.get< param, compflow, tag::gamma >()[e];
      tk::real r, p, u, v, w, rE;
      if (x<0.5) {
        // density
        r = 1.0;
        // pressure
        p = 1.0;
        // velocity
        u = 0.0;
        v = 0.0;
        w = 0.0;
      }
      else {
        // density
        r = 0.125;
        // pressure
        p = 0.1;
        // velocity
        u = 0.0;
        v = 0.0;
        w = 0.0;
      }
      // total specific energy
      rE = p/(g-1.0) + 0.5*r*(u*u + v*v + w*w);
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

    //! Compute and return source term for this problem
    //! \return Array of reals containing the source which is zero for this
    //!   problem
    static std::array< tk::real, 5 >
    src( tk::ctr::ncomp_type, tk::real, tk::real, tk::real, tk::real ) {
      return {{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
    }

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    static void side( std::unordered_set< int >& conf ) {
      using tag::param; using tag::compflow;

      for (const auto& s : g_inputdeck.get< param, compflow,
                                            tag::bcextrapolate >())
        for (const auto& i : s) conf.insert( std::stoi(i) );

      for (const auto& s : g_inputdeck.get< param, compflow, tag::bcsym >())
        for (const auto& i : s) conf.insert( std::stoi(i) );
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    static std::vector< std::string > fieldNames() {
      std::vector< std::string > n;
      n.push_back( "density_numerical" );
      //n.push_back( "density_analytical" );
      n.push_back( "x-velocity_numerical" );
      //n.push_back( "x-velocity_analytical" );
      //n.push_back( "err(u)" );
      n.push_back( "y-velocity_numerical" );
      //n.push_back( "y-velocity_analytical" );
      n.push_back( "z-velocity_numerical" );
      //n.push_back( "z-velocity_analytical" );
      n.push_back( "specific_total_energy_numerical" );
      //n.push_back( "specific_total_energy_analytical" );
      //n.push_back( "err(E)" );
      n.push_back( "pressure_numerical" );
      //n.push_back( "pressure_analytical" );
      return n;
    }

    //! Return field output going to file
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
//    //! \param[in] V Total mesh volume
//    //! \param[in] vol Nodal mesh volumes
//    //! \param[in] coord Mesh node coordinates
    //! \param[in] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    static std::vector< std::vector< tk::real > >
    fieldOutput( tk::ctr::ncomp_type e,
                 tk::ctr::ncomp_type offset,
                 tk::real,
                 tk::real /*V*/,
                 const std::vector< tk::real >& /*vol*/,
                 const std::array< std::vector< tk::real >, 3 >& /*coord*/,
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
      //const auto& x = coord[0];
      //const auto& y = coord[1];

      out.push_back( r );
      //out.push_back( std::vector< tk::real >( r.size(), 1.0 ) );

      std::vector< tk::real > u = ru;
      std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( u );
      //std::vector< tk::real > ua = ru;
      //for (std::size_t i=0; i<ua.size(); ++i)
      //  ua[i] = std::sin(M_PI*x[i]) * std::cos(M_PI*y[i]);
      //out.push_back( ua );

      //// error in x-velocity
      //auto err = u;
      //for (std::size_t i=0; i<u.size(); ++i)
      //   err[i] = std::pow( ua[i] - u[i], 2.0 ) * vol[i] / V;
      // out.push_back( err );

      std::vector< tk::real > v = rv;
      //std::vector< tk::real > va = rv;
      std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( v );
      //for (std::size_t i=0; i<va.size(); ++i)
      //  va[i] = -std::cos(M_PI*x[i]) * std::sin(M_PI*y[i]);
      //out.push_back( va );

      std::vector< tk::real > w = rw;
      //std::vector< tk::real > wa = rw;
      std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( w );
      //for (std::size_t i=0; i<wa.size(); ++i)
      //  wa[i] = 0.0;
      //out.push_back( wa );

      std::vector< tk::real > E = re;
      //std::vector< tk::real > Ea = re;
      //std::vector< tk::real > Pa( r.size(), 0.0 );
      std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( E );
      //for (std::size_t i=0; i<Ea.size(); ++i) {
      //  Pa[i] = 10.0 +
      //    r[i]/4.0*(std::cos(2.0*M_PI*x[i]) + std::cos(2.0*M_PI*y[i]));
      //  Ea[i] = Pa[i]/(g-1.0)/r[i] +
      //          0.5*(ua[i]*ua[i] + va[i]*va[i] + wa[i]*wa[i])/r[i];
      //}
      //out.push_back( Ea );

      //// error in total specific energy
      //for (std::size_t i=0; i<v.size(); ++i)
      //  err[i] = std::pow( Ea[i] - E[i], 2.0 ) * vol[i] / V;
      //out.push_back( err );

      std::vector< tk::real > P( r.size(), 0.0 );
      for (std::size_t i=0; i<P.size(); ++i)
        P[i] = (g-1.0)*r[i]*(E[i] - (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0);
      out.push_back( P );
      //out.push_back( Pa );

      return out;
   }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    static std::vector< std::string > names()
    { return { "r", "ru", "rv", "rw", "re" }; }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SOD_SHOCKTUBE; }
};

} // inciter::

#endif // CompFlowProblemSodShocktube_h
