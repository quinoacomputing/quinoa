// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/VorticalFlow.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Problem configuration for the single-material compressible flow
    equations
  \details   This file defines a Problem policy class for the single-material
    compressible flow equations, defined under PDE/CompFlow/. See
    PDE/CompFlow/Problem.h for general requirements on Problem policy classes
    for CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowProblemVorticalFlow_h
#define CompFlowProblemVorticalFlow_h

#include <string>
#include <unordered_set>

#include "Types.h"
#include "FunctionPrototypes.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! CompFlow system of PDEs problem: vortical flow
//! \see Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
//!   equations with relevance to Inertial Confinement Fusion", Journal of
//!   Computational Physics 267 (2014) 196-209.
class CompFlowProblemVorticalFlow {

  private:
    using ncomp_t = tk::ctr::ncomp_type;
    static constexpr ncomp_t m_ncomp = 5;    //!< Number of scalar components

  public:
    //! Evaluate analytical solution at (x,y,z) for all components
    //! \param[in] system Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \return Values of all components evaluated at (x,y,z)
    //! \note The function signature must follow tk::SolutionFn
    static tk::SolutionFn::result_type
    solution( ncomp_t system, ncomp_t ncomp, tk::real x, tk::real y, tk::real z,
              tk::real )
    {
      Assert( ncomp == m_ncomp, "Number of scalar components must be " +
                                std::to_string(m_ncomp) );
      IGNORE(ncomp);
      using tag::param; using tag::compflow;
      // manufactured solution parameters
      const auto& a =
        g_inputdeck.get< param, compflow, tag::alpha >()[ system ];
      const auto& b = g_inputdeck.get< param, compflow, tag::beta >()[ system ];
      const auto& p0 = g_inputdeck.get< param, compflow, tag::p0 >()[ system ];
      // ratio of specific heats
      tk::real g = g_inputdeck.get< param, compflow, tag::gamma >()[ system ];
      // velocity
      const tk::real ru = a*x - b*y;
      const tk::real rv = b*x + a*y;
      const tk::real rw = -2.0*a*z;
      // total specific energy
      const tk::real rE = (ru*ru+rv*rv+rw*rw)/2.0 + (p0-2.0*a*a*z*z)/(g-1.0);
      return {{ 1.0, ru, rv, rw, rE }};
    }

    //! \brief Evaluate the increment from t to t+dt of the analytical solution
    //!   at (x,y,z) for all components
    //! \return Increment in values of all components: all zero for this problem
    static std::vector< tk::real >
    solinc( ncomp_t, tk::real, tk::real, tk::real, tk::real, tk::real ) {
      return {{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
    }

    //! Compute and return source term for vortical flow manufactured solution
    //! \param[in] system Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \return Array of reals containing the source for all components
    //! \note The function signature must follow tk::SrcFn
    static tk::SrcFn::result_type
    src( ncomp_t system, tk::real x, tk::real y, tk::real z, tk::real ) {
      using tag::param; using tag::compflow;
      // manufactured solution parameters
      const auto& a =
        g_inputdeck.get< param, compflow, tag::alpha >()[ system ];
      const auto& b = g_inputdeck.get< param, compflow, tag::beta >()[ system ];
      // ratio of specific heats
      tk::real g = g_inputdeck.get< param, compflow, tag::gamma >()[ system ];
      // evaluate solution at x,y,z
      auto s = solution( system, m_ncomp, x, y, z, 0.0 );
      std::vector< tk::real > r( m_ncomp );
      // density source
      r[0] = 0.0;
      // momentum source
      r[1] = a*s[1]/s[0] - b*s[2]/s[0];
      r[2] = b*s[1]/s[0] + a*s[2]/s[0];
      r[3] = 0.0;
      // energy source
      r[4] = (r[1]*s[1] + r[2]*s[2])/s[0] + 8.0*a*a*a*z*z/(g-1.0);
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
      n.push_back( "y-velocity_numerical" );
      n.push_back( "y-velocity_analytical" );
      n.push_back( "z-velocity_numerical" );
      n.push_back( "z-velocity_analytical" );
      n.push_back( "specific_total_energy_numerical" );
      n.push_back( "specific_total_energy_analytical" );
      n.push_back( "pressure_numerical" );
      n.push_back( "pressure_analytical" );
      return n;
    }

    //! Return field output going to file
    //! \param[in] system Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] coord Mesh node coordinates
    //! \param[in] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    static std::vector< std::vector< tk::real > >
    fieldOutput( ncomp_t system,
                 ncomp_t offset,
                 tk::real,
                 tk::real,
                 const std::vector< tk::real >&,
                 const std::array< std::vector< tk::real >, 3 >& coord,
                 tk::Fields& U )
    {
      // manufactured solution parameters
      const auto& a =
        g_inputdeck.get< tag::param, tag::compflow, tag::alpha >()[system];
      const auto& b =
        g_inputdeck.get< tag::param, tag::compflow, tag::beta >()[system];
      const auto& p0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::p0 >()[system];
      // number of degree of freedom
      const std::size_t ndof = 
        g_inputdeck.get< tag::discr, tag::ndof >();
      // ratio of specific heats
      tk::real g =
        g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[system];

      std::vector< std::vector< tk::real > > out;
      const auto r  = U.extract( 0*ndof, offset );
      const auto ru = U.extract( 1*ndof, offset );
      const auto rv = U.extract( 2*ndof, offset );
      const auto rw = U.extract( 3*ndof, offset );
      const auto re = U.extract( 4*ndof, offset );

      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      out.push_back( r );
      out.push_back( std::vector< tk::real >( r.size(), 1.0 ) );

      std::vector< tk::real > u = ru;
      std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( u );
      for (std::size_t i=0; i<u.size(); ++i) u[i] = a*x[i] - b*y[i];
      out.push_back( u );

      std::vector< tk::real > v = rv;
      std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( v );
      for (std::size_t i=0; i<v.size(); ++i) v[i] = b*x[i] + a*y[i];
      out.push_back( v );

      std::vector< tk::real > w = rw;
      std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( w );
      for (std::size_t i=0; i<w.size(); ++i) w[i] = -2.0*a*z[i];
      out.push_back( w );

      std::vector< tk::real > E = re;
      std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( E );
      for (std::size_t i=0; i<E.size(); ++i)
        E[i] = 0.5*(u[i]*u[i] + v[i]*v[i] + w[i]*w[i]) +
               (p0 - 2.0*a*a*z[i]*z[i])/(g-1.0);
      out.push_back( E );

      std::vector< tk::real > P( r.size(), 0.0 );
      for (std::size_t i=0; i<P.size(); ++i)
        P[i] = (g-1.0)*(re[i] - r[i]*(u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0);
      out.push_back( P );
      for (std::size_t i=0; i<P.size(); ++i)
        P[i] = p0 - 2.0*a*a*z[i]*z[i];
      out.push_back( P );

      return out;
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    static std::vector< std::string > names()
    { return { "r", "ru", "rv", "rw", "re" }; }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::VORTICAL_FLOW; }
};

} // inciter::

#endif // CompFlowProblemVorticalFlow_h
