// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/UserDefined.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problem.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowProblemUserDefined_h
#define CompFlowProblemUserDefined_h

#include <string>
#include <unordered_set>

#include "Types.h"
#include "FunctionPrototypes.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! CompFlow system of PDEs problem: user defined
class CompFlowProblemUserDefined {

  private:
    using ncomp_t = tk::ctr::ncomp_type;
    using eq = tag::compflow;
    static constexpr ncomp_t m_ncomp = 5;    //!< Number of scalar components

  public:
    //! Evaluate initial condition solution at (x,y,z,t) for all components
    //! \param[in] system Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \param[in] t Physical time at which to evaluate the solution
    //! \return Values of all components evaluated at (x,y,z,t)
    //! \note The function signature must follow tk::SolutionFn
    static tk::SolutionFn::result_type
    solution( [[maybe_unused]] ncomp_t system,
              [[maybe_unused]] ncomp_t ncomp,
              [[maybe_unused]] tk::real x,
              [[maybe_unused]] tk::real y,
              [[maybe_unused]] tk::real z,
              [[maybe_unused]] tk::real t )
    {
      Assert( ncomp == m_ncomp, "Number of scalar components must be " +
                                std::to_string(m_ncomp) );
      return {{ 1.0, 0.0, 0.0, 1.0, 293.0 }};
    }

    //! \brief Evaluate the increment from t to t+dt of the analytical solution
    //!   at (x,y,z) for all components
    //! \return Increment in values of all components: all zero for now
    static std::array< tk::real, 5 >
    solinc( ncomp_t, tk::real, tk::real, tk::real, tk::real, tk::real ) {
      return {{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
    }

    //! Compute and return source term for Rayleigh-Taylor manufactured solution
    //! \details No-op for user-deefined problems.
    //! \return Array of reals containing the source for all components
    //! \note The function signature must follow tk::SrcFn
    static tk::SrcFn::result_type
    src( ncomp_t, ncomp_t, tk::real, tk::real, tk::real, tk::real )
    { return {{ 0.0, 0.0, 0.0, 0.0, 0.0 }}; }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    static std::vector< std::string > fieldNames( ncomp_t ) {
      std::vector< std::string > n;
      n.push_back( "density" );
      n.push_back( "x-velocity" );
      n.push_back( "y-velocity" );
      n.push_back( "z-velocity" );
      n.push_back( "specific total energy" );
      n.push_back( "pressure" );
      n.push_back( "temperature" );
      return n;
    }

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    static void side( std::unordered_set< int >& conf ) {
      using tag::param; using tag::bcdir;
      for (const auto& s : g_inputdeck.get< param, eq, bcdir >())
        conf.insert( std::stoi(s[0]) );
    }

    //! Return field output going to file
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    static std::vector< std::vector< tk::real > >
    fieldOutput( ncomp_t,
                 ncomp_t,
                 ncomp_t offset,
                 tk::real,
                 tk::real,
                 const std::vector< tk::real >&,
                 const std::array< std::vector< tk::real >, 3 >&,
                 tk::Fields& U )
    {
      std::vector< std::vector< tk::real > > out;
      const auto r = U.extract( 0, offset );
      const auto ru = U.extract( 1, offset );
      const auto rv = U.extract( 2, offset );
      const auto rw = U.extract( 3, offset );
      const auto re = U.extract( 4, offset );
      out.push_back( r );
      std::vector< tk::real > u = ru;
      std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( u );
      std::vector< tk::real > v = rv;
      std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( v );
      std::vector< tk::real > w = rw;
      std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( w );
      std::vector< tk::real > E = re;
      std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( E );
      std::vector< tk::real > p = r;
      tk::real g = g_inputdeck.get< tag::param, eq, tag::gamma >()[0];
      for (std::size_t i=0; i<p.size(); ++i)
        p[i] = (g-1.0)*r[i]*(E[i] - (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0);
      out.push_back( p );
      std::vector< tk::real > T = r;
      tk::real cv = g_inputdeck.get< tag::param, eq, tag::cv >()[0];
      for (std::size_t i=0; i<T.size(); ++i)
        T[i] = cv*(E[i] - (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0);
      out.push_back( T );
      return out;
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    static std::vector< std::string > names( ncomp_t )
    { return { "r", "ru", "rv", "rw", "re" }; }

   static ctr::ProblemType type() noexcept
   { return ctr::ProblemType::USER_DEFINED; }
};
} // inciter::

#endif // CompFlowProblemUserDefined_h
