// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/UserDefined.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the multi-material compressible flow
    equations
  \details   This file defines a Problem policy class for the multi-material
    compressible flow equations, defined under PDE/MultiMat. See
    PDE/MultiMat/Problem.h for general requirements on Problem policy
    classes for MultiMat.
*/
// *****************************************************************************
#ifndef MultiMatProblemUserDefined_h
#define MultiMatProblemUserDefined_h

#include <string>

#include "Types.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Problem.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! MultiMat system of PDEs problem: user defined
class MultiMatProblemUserDefined {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::multimat;

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
              [[maybe_unused]] tk::real t,
              int& )
    {
      Assert( ncomp == ncomp, "Number of scalar components must be " +
                              std::to_string(ncomp) );
      return {{ 1.0, 0.0, 0.0, 1.0, 293.0 }};
    }

    //! Compute and return source term for Rayleigh-Taylor manufactured solution
    //! \details No-op for user-deefined problems.
    static tk::MultiMatSrcFn::result_type
    src( ncomp_t, ncomp_t ncomp, tk::real, tk::real, tk::real, tk::real )
    { std::vector< tk::real > s( ncomp, 0.0 ); }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    static std::vector< std::string > names( ncomp_t )
    { return { "r", "ru", "rv", "rw", "re" }; }

   static ctr::ProblemType type() noexcept
   { return ctr::ProblemType::USER_DEFINED; }
};
} // inciter::

#endif // MultiMatProblemUserDefined_h
