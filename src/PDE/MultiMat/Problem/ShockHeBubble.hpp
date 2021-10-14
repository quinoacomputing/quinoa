// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/ShockHeBubble.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for gas impact
  \details   This file defines a policy class for the multi-material
    compressible flow equations, defined in PDE/MultiMat/DGMultiMat.hpp.
    See PDE/MultiMat/Problem.hpp for general requirements on Problem policy
    classes for MultiMat.
*/
// *****************************************************************************
#ifndef MultiMatProblemShockHeBubble_h
#define MultiMatProblemShockHeBubble_h

#include <string>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"

namespace inciter {

//! MultiMat system of PDEs problem: Shock He-bubble interaction
//! \see Quirk, J. J., & Karni, S. (1996). On the dynamics of a shock–bubble         
//!   interaction. Journal of Fluid Mechanics, 318, 129-163.
class MultiMatProblemShockHeBubble {

  protected:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::multimat;

  public:
    //! Initialize numerical solution
    static tk::InitializeFn::result_type
    initialize( ncomp_t system, ncomp_t ncomp, tk::real x, tk::real y,
                tk::real z, tk::real );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static std::vector< tk::real >
    analyticSolution( ncomp_t system, ncomp_t ncomp,  tk::real x, tk::real y,
                      tk::real z, tk::real t )
    { return initialize( system, ncomp, x, y, z, t ); }

    //! Compute and return source term for this problem
    static tk::MultiMatSrcFn::result_type
    src( ncomp_t, tk::real, tk::real, tk::real, tk::real,
         tk::real& r, tk::real& ru, tk::real& rv, tk::real& rw, tk::real& re )
    { r = ru = rv = rw = re = 0.0; }

    //! Return names of integral variables to be output to diagnostics file
    static std::vector< std::string > names( ncomp_t );

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SHOCK_HEBUBBLE; }
};

} // inciter::

#endif // MultiMatProblemShockHeBubble_h
