// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/UnderwaterEx.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for underwater explosion
  \details   This file defines a policy class for the multi-material
    compressible flow equations, defined in PDE/MultiMat/DGMultiMat.hpp.
    See PDE/MultiMat/Problem.hpp for general requirements on Problem policy
    classes for MultiMat.
*/
// *****************************************************************************
#ifndef MultiMatProblemUnderwaterEx_h
#define MultiMatProblemUnderwaterEx_h

#include <string>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"

namespace inciter {

//! MultiMat system of PDEs problem: Underwater explosion problem
//! \see Chiapolino, A., Saurel, R., & Nkonga, B. (2017). Sharpening diffuse
//!   interfaces with compressible fluids on unstructured meshes. Journal of
//!   Computational Physics, 340, 389-417.
class MultiMatProblemUnderwaterEx {

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
    { return ctr::ProblemType::UNDERWATER_EX; }
};

} // inciter::

#endif // MultiMatProblemUnderwaterEx_h
