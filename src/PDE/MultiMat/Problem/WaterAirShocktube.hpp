// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/WaterAirShocktube.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for Water-Air shock-tube
  \details   This file defines a policy class for the multi-material
    compressible flow equations, defined in PDE/MultiMat/MultiMat.hpp.
    See PDE/MultiMat/Problem.hpp for general requirements on Problem policy
    classes for MultiMat.
*/
// *****************************************************************************
#ifndef MultiMatProblemWaterAirShocktube_h
#define MultiMatProblemWaterAirShocktube_h

#include <string>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"

namespace inciter {

//! MultiMat system of PDEs problem: Water-Air shock-tube
//! \see Chiapolino, A., Saurel, R., & Nkonga, B. (2017). Sharpening diffuse
//!   interfaces with compressible fluids on unstructured meshes. Journal of
//!   Computational Physics, 340, 389-417.
class MultiMatProblemWaterAirShocktube {

  protected:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::multimat;

  public:
    //! Initialize numerical solution
    static tk::InitializeFn::result_type
    initialize( ncomp_t system, ncomp_t ncomp, tk::real x, tk::real, tk::real,
                tk::real );

    static std::vector< tk::real >
    analyticSolution( ncomp_t system, ncomp_t ncomp,  tk::real x, tk::real y,
                      tk::real z, tk::real t )
    { return initialize( system, ncomp, x, y, z, t ); }

    //! Compute and return source term for Rayleigh-Taylor manufactured solution
    //! Compute and return source term for this problem
    static tk::MultiMatSrcFn::result_type
    src( ncomp_t, ncomp_t ncomp, tk::real, tk::real, tk::real, tk::real )
    { std::vector< tk::real > s( ncomp, 0.0 ); }

    //! Return names of integral variables to be output to diagnostics file
    static std::vector< std::string > names( ncomp_t );

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::WATERAIR_SHOCKTUBE; }
};

} // inciter::

#endif // MultiMatProblemWaterAirShocktube_h
