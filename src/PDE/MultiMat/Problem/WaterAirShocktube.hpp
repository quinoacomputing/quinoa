// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/WaterAirShocktube.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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
#include <unordered_set>

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
    //! Evaluate analytical solution at (x,y,0) for all components
    static tk::SolutionFn::result_type
    solution( ncomp_t system, ncomp_t ncomp, tk::real x, tk::real, tk::real,
              tk::real );

    //! \brief Evaluate the increment from t to t+dt of the analytical solution
    //!   at (x,y,z) for all components
    static std::vector< tk::real >
    solinc( ncomp_t system, ncomp_t ncomp, tk::real x, tk::real y, tk::real z,
      tk::real t, tk::real dt );

    //! Compute and return source term for this problem
    static tk::SrcFn::result_type
    src( ncomp_t, ncomp_t, tk::real, tk::real, tk::real, tk::real );

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    static void side( std::unordered_set< int >& conf );

    //! Return field names to be output to file
    static std::vector< std::string > fieldNames( ncomp_t );

    //! Return field output going to file
    static std::vector< std::vector< tk::real > >
    fieldOutput( ncomp_t system,
                 ncomp_t ncomp,
                 ncomp_t offset,
                 tk::real,
                 tk::real /*V*/,
                 const std::vector< tk::real >& /*vol*/,
                 const std::array< std::vector< tk::real >, 3 >& /*coord*/,
                 tk::Fields& U );

    //! Return names of integral variables to be output to diagnostics file
    static std::vector< std::string > names( ncomp_t );

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::WATERAIR_SHOCKTUBE; }
};

} // inciter::

#endif // MultiMatProblemWaterAirShocktube_h
