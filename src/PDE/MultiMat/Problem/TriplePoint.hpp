// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/TriplePoint.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for triple point
  \details   This file defines a policy class for the multi-material
    compressible flow equations, defined in PDE/MultiMat/MultiMat.hpp.
    See PDE/MultiMat/Problem.hpp for general requirements on Problem policy
    classes for MultiMat.
*/
// *****************************************************************************
#ifndef MultiMatProblemTriplePoint_h
#define MultiMatProblemTriplePoint_h

#include <string>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"

namespace inciter {

//! MultiMat system of PDEs problem: Triple point problem
//! \see Galera, S., Maire, P. H., & Breil, J. (2010). A two-dimensional
//!   unstructured cell-centered multi-material ALE scheme using VOF interface
//!   reconstruction. Journal of Computational Physics, 229(16), 5755-5787.
class MultiMatProblemTriplePoint {

  protected:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::multimat;

  public:
    //! Evaluate analytical solution at (x,y,0) for all components
    static tk::SolutionFn::result_type
    solution( ncomp_t system, ncomp_t ncomp, tk::real x, tk::real y, tk::real,
              tk::real, int& );

    //! Compute and return source term for this problem
    static tk::SrcFn::result_type
    src( ncomp_t, ncomp_t, tk::real, tk::real, tk::real, tk::real );

    //! Return field names to be output to file
    static std::vector< std::string > fieldNames( ncomp_t );

    //! Return field output going to file
    static std::vector< std::vector< tk::real > >
    fieldOutput( ncomp_t system,
                 ncomp_t ncomp,
                 ncomp_t offset,
                 std::size_t,
                 tk::real,
                 tk::real /*V*/,
                 const std::vector< tk::real >& /*vol*/,
                 const std::array< std::vector< tk::real >, 3 >& /*coord*/,
                 tk::Fields& U,
                 const tk::Fields& P );

    //! Return names of integral variables to be output to diagnostics file
    static std::vector< std::string > names( ncomp_t );

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::TRIPLE_POINT; }
};

} // inciter::

#endif // MultiMatProblemTriplePoint_h
