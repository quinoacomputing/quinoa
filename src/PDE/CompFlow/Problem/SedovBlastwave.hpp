// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/SedovBlastwave.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for Sedov's blastwave
  \details   This file defines a policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problem.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowProblemSedovBlastwave_h
#define CompFlowProblemSedovBlastwave_h

#include <string>
#include <unordered_set>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"

namespace inciter {

//! CompFlow system of PDEs problem: Sedov blast-wave
class CompFlowProblemSedovBlastwave {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::compflow;

  public:
    //! Evaluate analytical solution at (x,y,0) for all components
    static tk::SolutionFn::result_type
    solution( ncomp_t system, ncomp_t ncomp,
              tk::real x, tk::real y, tk::real, tk::real );

    //! \brief Evaluate the increment from t to t+dt of the analytical solution
    //!   at (x,y,z) for all components
    std::vector< tk::real >
    solinc( ncomp_t system, ncomp_t ncomp,
            tk::real x, tk::real y, tk::real z, tk::real t, tk::real dt ) const;

    //! Compute and return source term for this problem
    static tk::SrcFn::result_type
    src( ncomp_t, ncomp_t, tk::real, tk::real, tk::real, tk::real );

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    void side( std::unordered_set< int >& conf ) const;

    //! Return field names to be output to file
    std::vector< std::string > fieldNames( ncomp_t ) const;

    //! Return field output going to file
    std::vector< std::vector< tk::real > >
    fieldOutput( ncomp_t system,
                 ncomp_t /*ncomp*/,
                 ncomp_t offset,
                 tk::real,
                 tk::real /*V*/,
                 const std::vector< tk::real >& /*vol*/,
                 const std::array< std::vector< tk::real >, 3 >& /*coord*/,
                 tk::Fields& U ) const;
 
    //! Return names of integral variables to be output to diagnostics file
    std::vector< std::string > names( ncomp_t ) const;

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SEDOV_BLASTWAVE; }
};

} // inciter::

#endif // CompFlowProblemSedovBlasttwave_h
