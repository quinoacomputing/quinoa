// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/SedovBlastwave.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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
#include "EoS/EoS_Base.hpp"

namespace inciter {

//! CompFlow system of PDEs problem: Sedov blast-wave
class CompFlowProblemSedovBlastwave {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::compflow;

  public:
    //! Initialize numerical solution
    static tk::InitializeFn::result_type
    initialize( ncomp_t system, ncomp_t, tk::real x, tk::real y,
                tk::real z, tk::real );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static tk::InitializeFn::result_type
    analyticSolution( ncomp_t system, ncomp_t, 
                      std::vector< EoS_Base* >, tk::real x, tk::real y,
                      tk::real z, tk::real );

    //! Compute and return source term for this problem
    //! \param[in,out] r Density source
    //! \param[in,out] ru X momentum source
    //! \param[in,out] rv Y momentum source
    //! \param[in,out] rw Z momentum source
    //! \param[in,out] re Specific total energy source
    //! \note The function signature must follow tk::SrcFn
    static tk::CompFlowSrcFn::result_type
    src( ncomp_t, tk::real, tk::real, tk::real, tk::real,
         tk::real& r, tk::real& ru, tk::real& rv, tk::real& rw, tk::real& re )
    { r = ru = rv = rw = re = 0.0; }

    //! Return analytic field names to be output to file
    std::vector< std::string > analyticFieldNames( ncomp_t ) const;

    //! Return names of integral variables to be output to diagnostics file
    std::vector< std::string > names( ncomp_t ) const;

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SEDOV_BLASTWAVE; }
};

} // inciter::

#endif // CompFlowProblemSedovBlasttwave_h
