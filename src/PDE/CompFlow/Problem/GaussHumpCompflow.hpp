// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/GaussHumpCompflow.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file declares a problem policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problem.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowProblemGaussHump_h
#define CompFlowProblemGaussHump_h

#include <string>
#include <unordered_set>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"
#include "EoS/EosVariant.hpp"

namespace inciter {

//! CompFlow system of PDEs problem: GaussHump
class CompFlowProblemGaussHump {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::compflow;

  public:
    //! Initialize numerical solution
    static tk::InitializeFn::result_type
    initialize( ncomp_t system, ncomp_t ncomp, const std::vector< EOS >&,
                tk::real x, tk::real y, tk::real, tk::real t );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static tk::InitializeFn::result_type
    analyticSolution( ncomp_t system, ncomp_t ncomp,
                      const std::vector< EOS >&, tk::real x, tk::real y,
                      tk::real, tk::real t );

    //! Compute and return source term for manufactured solution
    //! \param[in,out] sv Source term vector
    //! \note The function signature must follow tk::SrcFn
    static tk::SrcFn::result_type
    src( ncomp_t, ncomp_t, const std::vector< EOS >&, tk::real, tk::real,
         tk::real, tk::real, std::vector< tk::real >& sv )
    {
      Assert(sv.size() == 5, "Incorrect source vector size");
      sv[0] = sv[1] = sv[2] = sv[3] = sv[4] = 0.0;
    }

    //! Return analytic field names to be output to file
    std::vector< std::string > analyticFieldNames( ncomp_t ) const;

    //! Return names of integral variables to be output to diagnostics file
    std::vector< std::string > names( ncomp_t ) const;

    //! Assign prescribed velocity at a point
    static std::vector< std::array< tk::real, 3 > >
    prescribedVelocity( ncomp_t,
                        ncomp_t ncomp,
                        tk::real,
                        tk::real,
                        tk::real );

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::GAUSS_HUMP_COMPFLOW; }
};

} // inciter::

#endif // CompFlowProblemGaussHump_h
