// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/EquilInterfaceAdvect.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for equilibrium interface advection
  \details   This file defines a policy class for the multi-material
    compressible flow equations, defined in PDE/MultiMat/MultiMat.hpp.
    See PDE/MultiMat/Problem.hpp for general requirements on Problem policy
    classes for MultiMat.
*/
// *****************************************************************************
#ifndef MultiMatProblemEquilInterfaceAdvect_h
#define MultiMatProblemEquilInterfaceAdvect_h

#include <string>

#include "Types.hpp"
#include "Fields.hpp"
#include "Vector.hpp"
#include "FunctionPrototypes.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "EoS/EOS.hpp"

namespace inciter {

//! MultiMat system of PDEs problem: equilibrium interface advection
class MultiMatProblemEquilInterfaceAdvect {

  protected:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::multimat;

  public:
    //! Initialize numerical solution
    static tk::InitializeFn::result_type
    initialize( ncomp_t ncomp, const std::vector< EOS >&,
                tk::real x, tk::real, tk::real, tk::real );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static std::vector< tk::real >
    analyticSolution( ncomp_t ncomp,
                      const std::vector< EOS >& mat_blk, tk::real x,
                      tk::real y, tk::real z, tk::real t )
    { return initialize( ncomp, mat_blk, x, y, z, t ); }

    //! Compute and return source term for this problem
    static tk::SrcFn::result_type
    src( ncomp_t nmat, const std::vector< EOS >& mat_blk,
         tk::real x, tk::real y, tk::real z, tk::real t,
         std::vector< tk::real >& sv )
    {
      auto ncomp = 3*nmat+3;
      Assert(sv.size() == ncomp, "Incorrect source vector size");

      // solution at given location and time
      auto s = initialize(ncomp, mat_blk, x, y, z, t);
      tk::real rhob(0.0);
      for (std::size_t k=0; k<nmat; ++k) {
        rhob += s[densityIdx(nmat,k)];
      }
      std::array< tk::real, 3 > u0 {{s[momentumIdx(nmat,0)]/rhob,
        s[momentumIdx(nmat,1)]/rhob, s[momentumIdx(nmat,2)]/rhob}};
      tk::real umag2 = tk::dot(u0, u0);

      // source terms
      for (std::size_t i=0; i<3; ++i) {
        sv[momentumIdx(nmat,i)] = 0.0;
      }
      for (std::size_t k=0; k<nmat; ++k) {
        sv[volfracIdx(nmat,k)] = 0.0;
        sv[densityIdx(nmat,k)] = (u0[0]+u0[1]+u0[2]) * s[volfracIdx(nmat,k)];
        for (std::size_t i=0; i<3; ++i) {
          sv[momentumIdx(nmat,i)] += u0[i] * sv[densityIdx(nmat,k)];
        }
        sv[energyIdx(nmat,k)] = 0.5 * umag2 * sv[densityIdx(nmat,k)];
      }
    }

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::EQUILINTERFACE_ADVECT; }
};

} // inciter::

#endif // MultiMatProblemEquilInterfaceAdvect_h
