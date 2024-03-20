// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/InterfaceAdvection.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the multi-material compressible flow
    equations
  \details   This file defines a Problem policy class for the multi-material
    compressible flow equations, defined under PDE/MultiMat. See
    PDE/MultiMat/Problem.hpp for general requirements on Problem policy classes
    for MultiMat.
*/
// *****************************************************************************
#ifndef MultiMatProblemInterfaceAdvection_h
#define MultiMatProblemInterfaceAdvection_h

#include <string>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Problem.hpp"
#include "Inciter/InputDeck/New2InputDeck.hpp"
#include "EoS/EOS.hpp"

namespace inciter {

//! MultiMat system of PDEs problem: interface advection
//! \see Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
//!   equations with relevance to Inertial Confinement Fusion", Journal of
//!   Computational Physics 267 (2014) 196-209.
class MultiMatProblemInterfaceAdvection {

  private:
    using ncomp_t = tk::ncomp_t;
    using eq = newtag::multimat;

  public:
    //! Initialize numerical solution
    static tk::InitializeFn::result_type
    initialize( ncomp_t ncomp,
                const std::vector< EOS >&,
                tk::real x,
                tk::real y,
                tk::real /*z*/,
                tk::real t );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static std::vector< tk::real >
    analyticSolution( ncomp_t ncomp,
                      const std::vector< EOS >& mat_blk, tk::real x,
                      tk::real y, tk::real z, tk::real t )
    { return initialize( ncomp, mat_blk, x, y, z, t ); }

    //! Compute and return source term for interface advection
    static tk::SrcFn::result_type
    src( ncomp_t, const std::vector< EOS >&, tk::real, tk::real,
         tk::real, tk::real, std::vector< tk::real >& sv )
    {
      for (std::size_t i=0; i<sv.size(); ++i) {
        sv[i] = 0.0;
      }
    }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::INTERFACE_ADVECTION; }
};

} // inciter::

#endif // MultiMatProblemInterfaceAdvection_h
