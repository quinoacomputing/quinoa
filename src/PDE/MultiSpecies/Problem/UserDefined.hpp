// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Problem/UserDefined.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the multi-species compressible flow
    equations
  \details   This file defines a Problem policy class for the multi-species
    compressible flow equations, defined under PDE/MultiSpecies. See
    PDE/MultiSpecies/Problem.h for general requirements on Problem policy
    classes for MultiSpecies.
*/
// *****************************************************************************
#ifndef MultiSpeciesProblemUserDefined_h
#define MultiSpeciesProblemUserDefined_h

#include <string>

#include "Types.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Problem.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"
#include "EoS/EOS.hpp"

namespace inciter {

//! MultiSpecies system of PDEs problem: user defined
class MultiSpeciesProblemUserDefined {

  private:
    using ncomp_t = tk::ncomp_t;
    using eq = tag::multispecies;

  public:
    //! Initialize numerical solution
    static tk::InitializeFn::result_type
    initialize( ncomp_t ncomp, const std::vector< EOS >&,
                tk::real, tk::real, tk::real, tk::real );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static std::vector< tk::real >
    analyticSolution( ncomp_t ncomp,
                      const std::vector< EOS >& mat_blk, tk::real x,
                      tk::real y, tk::real z, tk::real t )
    { return initialize( ncomp, mat_blk, x, y, z, t ); }

    //! Compute and return source term for Rayleigh-Taylor manufactured solution
    //! \details No-op for user-deefined problems.
    static tk::SrcFn::result_type
    src( ncomp_t, const std::vector< EOS >&,tk::real, tk::real,
         tk::real, tk::real, std::vector< tk::real >& sv )
    {
      for (std::size_t i=0; i<sv.size(); ++i) {
        sv[i] = 0.0;
      }
    }

   static ctr::ProblemType type() noexcept
   { return ctr::ProblemType::USER_DEFINED; }
};
} // inciter::

#endif // MultiSpeciesProblemUserDefined_h
