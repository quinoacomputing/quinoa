// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Problem/MixingLayer.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for mixing layer
  \details   This file defines a policy class for the multi-species
    compressible flow equations, defined in PDE/MultiSpecies/MultiSpecies.hpp.
    See PDE/MultiSpecies/Problem.hpp for general requirements on Problem policy
    classes for MultiMat.
*/
// *****************************************************************************
#ifndef MultiSpeciesProblemMixingLayer_h
#define MultiSpeciesProblemMixingLayer_h

#include <string>

#include "Types.hpp"
#include "Fields.hpp"
#include "Vector.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Problem.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"
#include "EoS/EOS.hpp"


namespace inciter {

//! MultiSpecies system of PDEs problem: mixing layer
class MultiSpeciesProblemMixingLayer {

  protected:
    using ncomp_t = tk::ncomp_t;
    using eq = tag::multispecies;

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
    { return initialize( ncomp, mat_blk, x, y, z, t); }

    //! Compute and return source term for this problem
    static tk::SrcFn::result_type
    src( ncomp_t, const std::vector< EOS >&,
         tk::real, tk::real, tk::real, tk::real,
         std::vector< tk::real >& sv )
    {
      for (std::size_t i=0; i<sv.size(); ++i) {
        sv[i] = 0.0;
      }
    }

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::MIXING_LAYER; }
};

} // inciter::

#endif // MultiSpeciesProblemMixingLayer_h
