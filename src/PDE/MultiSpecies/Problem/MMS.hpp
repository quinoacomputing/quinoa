// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Problem/MMS.hpp
  \brief     Three-dimensional manufactured solution for multi-species flow
  \details   Problem policy for the single-species Danis-Yan traveling wave.
    Requires constant viscosity mu=1e-3 and a calorically perfect gas with
    R=0.4, cp=1.4, and zero internal-energy offset, so that e=T and gamma=1.4.
    The source uses Pr=0.71, matching the multi-species viscous discretization.
*/
// *****************************************************************************
#ifndef MultiSpeciesProblemMMS_h
#define MultiSpeciesProblemMMS_h

#include <vector>

#include "Types.hpp"
#include "Tags.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Problem.hpp"
#include "EoS/EOS.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"

namespace inciter {

//! MultiSpecies system of PDEs problem: manufactured traveling wave
class MultiSpeciesProblemMMS {

  private:
    using ncomp_t = tk::ncomp_t;
    using eq = tag::multispecies;

  public:
    //! Evaluate the exact conserved state at (x,y,z,t)
    static tk::InitializeFn::result_type
    initialize( ncomp_t ncomp, const std::vector< EOS >& mat_blk,
                tk::real x, tk::real y, tk::real z, tk::real t );

    //! Evaluate the exact conserved state with temperature appended
    static std::vector< tk::real >
    analyticSolution( ncomp_t ncomp, const std::vector< EOS >& mat_blk,
                      tk::real x, tk::real y, tk::real z, tk::real t );

    //! Evaluate the manufactured source for all conserved components
    static tk::SrcFn::result_type
    src( ncomp_t, const std::vector< EOS >& mat_blk,
         tk::real x, tk::real y, tk::real z, tk::real t,
         std::vector< tk::real >& sv );

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::MMS; }
};

} // inciter::

#endif // MultiSpeciesProblemMMS_h
