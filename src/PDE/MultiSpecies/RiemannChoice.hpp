// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/RiemannChoice.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register available Riemann solvers for multi-species compressible
             fluid dynamics.
  \details   Register available Riemann solvers for multi-species compressible
             fluid dynamics.
*/
// *****************************************************************************
#ifndef RiemannChoice_h
#define RiemannChoice_h

#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "Riemann/AUSMMultiSpecies.hpp"

namespace inciter {

  //! Get the Riemann solver function according to control file setup
  //! \param[in] flux Riemann solver from input deck
  //! \return Function pointer to the Riemann solver, must be of type
  //!   tk::RiemannFluxFn
  const static tk::RiemannFluxFn multispeciesRiemannSolver(ctr::FluxType flux)
  {
    tk::RiemannFluxFn fluxfn;

    if (flux == ctr::FluxType::AUSM) {
      fluxfn = AUSMMultiSpecies::flux;
    }
    else {
      Throw("Riemann solver not set up for multi-material PDEs.");
    }

    return fluxfn;
  }

} // inciter::

#endif // RiemannChoice_h
