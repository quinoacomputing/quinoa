// *****************************************************************************
/*!
  \file      src/PDE/EoS/StiffenedGas.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Stiffened-gas equation of state
  \details   This file defines functions for the stiffened gas equation of
             state for the compressible flow equations.
*/
// *****************************************************************************
#ifndef StiffenedGas_h
#define StiffenedGas_h

#include <cmath>
#include <iostream>
#include "Data.hpp"
#include "EoS_Base.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

class StiffenedGas: public EoS_Base {

  private:
    tk::real gamma, pstiff;

  public:
    // *****************************************************************************
    //  Constructor
    //! \param[in] gamma Ratio of specific heats
    //! \param[in] pstiff Stiffened pressure term
    //! \param[in] size_t Material index
    // *****************************************************************************
    StiffenedGas(tk::real x, tk::real y, std::size_t ) : gamma(x), pstiff(y)
    { }

    //! \brief Calculate pressure from the material density, momentum and total
    //!   energy using the stiffened-gas equation of state
    //! \param[in] system Equation system index
    //! \param[in] arho Material partial density (alpha_k * rho_k)
    //! \param[in] u X-velocity
    //! \param[in] v Y-velocity
    //! \param[in] w Z-velocity
    //! \param[in] arhoE Material total energy (alpha_k * rho_k * E_k)
    //! \param[in] alpha Material volume fraction. Default is 1.0, so that for the
    //!   single-material system, this argument can be left unspecified by the
    //!   calling code
    //! \param[in] imat Material-id who's EoS is required. Default is 0, so that
    //!   for the single-material system, this argument can be left unspecified by
    //!   the calling code
    //! \return Material partial pressure (alpha_k * p_k) calculated using the
    //!   stiffened-gas EoS
    tk::real eos_pressure( ncomp_t,
                           tk::real arho,
                           tk::real u,
                           tk::real v,
                           tk::real w,
                           tk::real arhoE,
                           tk::real alpha=1.0,
                           std::size_t imat=0 )
  {
    tk::real g = gamma;
    tk::real p_c = pstiff;

    tk::real partpressure = (arhoE - 0.5 * arho * (u*u + v*v + w*w) - alpha*p_c)
                            * (g-1.0) - alpha*p_c;

    // check partial pressure divergence
    if (!std::isfinite(partpressure)) {
      std::cout << "Material-id:      " << imat << std::endl;
      std::cout << "Volume-fraction:  " << alpha << std::endl;
      std::cout << "Partial density:  " << arho << std::endl;
      std::cout << "Total energy:     " << arhoE << std::endl;
      std::cout << "Velocity:         " << u << ", " << v << ", " << w
        << std::endl;
      Throw("Material-" + std::to_string(imat) + " has nan/inf partial pressure: "
        + std::to_string(partpressure) + ", material volume fraction: " +
        std::to_string(alpha));
    }

    return partpressure;
  }
};

} //inciter::

#endif // StiffenedGas_h

