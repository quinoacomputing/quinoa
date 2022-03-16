// *****************************************************************************
/*!
  \file      src/PDE/EoS/StiffenedGas.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Stiffened-gas equation of state
  \details   This file defines virtual functions for equations of state for the
    compressible flow equations.
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
    // Constructor
    StiffenedGas(tk::real x, tk::real y, std::size_t imat) : gamma(x), pstiff(y) {
    std::cout << "EOS - Stiffened Gas Initialization: material= " << imat << ", gamma= " << gamma << ", pstiff= " << pstiff << std::endl;
    }

    tk::real eos_pressure( ncomp_t system,
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
//    if (!std::isfinite(partpressure)) {
//      std::cout << "Material-id:      " << imat << std::endl;
//      std::cout << "Volume-fraction:  " << alpha << std::endl;
//      std::cout << "Partial density:  " << arho << std::endl;
//      std::cout << "Total energy:     " << arhoE << std::endl;
//      std::cout << "Velocity:         " << u << ", " << v << ", " << w
//        << std::endl;
//      Throw("Material-" + std::to_string(imat) + " has nan/inf partial pressure: "
//        + std::to_string(partpressure) + ", material volume fraction: " +
//        std::to_string(alpha));
//    }

    return partpressure;
  }
};

} //inciter::

#endif // StiffenedGas_h

