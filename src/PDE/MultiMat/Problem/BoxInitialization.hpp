// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/BoxInitialization.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     User-defined box initialization
  \details   This file defines functions for initializing solutions for
    compressible multi-material equations inside the user-defined box.
*/
// *****************************************************************************
#ifndef MultiMatBoxInitialization_h
#define MultiMatBoxInitialization_h

#include "Fields.hpp"
#include "EoS/EoS.hpp"
#include "EoS/StiffenedGas.hpp"
#include "Control/Inciter/Types.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

//! Set the solution in the user-defined IC box
void initializeBox( std::size_t system,
                    tk::real VRatio,
                    tk::real,
                    const inciter::ctr::box& icbox,
                    std::vector< tk::real >& s );

} //inciter::

#endif // MultiMatBoxInitialization_h
