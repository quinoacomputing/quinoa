// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/BoxInitialization.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     User-defined box initialization
  \details   This file defines functions for initializing solutions for
    compressible single-material equations inside the user-defined box.
*/
// *****************************************************************************
#ifndef CompFlowBoxInitialization_h
#define CompFlowBoxInitialization_h

#include "Fields.hpp"
#include "EoS/EoS.hpp"
#include "Control/Inciter/Types.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

//! Set the solution in the user-defined IC box
void initializeBox( std::size_t system,
                    tk::real VRatio,
                    tk::real t,
                    const inciter::ctr::icbox& icbox,
                    tk::real bgpreic,
                    tk::real cv,
                    std::vector< tk::real >& s );

} //inciter::

#endif // CompFlowBoxInitialization_h
