// *****************************************************************************
/*!
  \file      src/PDE/EoS/EoS_Base.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Equation of state base class
  \details   This file defines virtual functions for equations of state for the
    compressible flow equations.
*/
// *****************************************************************************
#ifndef EoS_Base_h
#define EoS_Base_h

#include <cmath>
#include "Data.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

class EoS_Base {
  public:
    virtual tk::real eos_pressure( ncomp_t,
                                   tk::real,
                                   tk::real,
                                   tk::real,
                                   tk::real,
                                   tk::real,
                                   tk::real=1.0,
                                   std::size_t=0 )=0;

    // Virtual destructor
    virtual ~EoS_Base(){}
};

} //inciter::

#endif // EoS_Base_h
