//******************************************************************************
/*!
  \file      src/Walker/Integrator.C
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 01:15:10 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Integrator advances differential equations
  \details   Integrator advances differential equations. There are a potentially
    large number of Integrator Charm++ chares created by Distributor. Each
    integrator gets a chunk of the full load and does the same: initializes and
    advances multiple ordinary or stochastic differential equations in time.
    Note that there is no spatial dependence, these equations describe spatially
    homogeneous processes.
*/
//******************************************************************************

#include "Integrator.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "integrator.def.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

