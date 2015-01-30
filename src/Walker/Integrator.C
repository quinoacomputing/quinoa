//******************************************************************************
/*!
  \file      src/Walker/Integrator.C
  \author    J. Bakosi
  \date      Thu 29 Jan 2015 09:32:38 AM MST
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

#include <Integrator.h>

#include <integrator.def.h>
