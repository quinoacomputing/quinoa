// *****************************************************************************
/*!
  \file      src/Control/Inciter/Components.h
  \author    J. Bakosi
  \date      Wed 06 Jul 2016 11:51:01 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Storage for number of components
  \details   Storage for number of components. This is part of the input deck
     stack and is thus populated during the control file parsing.
*/
// *****************************************************************************
#ifndef InciterComponents_h
#define InciterComponents_h

#include "SystemComponents.h"

namespace inciter {
namespace ctr {

//! Number of components of partial differential equations
using ncomps = tk::ctr::ncomponents<
  tag::advdiff,      std::vector< tk::ctr::ncomp_type >,
  tag::poisson,      std::vector< tk::ctr::ncomp_type >,
  tag::euler,        std::vector< tk::ctr::ncomp_type >
>;

} // ctr::
} // inciter::

#endif // InciterComponents_h
