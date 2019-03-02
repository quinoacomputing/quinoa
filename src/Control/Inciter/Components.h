// *****************************************************************************
/*!
  \file      src/Control/Inciter/Components.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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
  tag::transport,             std::vector< tk::ctr::ncomp_type >,
  tag::compflow,              std::vector< tk::ctr::ncomp_type >,
  tag::multimat,              std::vector< tk::ctr::ncomp_type >
>;

} // ctr::
} // inciter::

#endif // InciterComponents_h
