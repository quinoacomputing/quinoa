// *****************************************************************************
/*!
  \file      src/Control/Inciter/Components.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Storage for number of components
  \details   Storage for number of components. This is part of the input deck
     stack and is thus populated during control file parsing.
*/
// *****************************************************************************
#ifndef InciterComponents_h
#define InciterComponents_h

#include "SystemComponents.hpp"

namespace inciter {
namespace ctr {

//! Number of components storage for all systems of equations supported
using ncomps = tk::ctr::ncomponents<
   tag::transport
 , tag::compflow
 , tag::multimat
>;

} // ctr::
} // inciter::

#endif // InciterComponents_h
