//******************************************************************************
/*!
  \file      src/Control/Inciter/Components.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:53:06 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Storage for number of components
  \details   Storage for number of components. This is part of the input deck
     stack and is thus populated during the control file parsing.
*/
//******************************************************************************
#ifndef InciterComponents_h
#define InciterComponents_h

#include "SystemComponents.h"

namespace inciter {
namespace ctr {

//! Number of components of models and equations
using ncomps = tk::ctr::ncomponents<
  tag::position,     std::vector< tk::ctr::ncomp_type >,
  tag::mass,         std::vector< tk::ctr::ncomp_type >,
  tag::hydro,        std::vector< tk::ctr::ncomp_type >,
  tag::mix,          std::vector< tk::ctr::ncomp_type >,
  tag::frequency,    std::vector< tk::ctr::ncomp_type >
>;

} // ctr::
} // inciter::

#endif // InciterComponents_h
