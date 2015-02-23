//******************************************************************************
/*!
  \file      src/Control/Inciter/Components.h
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 08:01:45 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Storage for number of components
  \details   Storage for number of components. This is part of the input deck
     stack and is thus populated during the control file parsing.
*/
//******************************************************************************
#ifndef InciterComponents_h
#define InciterComponents_h

#include <Components.h>

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
