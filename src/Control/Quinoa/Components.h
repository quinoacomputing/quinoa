//******************************************************************************
/*!
  \file      src/Control/Quinoa/Components.h
  \author    J. Bakosi
  \date      Wed 14 Jan 2015 10:49:00 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Storage for number of components
  \details   Storage for number of components. This is part of the input deck
     stack and is thus populated during the control file parsing.
*/
//******************************************************************************
#ifndef QuinoaComponents_h
#define QuinoaComponents_h

#include <Components.h>

namespace quinoa {
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
} // quinoa::

#endif // QuinoaComponents_h
