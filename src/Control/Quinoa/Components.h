//******************************************************************************
/*!
  \file      src/Control/Quinoa/Components.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 10:50:24 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Storage for number of components
  \details   Storage for number of components
*/
//******************************************************************************
#ifndef QuinoaComponents_h
#define QuinoaComponents_h

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/remove.hpp>
#include <boost/mpl/at.hpp>

#include <make_list.h>

#include <Components.h>

namespace quinoa {
namespace ctr {

//! Number of components of models and equations
using ncomps = tk::ctr::ncomponents<
  tag::position,     std::vector< unsigned int >, //!< Position models
  tag::mass,         std::vector< unsigned int >, //!< Mass models
  tag::hydro,        std::vector< unsigned int >, //!< Hydro models
  tag::mix,          std::vector< unsigned int >, //!< Material mix models
  tag::frequency,    std::vector< unsigned int >  //!< Turbulent frequency models
>;

} // ctr::
} // quinoa::

#endif // QuinoaComponents_h
