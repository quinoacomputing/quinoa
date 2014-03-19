//******************************************************************************
/*!
  \file      src/Control/Gmsh2Exo/Types.h
  \author    J. Bakosi
  \date      Wed Mar 19 10:13:36 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for Gmsh2Exo's parsers
  \details   Types for Gmsh2Exo's parsers
*/
//******************************************************************************
#ifndef Gmsh2ExoTypes_h
#define Gmsh2ExoTypes_h

#include <Gmsh2Exo/Tags.h>

namespace gmsh2exo {
//! control and parsing
namespace ctr {

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::control,   std::string                 //!< Control filename
>;

//! PEGTL location type to use throughout all of Gmsh2Exo's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // gmsh2exo::

#endif // Gmsh2ExoTypes_h
