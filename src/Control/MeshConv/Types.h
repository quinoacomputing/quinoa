//******************************************************************************
/*!
  \file      src/Control/MeshConv/Types.h
  \author    J. Bakosi
  \date      Tue 08 Apr 2014 09:04:37 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for MeshConv's parsers
  \details   Types for MeshConv's parsers
*/
//******************************************************************************
#ifndef MeshConvTypes_h
#define MeshConvTypes_h

#include <pegtl.hh>

#include <MeshConv/Tags.h>

namespace meshconv {
//! control and parsing
namespace ctr {

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::input,       std::string,  //!< Input filename
  tag::output,      std::string   //!< Output filename
>;

//! PEGTL location type to use throughout all of MeshConv's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // meshconv::

#endif // MeshConvTypes_h
