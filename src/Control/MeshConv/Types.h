//******************************************************************************
/*!
  \file      src/Control/MeshConv/Types.h
  \author    J. Bakosi
  \date      Mon 14 Jul 2014 09:32:48 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Types for MeshConv's parsers
  \details   Types for MeshConv's parsers
*/
//******************************************************************************
#ifndef MeshConvTypes_h
#define MeshConvTypes_h

#include <pegtl.hh>

#include <tkTags.h>
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
