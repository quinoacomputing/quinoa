//******************************************************************************
/*!
  \file      src/Control/MeshConv/Types.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 09:19:09 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Types for MeshConv's parsers
  \details   Types for MeshConv's parsers
*/
//******************************************************************************
#ifndef MeshConvTypes_h
#define MeshConvTypes_h

#include <Tags.h>
#include <Keyword.h>

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
