//******************************************************************************
/*!
  \file      src/Control/MeshConv/Types.h
  \author    J. Bakosi
  \date      Tue 24 Nov 2015 08:18:52 AM MST
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Types for MeshConv's parsers
  \details   Types for UnitTest's parsers. This file defines the components of
    the tagged tuple that stores heteroegeneous objects in a hierarchical way.
    These components are therefore part of the grammar stack that is filled
    during command-line argument parsing.
*/
//******************************************************************************
#ifndef MeshConvTypes_h
#define MeshConvTypes_h

#include "TaggedTuple.h"
#include "Tags.h"
#include "Keyword.h"

namespace meshconv {
namespace ctr {

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::input,           std::string,    //!< Input filename
  tag::output,          std::string     //!< Output filename
>;

//! PEGTL location type to use throughout all of MeshConv's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // meshconv::

#endif // MeshConvTypes_h
