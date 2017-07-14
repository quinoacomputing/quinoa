// *****************************************************************************
/*!
  \file      src/Control/MeshConv/Types.h
  \author    A. Pakki
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Types for FileConv's parsers
  \details   Types for UnitTest's parsers. This file defines the components of
    the tagged tuple that stores heteroegeneous objects in a hierarchical way.
    These components are therefore part of the grammar stack that is filled
    during command-line argument parsing.
*/
// *****************************************************************************
#ifndef FileConvTypes_h
#define FileConvTypes_h

#include "TaggedTuple.h"
#include "Tags.h"
#include "Keyword.h"

namespace fileconv {
namespace ctr {

using namespace tao;

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::input,           std::string,    //!< Input filename
  tag::output,          std::string     //!< Output filename
>;

//! PEGTL location/position type to use throughout all of MeshConv's parsers
using Location = pegtl::position;

} // ctr::
} // meshconv::

#endif // FileConvTypes_h
