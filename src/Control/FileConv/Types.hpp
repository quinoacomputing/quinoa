// *****************************************************************************
/*!
  \file      src/Control/MeshConv/Types.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Types for FileConv's parsers
  \details   Types for UnitTest's parsers. This file defines the components of
    the tagged tuple that stores heteroegeneous objects in a hierarchical way.
    These components are therefore part of the grammar stack that is filled
    during command-line argument parsing.
*/
// *****************************************************************************
#ifndef FileConvTypes_h
#define FileConvTypes_h

#include "TaggedTuple.hpp"
#include "Tags.hpp"
#include "Keyword.hpp"

namespace fileconv {
namespace ctr {

using namespace tao;

//! IO parameters storage
using ios = tk::TaggedTuple< brigand::list<
    tag::input,           std::string     //!< Input filename
  , tag::output,          std::string     //!< Output filename
> >;

//! PEGTL location/position type to use throughout all of MeshConv's parsers
using Location = pegtl::position;

} // ctr::
} // meshconv::

#endif // FileConvTypes_h
