// *****************************************************************************
/*!
  \file      src/Control/MeshConv/Types.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Types for MeshConv's parsers
  \details   Types for UnitTest's parsers. This file defines the components of
    the tagged tuple that stores heteroegeneous objects in a hierarchical way.
    These components are therefore part of the grammar stack that is filled
    during command-line argument parsing.
*/
// *****************************************************************************
#ifndef MeshConvTypes_h
#define MeshConvTypes_h

#include "TaggedTuple.hpp"
#include "Tags.hpp"
#include "Keyword.hpp"

namespace meshconv {
namespace ctr {

using namespace tao;

//! IO parameters storage
using ios = tk::TaggedTuple< brigand::list<
    tag::nrestart,  int                             //!< Number of restarts
  , tag::input,     std::string                     //!< Input filename
  , tag::output,    std::string                     //!< Output filename
  , tag::screen,    kw::screen::info::expect::type  //!< Screen output filename
> >;

//! PEGTL location/position type to use throughout all of MeshConv's parsers
using Location = pegtl::position;

} // ctr::
} // meshconv::

#endif // MeshConvTypes_h
