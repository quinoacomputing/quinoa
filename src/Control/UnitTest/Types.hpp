// *****************************************************************************
/*!
  \file      src/Control/UnitTest/Types.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Types for UnitTest's parsers
  \details   Types for UnitTest's parsers. This file defines the components of
    the tagged tuple that stores heteroegeneous objects in a hierarchical way.
    These components are therefore part of the grammar stack that is filled
    during command-line argument parsing.
*/
// *****************************************************************************
#ifndef UnitTestTypes_h
#define UnitTestTypes_h

#include "Tags.hpp"
#include "Keywords.hpp"

namespace unittest {
namespace ctr {

using namespace tao;

//! PEGTL location/position type to use throughout all of UnitTest's parsers
using Location = pegtl::position;

//! IO parameters storage
using ios = tk::TaggedTuple< brigand::list<
    tag::nrestart,  int                             //!< Number of restarts
  , tag::screen,    kw::screen::info::expect::type  //!< Screen output filename
> >;


} // ctr::
} // unittest::

#endif // UnitTestTypes_h
