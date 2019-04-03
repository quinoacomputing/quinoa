// *****************************************************************************
/*!
  \file      src/Control/UnitTest/Types.hpppp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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
#include "Keyword.hpp"

namespace unittest {
namespace ctr {

using namespace tao;

//! PEGTL location/position type to use throughout all of UnitTest's parsers
using Location = pegtl::position;

} // ctr::
} // unittest::

#endif // UnitTestTypes_h
