//******************************************************************************
/*!
  \file      src/Control/UnitTest/Types.h
  \author    J. Bakosi
  \date      Sat 17 Jan 2015 06:48:40 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Types for UnitTest's parsers
  \details   Types for UnitTest's parsers. This file defines the components of
    the tagged tuple that stores heteroegeneous objects in a hierarchical way.
    These components are therefore part of the grammar stack that is filled
    during command0line argument parsing.
*/
//******************************************************************************
#ifndef UnitTestTypes_h
#define UnitTestTypes_h

#include <Tags.h>
#include <Keyword.h>

namespace unittest {
namespace ctr {

//! PEGTL location type to use throughout all of UnitTest's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // unittest::

#endif // UnitTestTypes_h
