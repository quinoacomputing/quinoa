//******************************************************************************
/*!
  \file      src/Control/RegTest/Types.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:15:06 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Types for RegTest's parsers
  \details   Types for RegTest's parsers. This file defines the components of
    the tagged tuple that stores heteroegeneous objects in a hierarchical way.
    These components are therefore part of the grammar stack that is filled
    during command-line argument parsing.
*/
//******************************************************************************
#ifndef RegTestTypes_h
#define RegTestTypes_h

#include "Tags.h"
#include "Keyword.h"

namespace regtest {
namespace ctr {

//! PEGTL location type to use throughout all of RegTest's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // regtest::

#endif // RegTestTypes_h
