//******************************************************************************
/*!
  \file      src/Control/UnitTest/Types.h
  \author    J. Bakosi
  \date      Fri 20 Mar 2015 11:34:07 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Types for UnitTest's parsers
  \details   Types for UnitTest's parsers. This file defines the components of
    the tagged tuple that stores heteroegeneous objects in a hierarchical way.
    These components are therefore part of the grammar stack that is filled
    during command-line argument parsing.
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
