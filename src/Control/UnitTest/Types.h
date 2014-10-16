//******************************************************************************
/*!
  \file      src/Control/UnitTest/Types.h
  \author    J. Bakosi
  \date      Thu 24 Jul 2014 11:25:18 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Types for UnitTest's parsers
  \details   Types for UnitTest's parsers
*/
//******************************************************************************
#ifndef UnitTestTypes_h
#define UnitTestTypes_h

#include <pegtl.hh>

#include <tkTags.h>
#include <UnitTest/Tags.h>

namespace unittest {
//! control and parsing
namespace ctr {

//! PEGTL location type to use throughout all of UnitTest's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // unittest::

#endif // UnitTestTypes_h
