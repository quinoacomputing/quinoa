//******************************************************************************
/*!
  \file      src/Control/UnitTest/Types.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 09:20:20 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Types for UnitTest's parsers
  \details   Types for UnitTest's parsers
*/
//******************************************************************************
#ifndef UnitTestTypes_h
#define UnitTestTypes_h

#include <Tags.h>
#include <Keyword.h>

namespace unittest {
//! control and parsing
namespace ctr {

//! PEGTL location type to use throughout all of UnitTest's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // unittest::

#endif // UnitTestTypes_h
