//******************************************************************************
/*!
  \file      src/Control/UnitTest/Types.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 07:19:22 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::input,       std::string,  //!< Input filename
  tag::output,      std::string   //!< Output filename
>;

//! PEGTL location type to use throughout all of UnitTest's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // unittest::

#endif // UnitTestTypes_h
