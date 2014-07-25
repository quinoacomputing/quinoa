//******************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/Keywords.h
  \author    J. Bakosi
  \date      Thu 24 Jul 2014 11:26:11 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     UnitTest's command line keywords
  \details   All keywords recognized by UnitTest's command line parser. The
  keywords are defined by specializing struct 'keyword', defined in
  Control/Keyword.h. Introducing a new keyword requires a more human readable
  (but still short) name as well as a short, few-line, help-like description.
*/
//******************************************************************************
#ifndef UnitTestCmdLineKeywords_h
#define UnitTestCmdLineKeywords_h

//! Signal to compiler that we are building a list of keywords. This is used by
//! the inline includes, such as *Keywords.h, below (if any) to make sure they
//! get included in the correct namespace and not polluting the global one.
#define Keywords

#include <Keyword.h>
#include <CmdLineBaseKeywords.h>

namespace unittest {
//! List of keywords the parser understands
namespace kw {

using namespace pegtl::ascii;
using tk::kw::cmdline_keyword;

} // kw::
} // unittest::

#undef Keywords

#endif // UnitTestCmdLineKeywords_h
