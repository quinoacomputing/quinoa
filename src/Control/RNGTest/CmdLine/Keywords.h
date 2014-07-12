//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/Keywords.h
  \author    J. Bakosi
  \date      Sat 12 Jul 2014 06:44:25 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     RNGTest's command line keywords
  \details   All keywords recognized by RNGTest's command line parser. The
  keywords are defined by specializing struct 'keyword', defined in
  Control/Keyword.h. Introducing a new keyword requires a more human readable
  (but still short) name as well as a short, few-line, help-like description.
*/
//******************************************************************************
#ifndef RNGTestCmdLineKeywords_h
#define RNGTestCmdLineKeywords_h

//! Signal to compiler that we are building a list of keywords. This is used by
//! the inline includes, such as *Keywords.h, below (if any) to make sure they
//! get included in the correct namespace and not polluting the global one.
#define Keywords

#include <Keyword.h>
#include <CmdLineBaseKeywords.h>

namespace rngtest {
//! List of keywords the parser understands
namespace kw {

using namespace pegtl::ascii;
using tk::kw::cmdline_keyword;

// Keyword 'control', cmdline '--control' with alias '-c'
struct control_info {
  static const char* name() { return "control"; }
  static const char* help() { return
    "This option is used to define the control file containing user input.";
  }
};
using control = cmdline_keyword< control_info, c, c,o,n,t,r,o,l >;

} // kw::
} // rngtest::

#undef Keywords

#endif // RNGTestCmdLineKeywords_h
