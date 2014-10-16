//******************************************************************************
/*!
  \file      src/Control/CmdLineBaseKeywords.h
  \author    J. Bakosi
  \date      Sat 12 Jul 2014 11:10:11 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Basic keywords recognized by all command line parsers
  \details   Basic keywords recognized by all command line parsers
*/
//******************************************************************************
#ifndef Keywords
#error "CmdLineBaseKeywords.h should only be included within a *Keywords.h"
#endif

namespace tk {
namespace kw {

using namespace pegtl::ascii;

// Keyword 'verbose', cmdline '--verbose' with alias '-v'
struct verbose_info {
  static const char* name() { return "verbose"; }
  static const char* help() { return
    "This option is used to pick verbose output.";
  }
};
using verbose = cmdline_keyword< verbose_info, v, v,e,r,b,o,s,e >;

} // kw::
} // tk::
