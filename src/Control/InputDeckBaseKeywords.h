//******************************************************************************
/*!
  \file      src/Control/InputDeckBaseKeywords.h
  \author    J. Bakosi
  \date      Sat 12 Jul 2014 06:39:24 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Basic keywords recognized by all input deck parsers
  \details   Basic keywords recognized by all input deck parsers
*/
//******************************************************************************
#ifndef Keywords
#error "InputDeckBaseKeywords.h should only be included within a *Keywords.h"
#endif

namespace tk {
namespace kw {

using namespace pegtl::ascii;

// Keyword 'title'
struct title_info {
  static const char* name() { return "Analysis title"; }
  static const char* help() { return
    "The analysis title may be specified in the input file using the 'title' "
    "keyword. The 'title' keyword must be followed by a doubly-quoted string "
    "specifying the analysis title. Example: title \"Example problem\".";
  }
};
using title = keyword< title_info, t,i,t,l,e >;

// Keyword 'end'
struct end_info {
  static const char* name() { return "End of block"; }
  static const char* help() { return
    "The end of a block is given by the 'end' keyword.";
  }
};
using end = keyword< end_info, e,n,d >;

} // kw::
} // tk::
