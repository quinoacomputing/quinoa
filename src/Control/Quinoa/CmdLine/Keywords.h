//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Keywords.h
  \author    J. Bakosi
  \date      Sun 10 Nov 2013 10:57:20 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's command line keywords
  \details   All keywords recognized by Quinoa's command line parser. The
  keywords are defined by specializing struct 'keyword', defined in
  Control/Keyword.h. Introducing a new keyword requires a more human readable
  (but still short) name as well as a short, few-line, help-like description.
*/
//******************************************************************************
#ifndef QuinoaCmdLineKeywords_h
#define QuinoaCmdLineKeywords_h

//! Signal to compiler that we are building a list of keywords. This is used by
//! the inline includes, such as *Keywords.h, below (if any) to make sure they
//! get included in the correct namespace and not polluting the global one.
#define Keywords

#include <Keyword.h>

namespace quinoa {
//! List of keywords the parser understands
namespace kw {

using namespace pegtl::ascii;
using tk::kw::cmdline_keyword;
using tk::kw::undefined_info;

// Keyword 'control', cmdline '--control' with alias '-c'
using control = cmdline_keyword<undefined_info, c, c,o,n,t,r,o,l>;

// Keyword 'input', cmdline '--input' with alias '-i'
using input = cmdline_keyword<undefined_info, i, i,n,p,u,t>;

// Keyword 'output', cmdline '--output' with alias '-o'
using output = cmdline_keyword<undefined_info, o, o,u,t,p,u,t>;

// Keyword 'pdf', cmdline '--pdf' with alias '-p'
using pdf = cmdline_keyword<undefined_info, p, p,d,f >;

// Keyword 'glob', cmdline '--glob' with alias '-g'
using glob = cmdline_keyword<undefined_info, g, g,l,o,b >;

// Keyword 'stat', cmdline '--stat' with alias '-s'
using stat = cmdline_keyword<undefined_info, s, s,t,a,t >;

} // kw::
} // quinoa::

#undef Keywords

#endif // QuinoaCmdLineKeywords_h
