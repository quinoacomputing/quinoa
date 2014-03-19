//******************************************************************************
/*!
  \file      src/Control/Gmsh2Exo/CmdLine/Keywords.h
  \author    J. Bakosi
  \date      Wed Mar 19 10:16:15 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh2Exo's command line keywords
  \details   All keywords recognized by Gmsh2Exo's command line parser. The
  keywords are defined by specializing struct 'keyword', defined in
  Control/Keyword.h. Introducing a new keyword requires a more human readable
  (but still short) name as well as a short, few-line, help-like description.
*/
//******************************************************************************
#ifndef Gmsh2ExoCmdLineKeywords_h
#define Gmsh2ExoCmdLineKeywords_h

//! Signal to compiler that we are building a list of keywords. This is used by
//! the inline includes, such as *Keywords.h, below (if any) to make sure they
//! get included in the correct namespace and not polluting the global one.
#define Keywords

#include <Keyword.h>

namespace gmsh2exo {
//! List of keywords the parser understands
namespace kw {

using namespace pegtl::ascii;
using tk::kw::cmdline_keyword;
using tk::kw::undefined_info;

// Keyword 'control', cmdline '--control' with alias '-c'
using control = cmdline_keyword<undefined_info, c, c,o,n,t,r,o,l>;

} // kw::
} // gmsh2exo::

#undef Keywords

#endif // Gmsh2ExoCmdLineKeywords_h
