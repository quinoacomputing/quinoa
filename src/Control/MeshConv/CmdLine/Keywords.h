//******************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/Keywords.h
  \author    J. Bakosi
  \date      Tue 08 Apr 2014 09:26:49 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshConv's command line keywords
  \details   All keywords recognized by MeshConv's command line parser. The
  keywords are defined by specializing struct 'keyword', defined in
  Control/Keyword.h. Introducing a new keyword requires a more human readable
  (but still short) name as well as a short, few-line, help-like description.
*/
//******************************************************************************
#ifndef MeshConvCmdLineKeywords_h
#define MeshConvCmdLineKeywords_h

//! Signal to compiler that we are building a list of keywords. This is used by
//! the inline includes, such as *Keywords.h, below (if any) to make sure they
//! get included in the correct namespace and not polluting the global one.
#define Keywords

#include <Keyword.h>

namespace meshconv {
//! List of keywords the parser understands
namespace kw {

using namespace pegtl::ascii;
using tk::kw::cmdline_keyword;
using tk::kw::undefined_info;

// Keyword 'input', cmdline '--input' with alias '-i'
using input = cmdline_keyword<undefined_info, i, i,n,p,u,t>;

// Keyword 'output', cmdline '--output' with alias '-o'
using output = cmdline_keyword<undefined_info, o, o,u,t,p,u,t>;

} // kw::
} // meshconv::

#undef Keywords

#endif // MeshConvCmdLineKeywords_h
