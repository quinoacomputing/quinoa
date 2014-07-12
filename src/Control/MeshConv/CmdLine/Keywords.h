//******************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/Keywords.h
  \author    J. Bakosi
  \date      Sun 22 Jun 2014 10:36:40 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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
#include <BaseKeywords.h>

namespace meshconv {
//! List of keywords the parser understands
namespace kw {

using namespace pegtl::ascii;
using tk::kw::cmdline_keyword;

// Keyword 'input', cmdline '--input' with alias '-i'
struct input_info {
  static const char* name() { return "input"; }
  static const char* help() { return
    "This option is used to define the input file.";
  }
};
using input = cmdline_keyword< input_info, i, i,n,p,u,t >;

// Keyword 'output', cmdline '--output' with alias '-o'
struct output_info {
  static const char* name() { return "output"; }
  static const char* help() { return
    "This option is used to define the output file.";
  }
};
using output = cmdline_keyword< output_info, o, o,u,t,p,u,t >;

} // kw::
} // meshconv::

#undef Keywords

#endif // MeshConvCmdLineKeywords_h
