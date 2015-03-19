//******************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/Parser.C
  \author    J. Bakosi
  \date      Thu 19 Mar 2015 11:49:11 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     UnitTest's comamnd line parser
  \details   This file defines the command-line argument parser for the unit
     test suite, UnitTest.
*/
//******************************************************************************
#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <charm++.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <Config.h>
#include <UnitTest/CmdLine/Parser.h>
#include <UnitTest/CmdLine/Grammar.h>

namespace tk {
namespace grm {

tk::Print g_print;

} // grm::
} // tk::

using unittest::CmdLineParser;

CmdLineParser::CmdLineParser( int argc,
                              char** argv,
                              const tk::Print& print,
                              ctr::CmdLine& cmdline,
                              bool& helped ) :
  StringParser( argc, argv )
//******************************************************************************
//  Contructor: parse the command line for UnitTest
//! \param[in] argc Number of C-style character arrays in argv
//! \param[in] argv C-style character array of character arrays
//! \param[in] print Pretty printer
//! \param[inout] helped Boolean indicating if command-line help was requested
//! \param[inout] cmdline Command-line stack where data is stored from parsing
//! \author  J. Bakosi
//******************************************************************************
{
  // Create PEGTL string input from std::string (i.e. concatenated argv[])
  pegtl::string_input< ctr::Location > input( m_string );

  // Create PEGTLCmdLine to store parsed command line data which derives from
  // CmdLine and has location() used during parsing
  cmd::PEGTLCmdLine cmd( input );

  // Reset parser's output stream to that of print's. This is so that mild
  // warnings emitted during parsing can be output using the pretty printer.
  // Usually, errors and warnings are simply accumulated during parsing and
  // printed during diagnostics after the parser has finished. However, in some
  // special cases we can provide a more user-friendly message right during
  // parsing since there is more information available to construct a more
  // sensible message. This is done in e.g., tk::grm::store_option. Resetting
  // the global g_print, to that of passed in as the constructor argument allows
  // not to have to create a new pretty printer, but use the existing one.
  tk::grm::g_print.reset( print.save() );

  // Parse command line string by populating the underlying tagged tuple:
  // basic_parse() below gives debug info during parsing, use it for debugging
  // the parser itself, i.e., when modifying the grammar, otherwise, use
  // dummy_parse() to compile faster
  pegtl::dummy_parse< cmd::read_string >( input, cmd );

  // Echo errors and warnings accumulated during parsing
  diagnostics( print, cmd.get< tag::error >() );

  // Strip command line (and its underlying tagged tuple) from PEGTL instruments
  // and transfer it out
  cmdline = std::move( cmd );

  // If we got here, the parser has succeeded
  print.item("Parsed command line", "success");

  // Will exit in main chare constructor if any help was output
  if (cmdline.get< tag::help >() ||           // help on all cmdline args
      !cmdline.get< tag::helpkw >().keyword.empty()) // help on a keyword
    helped = true;
  else
    helped = false;
}
