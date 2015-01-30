//******************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/Parser.C
  \author    J. Bakosi
  \date      Sun 18 Jan 2015 07:10:19 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     UnitTest's comamnd line parser
  \details   This file defines the command-line argument parser for the unit
     test suite, UnitTest.
*/
//******************************************************************************

#include <charm++.h>

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
                              ctr::CmdLine& cmdline ) :
  StringParser( argc, argv )
//******************************************************************************
//  Contructor: parse the command line for UnitTest
//! \param[in] argc Number of C-style character arrays in argv
//! \param[in] argv C-style character array of character arrays
//! \param[in] print Pretty printer
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

  // Print out help on all command-line arguments if the help was requested
  const auto helpcmd = cmdline.get< tag::help >();
  if (helpcmd)
    print.help< tk::QUIET >( UNITTEST_EXECUTABLE, cmdline.get< tag::cmdinfo >(),
                             "Command-line Parameters:", "-" );

  // Print out verbose help for a single keyword if requested
  const auto helpkw = cmdline.get< tag::helpkw >();
  if (!helpkw.keyword.empty())
    print.helpkw< tk::QUIET >( UNITTEST_EXECUTABLE, helpkw );

  // Immediately exit if any help was output
  if (helpcmd || !helpkw.keyword.empty()) CkExit();
}
