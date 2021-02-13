// *****************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/Parser.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     UnitTest's command line parser
  \details   This file defines the command-line argument parser for the unit
     test suite, UnitTest.
*/
// *****************************************************************************

#include "NoWarning/pegtl.hpp"

#include "Print.hpp"
#include "UnitTest/Types.hpp"
#include "UnitTest/CmdLine/Parser.hpp"
#include "UnitTest/CmdLine/Grammar.hpp"
#include "UnitTest/CmdLine/CmdLine.hpp"

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
// *****************************************************************************
//  Contructor: parse the command line for UnitTest
//! \param[in] argc Number of C-style character arrays in argv
//! \param[in] argv C-style character array of character arrays
//! \param[in] print Pretty printer
//! \param[inout] cmdline Command-line stack where data is stored from parsing
//! \param[inout] helped Boolean indicating if command-line help was requested
// *****************************************************************************
{
  // Create CmdLine (a tagged tuple) to store parsed input
  ctr::CmdLine cmd;

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

  // Parse command line string by populating the underlying tagged tuple
  tao::pegtl::memory_input<> in( m_string, "command line" );
  tao::pegtl::parse< cmd::read_string, tk::grm::action >( in, cmd );

  // Echo errors and warnings accumulated during parsing
  diagnostics( print, cmd.get< tag::error >() );

  // Strip command line (and its underlying tagged tuple) from PEGTL instruments
  // and transfer it out
  cmdline = std::move( cmd );

  // If we got here, the parser has succeeded
  print.item("Parsed command line", "success");

  // Print out help on all command-line arguments if requested
  const auto helpcmd = cmdline.get< tag::help >();
  if (helpcmd) {
    print.help< tk::QUIET >( tk::unittest_executable(),
                             cmdline.get< tag::cmdinfo >(),
                             "Command-line Parameters:", "-" );
   print.mandatory< tk::QUIET >( "None of the arguments are mandatory." );
  }

  // Print out verbose help for a single keyword if requested
  const auto helpkw = cmdline.get< tag::helpkw >();
  if (!helpkw.keyword.empty())
    print.helpkw< tk::QUIET >( tk::unittest_executable(), helpkw );

  // Print out version information if it was requested
  const auto version = cmdline.get< tag::version >();
  if (version)
    print.version< tk::QUIET >( tk::unittest_executable(),
                                tk::quinoa_version(),
                                tk::git_commit(),
                                tk::copyright() );

  // Print out license information if it was requested
  const auto license = cmdline.get< tag::license >();
  if (license)
    print.license< tk::QUIET >( tk::unittest_executable(), tk::license() );

  // Will exit in main chare constructor if any help was output
  if (cmdline.get< tag::help >() ||           // help on all cmdline args
      !cmdline.get< tag::helpkw >().keyword.empty() || // help on a keyword
      version || license)                     // version or license output
    helped = true;
  else
    helped = false;
}
