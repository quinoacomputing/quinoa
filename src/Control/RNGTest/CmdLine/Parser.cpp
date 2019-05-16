// *****************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/Parser.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     RNGTest's command line parser
  \details   This file defines the command-line argument parser for the random
    number generator test suite, RNGTest.
*/
// *****************************************************************************

#include <map>
#include <ostream>
#include <string>
#include <type_traits>

#include "NoWarning/pegtl.hpp"

#include "NoWarning/charm.hpp"
#include "QuinoaConfig.hpp"
#include "Exception.hpp"
#include "Print.hpp"
#include "Keywords.hpp"
#include "HelpFactory.hpp"
#include "RNGTest/Types.hpp"
#include "RNGTest/CmdLine/Parser.hpp"
#include "RNGTest/CmdLine/Grammar.hpp"
#include "RNGTest/CmdLine/CmdLine.hpp"
#include "RNGTest/InputDeck/InputDeck.hpp"

namespace tk {
namespace grm {

tk::Print g_print;

} // grm::
} // tk::

namespace rngtest {

extern ctr::InputDeck g_inputdeck;

} // rngtest::

using rngtest::CmdLineParser;

CmdLineParser::CmdLineParser( int argc,
                              char** argv,
                              const tk::Print& print,
                              ctr::CmdLine& cmdline ) :
  StringParser( argc, argv )
// *****************************************************************************
//  Contructor: parse the command line for RNGTest
//! \param[in] argc Number of C-style character arrays in argv
//! \param[in] argv C-style character array of character arrays
//! \param[in] print Pretty printer
//! \param[in,out] cmdline Command-line stack where data is stored from parsing
// *****************************************************************************
{
  // Create CmdLine (a tagged tuple) to store parsed input
  ctr::CmdLine cmd( g_inputdeck.get< tag::cmd, tag::ctrinfo >() );

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
  tao::pegtl::memory_input<> in( m_string, "command line" );
  tao::pegtl::parse< cmd::read_string, tk::grm::action >( in, cmd );

  // Echo errors and warnings accumulated during parsing
  diagnostics( print, cmd.get< tag::error >() );

  // Strip command line (and its underlying tagged tuple) from PEGTL instruments
  // and transfer it out
  cmdline = std::move( cmd );

  // If we got here, parser succeeded
  print.item( "Parsed command line", "success" );

  // Print out help on all command-line arguments if the executable was invoked
  // without arguments or the help was requested
  const auto helpcmd = cmdline.get< tag::help >();
  if (argc == 1 || helpcmd)
    print.help< tk::QUIET >( tk::rngtest_executable(),
                             cmdline.get< tag::cmdinfo >(),
                             "Command-line Parameters:", "-" );

  // Print out help on all control file keywords if they were requested
  const auto helpctr = cmdline.get< tag::helpctr >();
  if (helpctr)
    print.help< tk::QUIET >( tk::rngtest_executable(),
                             cmdline.get< tag::ctrinfo >(),
                             "Control File Keywords:" );

  // Print out verbose help for a single keyword if requested
  const auto helpkw = cmdline.get< tag::helpkw >();
  if (!helpkw.keyword.empty())
    print.helpkw< tk::QUIET >( tk::rngtest_executable(), helpkw );

  // Immediately exit if any help was output or was called without any argument
  // with zero exit code
  if (argc == 1 || helpcmd || helpctr || !helpkw.keyword.empty()) CkExit();

  // Make sure mandatory arguments are set
  auto alias = kw::control().alias();
  ErrChk( !(cmdline.get< tag::io, tag::control >().empty()),
          "Mandatory control file not specified. "
          "Use '--" + kw::control().string() + " <filename>'" +
          ( alias ? " or '-" + *alias + " <filename>'" : "" ) + '.' );
}