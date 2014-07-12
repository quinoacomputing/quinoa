//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Parser.C
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:13:12 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's comamnd line parser
  \details   Quinoa's comamnd line parser
*/
//******************************************************************************

#include <make_unique.h>

#include <Quinoa/CmdLine/Parser.h>
#include <Quinoa/CmdLine/Grammar.h>

using quinoa::CmdLineParser;

CmdLineParser::CmdLineParser(int argc, char** argv,
                             const tk::Print& print,
                             std::unique_ptr< ctr::CmdLine >& cmdline) :
  StringParser( argc, argv )
//******************************************************************************
//  Contructor: parse command line
//! \author  J. Bakosi
//******************************************************************************
{
  // Create PEGTL string input from std::string (i.e. concatenated argv[])
  pegtl::string_input< ctr::Location > input( m_string );

  // Create std::unique_ptr behind which to store parsed command line data:
  // PEGTLCmdLine derives from CmdLine and has location() used during parsing
  std::unique_ptr< cmd::PEGTLCmdLine >
    pcl( tk::make_unique< cmd::PEGTLCmdLine >( input ) );

  // Parse command line string by populating the underlying tagged tuple:
  // basic_parse() below gives debug info during parsing, use it for debugging
  // the parser itself, i.e., when modifying the grammar, otherwise, use
  // dummy_parse() to compile faster
  pegtl::dummy_parse< cmd::read_string >( input, *pcl );

  // Strip cmdline (and its underlying tagged tuple) from PEGTL instruments by
  // creating a unique_ptr to the base class (CmdLine) and transfer it out
  cmdline = std::unique_ptr< ctr::CmdLine >( std::move(pcl) );

  // If we got here, parser succeeded
  print.item("Parsed command line", "success");

  // Make sure mandatory arguments are set
  ErrChk( !(cmdline->get<tag::io, tag::control>().empty()),
          "Mandatory control file not specified. "
          "Use '--" + kw::control().string() + " <filename>' or '-" +
          kw::control().alias() + " <filename>'.");
}
