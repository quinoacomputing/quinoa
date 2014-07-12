//******************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/Parser.C
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:13:05 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     MeshConv's comamnd line parser
  \details   MeshConv's comamnd line parser
*/
//******************************************************************************

#include <make_unique.h>

#include <MeshConv/CmdLine/Parser.h>
#include <MeshConv/CmdLine/Grammar.h>

namespace tk {
namespace grm {

tk::Print g_print;

} // grm::
} // tk::

using meshconv::CmdLineParser;

CmdLineParser::CmdLineParser( int argc, char** argv,
                              const tk::Print& print,
                              ctr::CmdLine& cmdline ) :
  StringParser( argc, argv )
//******************************************************************************
//  Contructor: parse command line
//! \author  J. Bakosi
//******************************************************************************
{
  // Create PEGTL string input from std::string (i.e. concatenated argv[])
  pegtl::string_input< ctr::Location > input( m_string );

  // Create PEGTLCmdLine to store parsed command line data which derives from
  // CmdLine and has location() used during parsing
  cmd::PEGTLCmdLine cmd( input );

  // Reset parser's output stream to that of print's
  tk::grm::g_print.reset( print.save() );

  // Parse command line string by populating the underlying tagged tuple:
  // basic_parse() below gives debug info during parsing, use it for debugging
  // the parser itself, i.e., when modifying the grammar, otherwise, use
  // dummy_parse() to compile faster
  pegtl::dummy_parse< cmd::read_string >( input, cmd );

  // Strip command line (and its underlying tagged tuple) from PEGTL instruments
  // and transfer it out
  cmdline = cmd;

  // If we got here, parser succeeded
  print.item("Parsed command line", "success");

  // Make sure mandatory arguments are set
  ErrChk( !(cmdline.get< tag::io, tag::input >().empty()),
          "Mandatory input file not specified. "
          "Use '--" + kw::input().string() + " <filename>' or '-" +
          kw::input().alias() + " <filename>'.");
  ErrChk( !(cmdline.get< tag::io, tag::output >().empty()),
          "Mandatory output file not specified. "
          "Use '--" + kw::output().string() + " <filename>' or '-" +
          kw::output().alias() + " <filename>'.");
}
