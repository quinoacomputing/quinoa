// *****************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/Grammar.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     UnitTest's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef UnitTestCmdLineGrammar_h
#define UnitTestCmdLineGrammar_h

#include "CommonGrammar.h"
#include "Keywords.h"

namespace unittest {
//! UnitTest command line grammar definition
namespace cmd {

  //! \brief Specialization of tk::grm::use for UnitTest's command line parser
  //! \author J. Bakosi
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords >;

  // UnitTest's CmdLine state

  // UnitTest's CmdLine grammar

  //! \brief Match and set verbose switch (i.e., verbose or quiet output)
  //! \author J. Bakosi
  struct verbose :
         tk::grm::process_cmd_switch< use< kw::verbose >,
                                      tag::verbose > {};

  //! \brief Match help on command-line parameters
  //! \author J. Bakosi
  struct help :
         tk::grm::process_cmd_switch< use< kw::help >,
                                      tag::help > {};

  //! \brief Match help on a command-line keyword
  //! \author J. Bakosi
  struct helpkw :
         tk::grm::process_cmd< use< kw::helpkw >,
                               tk::grm::helpkw,
                               pegtl::alnum > {};

  //! \brief Match test group name(s) and only run those
  //! \author J. Bakosi
  struct group :
         tk::grm::process_cmd< use< kw::group >,
                               tk::grm::Store< tag::group > > {};

  //! \brief Match all command line keywords
  //! \author J. Bakosi
  struct keywords :
         pegtl::sor< verbose, help, helpkw, group > {};

  //! \brief Grammar entry point: parse keywords until end of string
  //! \author J. Bakosi
  struct read_string :
         tk::grm::read_string< keywords > {};

} // cmd::
} // unittest::

#endif // UnitTestCmdLineGrammar_h
