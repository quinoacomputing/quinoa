//******************************************************************************
/*!
  \file      src/Control/Inciter/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 08:53:34 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef InciterCmdLineGrammar_h
#define InciterCmdLineGrammar_h

#include <Grammar.h>
#include <PEGTLParsed.h>
#include <Keywords.h>

namespace inciter {
//! Inciter command line grammar definition
namespace cmd {

  //! PEGTLParsed type specialized to Inciter's command line parser
  //! \details PEGTLCmdLine is practically CmdLine equipped with PEGTL location
  //!    information so the location can be tracked during parsing.
  //! \author J. Bakosi
  using PEGTLCmdLine =
    tk::ctr::PEGTLParsed< ctr::CmdLine, pegtl::string_input< ctr::Location > >;

  //! \brief Specialization of tk::grm::use for Inciter's command line parser
  //! \author J. Bakosi
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords >;

  // Inciter's CmdLine state

  //! Everything is stored in Stack during parsing
  //! \author J. Bakosi
  using Stack = PEGTLCmdLine;

  // Inciter's CmdLine grammar

  //! \brief Match and set verbose switch (i.e., verbose or quiet output)
  //! \author J. Bakosi
  struct verbose :
         tk::grm::process_cmd_switch< Stack,
                                      use< kw::verbose >,
                                      tag::verbose > {};

  //! \brief Match and set virtualization parameter
  //! \author J. Bakosi
  struct virtualization :
         tk::grm::process_cmd< Stack,
                               use< kw::virtualization >,
                               tk::grm::Store< Stack, tag::virtualization >,
                               tk::grm::number > {};

  //! \brief Match and set io parameter
  //! \author J. Bakosi
  template< typename keyword, typename io_tag >
  struct io :
         tk::grm::process_cmd< Stack,
                               keyword,
                               tk::grm::Store< Stack, tag::io, io_tag > > {};

  //! \brief Match help on command-line parameters
  //! \author J. Bakosi
  struct help :
         tk::grm::process_cmd_switch< Stack,
                                      use< kw::help >,
                                      tag::help > {};


  //! \brief Match help on control file keywords
  //! \author J. Bakosi
  struct helpctr :
         tk::grm::process_cmd_switch< Stack,
                                      use< kw::helpctr >,
                                      tag::helpctr > {};

  //! \brief Match help on a command-line or control file keyword
  //! \author J. Bakosi
  struct helpkw :
         tk::grm::process_cmd< Stack,
                               use< kw::helpkw >,
                               tk::grm::helpkw< Stack >,
                               pegtl::alnum > {};

  //! \brief Match all command line keywords
  //! \author J. Bakosi
  struct keywords :
         pegtl::sor< verbose,
                     virtualization,
                     help,
                     helpctr,
                     helpkw,
                     io< use< kw::control >, tag::control >,
                     io< use< kw::input >, tag::input >,
                     io< use< kw::output >, tag::output > > {};

  //! \brief Grammar entry point: parse keywords until end of string
  //! \author J. Bakosi
  struct read_string :
         tk::grm::read_string< Stack, keywords > {};

} // cmd::
} // inciter::

#endif // InciterCmdLineGrammar_h
