//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Fri 16 Jan 2015 06:14:18 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Quinoa's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef QuinoaCmdLineGrammar_h
#define QuinoaCmdLineGrammar_h

#include <Grammar.h>
#include <PEGTLParsed.h>
#include <Keywords.h>

namespace quinoa {
//! Quinoa command line grammar definition
namespace cmd {

  //! PEGTLParsed type specialized to Quinoa's command line parser
  //! \details PEGTLCmdLine is practically CmdLine equipped with PEGTL location
  //!    information so the location can be tracked during parsing.
  //! \author J. Bakosi
  using PEGTLCmdLine =
    tk::ctr::PEGTLParsed< ctr::CmdLine, pegtl::string_input< ctr::Location > >;

  //! \brief Specialization of tk::grm::use for Quinoa's command line parser
  //! \author J. Bakosi
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords >;

  // Quinoa's CmdLine state

  //! Everything is stored in Stack during parsing
  //! \author J. Bakosi
  using Stack = PEGTLCmdLine;

  // Quinoa's CmdLine grammar

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
                     io< use< kw::output >, tag::output >,
                     io< use< kw::pdf >, tag::pdf >,
                     io< use< kw::glob >, tag::glob >,
                     io< use< kw::stat >, tag::stat > > {};

  //! \brief Grammar entry point: parse keywords until end of string
  //! \author J. Bakosi
  struct read_string :
         tk::grm::read_string< Stack, keywords > {};

} // cmd::
} // quinoa::

#endif // QuinoaCmdLineGrammar_h
