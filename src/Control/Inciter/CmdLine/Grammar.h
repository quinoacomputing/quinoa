// *****************************************************************************
/*!
  \file      src/Control/Inciter/CmdLine/Grammar.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Inciter's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef InciterCmdLineGrammar_h
#define InciterCmdLineGrammar_h

#include "CommonGrammar.h"
#include "Keywords.h"

namespace inciter {
//! Inciter command line grammar definition
namespace cmd {

  //! \brief Specialization of tk::grm::use for Inciter's command line parser
  //! \author J. Bakosi
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords >;

  // Inciter's CmdLine grammar

  //! \brief Match and set verbose switch (i.e., verbose or quiet output)
  //! \author J. Bakosi
  struct verbose :
         tk::grm::process_cmd_switch< use< kw::verbose >,
                                      tag::verbose > {};

  //! \brief Match and set benchmark switch (i.e., benchmark mode)
  //! \author J. Bakosi
  struct benchmark :
         tk::grm::process_cmd_switch< use< kw::benchmark >,
                                      tag::benchmark > {};

  //! \brief Match and set feedback switch (i.e., feedback mode)
  //! \author J. Bakosi
  struct feedback :
         tk::grm::process_cmd_switch< use< kw::feedback >,
                                      tag::feedback > {};

  //! \brief Match and set virtualization parameter
  //! \author J. Bakosi
  struct virtualization :
         tk::grm::process_cmd< use< kw::virtualization >,
                               tk::grm::Store< tag::virtualization >,
                               tk::grm::number > {};

  //! \brief Match and set io parameter
  //! \author J. Bakosi
  template< typename keyword, typename io_tag >
  struct io :
         tk::grm::process_cmd< keyword,
                               tk::grm::Store< tag::io, io_tag > > {};

  //! \brief Match help on command-line parameters
  //! \author J. Bakosi
  struct help :
         tk::grm::process_cmd_switch< use< kw::help >,
                                      tag::help > {};


  //! \brief Match help on control file keywords
  //! \author J. Bakosi
  struct helpctr :
         tk::grm::process_cmd_switch< use< kw::helpctr >,
                                      tag::helpctr > {};

  //! \brief Match help on a command-line or control file keyword
  //! \author J. Bakosi
  struct helpkw :
         tk::grm::process_cmd< use< kw::helpkw >,
                               tk::grm::helpkw,
                               pegtl::alnum > {};

  //! \brief Match all command line keywords
  //! \author J. Bakosi
  struct keywords :
         pegtl::sor< verbose,
                     benchmark,
                     feedback,
                     virtualization,
                     help,
                     helpctr,
                     helpkw,
                     io< use< kw::control >, tag::control >,
                     io< use< kw::input >, tag::input >,
                     io< use< kw::output >, tag::output >,
                     io< use< kw::diagnostics >, tag::diag > > {};

  //! \brief Grammar entry point: parse keywords until end of string
  //! \author J. Bakosi
  struct read_string :
         tk::grm::read_string< keywords > {};

} // cmd::
} // inciter::

#endif // InciterCmdLineGrammar_h
