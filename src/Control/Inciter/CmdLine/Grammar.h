// *****************************************************************************
/*!
  \file      src/Control/Inciter/CmdLine/Grammar.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

  using namespace tao;

  //! Specialization of tk::grm::use for Inciter's command line parser
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords >;

  // Inciter's CmdLine grammar

  //! Match and set verbose switch (i.e., verbose or quiet output)
  struct verbose :
         tk::grm::process_cmd_switch< use< kw::verbose >,
                                      tag::verbose > {};

  //! Match and set chare state switch
  struct charestate :
         tk::grm::process_cmd_switch< use< kw::charestate >,
                                      tag::chare > {};

  //! Match and set non-blocking (migration) switch
  struct nonblocking :
         tk::grm::process_cmd_switch< use< kw::nonblocking >,
                                      tag::nonblocking > {};


  //! Match and set benchmark switch (i.e., benchmark mode)
  struct benchmark :
         tk::grm::process_cmd_switch< use< kw::benchmark >,
                                      tag::benchmark > {};

  //! Match and set feedback switch (i.e., feedback mode)
  struct feedback :
         tk::grm::process_cmd_switch< use< kw::feedback >,
                                      tag::feedback > {};

  //! Match and set virtualization parameter
  struct virtualization :
         tk::grm::process_cmd< use< kw::virtualization >,
                               tk::grm::Store< tag::virtualization >,
                               tk::grm::number > {};

  //! Match and set io parameter
  template< typename keyword, typename io_tag >
  struct io :
         tk::grm::process_cmd< keyword,
                               tk::grm::Store< tag::io, io_tag > > {};

  //! Match help on command-line parameters
  struct help :
         tk::grm::process_cmd_switch< use< kw::help >,
                                      tag::help > {};


  //! Match help on control file keywords
  struct helpctr :
         tk::grm::process_cmd_switch< use< kw::helpctr >,
                                      tag::helpctr > {};

  //! Match help on a command-line or control file keyword
  struct helpkw :
         tk::grm::process_cmd< use< kw::helpkw >,
                               tk::grm::helpkw,
                               pegtl::alnum > {};

  //! Match help on control file keywords
  struct quiescence :
         tk::grm::process_cmd_switch< use< kw::quiescence >,
                                      tag::quiescence > {};

  //! Match all command line keywords
  struct keywords :
         pegtl::sor< verbose,
                     charestate,
                     nonblocking,
                     benchmark,
                     feedback,
                     virtualization,
                     help,
                     helpctr,
                     helpkw,
                     quiescence,
                     io< use< kw::control >, tag::control >,
                     io< use< kw::input >, tag::input >,
                     io< use< kw::output >, tag::output >,
                     io< use< kw::diagnostics >, tag::diag > > {};

  //! Grammar entry point: parse keywords until end of string
  struct read_string :
         tk::grm::read_string< keywords > {};

} // cmd::
} // inciter::

#endif // InciterCmdLineGrammar_h
