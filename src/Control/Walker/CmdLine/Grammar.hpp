// *****************************************************************************
/*!
  \file      src/Control/Walker/CmdLine/Grammar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Walker's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef WalkerCmdLineGrammar_h
#define WalkerCmdLineGrammar_h

#include "CommonGrammar.hpp"
#include "Keywords.hpp"

namespace walker {
//! Walker command line grammar definition
namespace cmd {

  using namespace tao;

  //! \brief Specialization of tk::grm::use for Walker's command line parser
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords::set >;

  // Walker's CmdLine grammar

  //! verbose (i.e., verbose or quiet output)
  struct verbose :
         tk::grm::process_cmd_switch< use, kw::verbose,
                                      tag::verbose > {};

  //! Match and set chare state switch
  struct charestate :
         tk::grm::process_cmd_switch< use, kw::charestate,
                                      tag::chare > {};

  //! virtualization parameter
  struct virtualization :
         tk::grm::process_cmd< use, kw::virtualization,
                               tk::grm::Store< tag::virtualization >,
                               tk::grm::number,
                               tag::virtualization > {};

  //! io parameter
  template< typename keyword, typename io_tag >
  struct io :
         tk::grm::process_cmd< use, keyword,
                               tk::grm::Store< tag::io, io_tag >,
                               pegtl::any,
                               tag::io, io_tag > {};

  //! help on control file keywords
  struct helpctr :
         tk::grm::process_cmd_switch< use, kw::helpctr,
                                      tag::helpctr > {};

  //! help on command-line parameters
  struct help :
         tk::grm::process_cmd_switch< use, kw::help, tag::help > {};

  //! help on a command-line keyword
  struct helpkw :
         tk::grm::process_cmd< use, kw::helpkw,
                               tk::grm::helpkw,
                               pegtl::alnum,
                               tag::discr /* = unused */ > {};

  //! Match switch on quiescence
  struct quiescence :
         tk::grm::process_cmd_switch< use, kw::quiescence,
                                      tag::quiescence > {};

  //! Match switch on trace output
  struct trace :
         tk::grm::process_cmd_switch< use, kw::trace,
                                      tag::trace > {};

  //! Match switch on version output
  struct version :
         tk::grm::process_cmd_switch< use, kw::version,
                                      tag::version > {};

  //! Match switch on license output
  struct license :
         tk::grm::process_cmd_switch< use, kw::license,
                                      tag::license > {};

  //! command line keywords
  struct keywords :
         pegtl::sor< verbose,
                     charestate,
                     help,
                     helpctr,
                     helpkw,
                     virtualization,
                     quiescence,
                     trace,
                     version,
                     license,
                     io< kw::control, tag::control >,
                     io< kw::pdf, tag::pdf >,
                     io< kw::stat, tag::stat >,
                     io< kw::screen, tag::screen >,
                     io< kw::particles, tag::particles > > {};

  //! entry point: parse keywords and until end of string
  struct read_string :
         tk::grm::read_string< keywords > {};

} // cmd::
} // walker::

#endif // WalkerCmdLineGrammar_h
