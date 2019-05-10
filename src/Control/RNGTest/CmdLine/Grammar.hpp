// *****************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/Grammar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     RNGTest's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef RNGTestCmdLineGrammar_h
#define RNGTestCmdLineGrammar_h

#include "CommonGrammar.hpp"
#include "Keywords.hpp"

namespace rngtest {
//! RNGTest command line grammar definition
namespace cmd {

  using namespace tao;

  //! \brief Specialization of tk::grm::use for RNGTest's command line parser
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords::set >;

  // RNGTest's CmdLine grammar

  //! \brief Match and set verbose switch (i.e., verbose or quiet output)
  struct verbose :
         tk::grm::process_cmd_switch< use, kw::verbose,
                                      tag::verbose > {};

  //! Match and set chare state switch
  struct charestate :
         tk::grm::process_cmd_switch< use, kw::charestate,
                                      tag::chare > {};

  //! \brief Match and set control (i.e., input deck) file name
  struct control :
         tk::grm::process_cmd< use, kw::control,
                               tk::grm::Store< tag::io, tag::control >,
                               pegtl::any,
                               tag::io, tag::control > {};

  //! \brief Match help on control file keywords
  struct helpctr :
         tk::grm::process_cmd_switch< use, kw::helpctr,
                                      tag::helpctr > {};

  //! \brief Match help on command-line parameters
  struct help :
         tk::grm::process_cmd_switch< use, kw::help,
                                      tag::help > {};

  //! \brief Match help on a command-line keyword
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

  //! Match all command line keywords
  struct keywords :
         pegtl::sor< verbose,
                     charestate,
                     control,
                     help,
                     helpctr,
                     helpkw,
                     quiescence,
                     trace > {};

  //! \brief Grammar entry point: parse keywords until end of string
  struct read_string :
         tk::grm::read_string< keywords > {};

} // cmd::
} // rngtest::

#endif // RNGTestCmdLineGrammar_h
