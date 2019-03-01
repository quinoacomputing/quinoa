// *****************************************************************************
/*!
  \file      src/Control/FileConv/CmdLine/Grammar.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     FileConv's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef FileConvCmdLineGrammar_h
#define FileConvCmdLineGrammar_h

#include "CommonGrammar.h"
#include "Keywords.h"

namespace fileconv {
//! File converter command line grammar definition
namespace cmd {

  using namespace tao;

  //! \brief Specialization of tk::grm::use for FileConv's command line parser
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords >;

  // FileConv's CmdLine state

  // FileConv's CmdLine grammar

  //! brief Match and set verbose switch (i.e., verbose or quiet output)
  struct verbose :
         tk::grm::process_cmd_switch< use< kw::verbose >, tag::verbose > {};


  //! Match and set chare state switch
  struct charestate :
         tk::grm::process_cmd_switch< use< kw::charestate >,
                                      tag::chare > {};

  //! \brief Match and set io parameter
  template< typename keyword, typename io_tag >
  struct io :
         tk::grm::process_cmd< keyword, tk::grm::Store< tag::io, io_tag > > {};

  //! \brief Match help on command-line parameters
  struct help :
         tk::grm::process_cmd_switch< use< kw::help >, tag::help > {};

  //! \brief Match help on a single command-line or control file keyword
  struct helpkw :
         tk::grm::process_cmd< use< kw::helpkw >,
                               tk::grm::helpkw,
                               pegtl::alnum > {};

  //! Match help on control file keywords
  struct quiescence :
         tk::grm::process_cmd_switch< use< kw::quiescence >,
                                      tag::quiescence > {};

  //! \brief Match all command line keywords
  struct keywords :
         pegtl::sor< verbose,
                     charestate,
                     help,
                     helpkw,
                     quiescence,
                     io< use< kw::input >, tag::input >,
                     io< use< kw::output >, tag::output > > {};

  //! \brief Grammar entry point: parse keywords until end of string
  struct read_string :
         tk::grm::read_string< keywords > {};

} // cmd::
} // fileconv::

#endif // FileConvCmdLineGrammar_h
