// *****************************************************************************
/*!
  \file      src/Control/FileDiff/CmdLine/Grammar.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     FileDiff's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef FileDiffCmdLineGrammar_h
#define FileDiffCmdLineGrammar_h

#include "CommonGrammar.h"
#include "Keywords.h"

namespace filediff {
//! File converter command line grammar definition
namespace cmd {

  using namespace tao;

  //! \brief Specialization of tk::grm::use for FileDiff's command line parser
  //! \author J. Bakosi
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords >;

  // FileDiff's CmdLine state

  // FileDiff's CmdLine grammar

  //! brief Match and set verbose switch (i.e., verbose or quiet output)
  struct verbose :
         tk::grm::process_cmd_switch< use< kw::verbose >, tag::verbose > {};

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

  //! \brief Match all command line keywords
  struct keywords :
         pegtl::sor< verbose,
                     help,
                     helpkw,
                     io< use< kw::input >, tag::input >,
                     io< use< kw::output >, tag::output > > {};

  //! \brief Grammar entry point: parse keywords until end of string
  struct read_string :
         tk::grm::read_string< keywords > {};

} // cmd::
} // filediff::

#endif // FileDiffCmdLineGrammar_h
