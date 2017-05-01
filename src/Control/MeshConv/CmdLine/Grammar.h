// *****************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/Grammar.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     MeshConv's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef MeshConvCmdLineGrammar_h
#define MeshConvCmdLineGrammar_h

#include "CommonGrammar.h"
#include "Keywords.h"

namespace meshconv {
//! Mesh converter command line grammar definition
namespace cmd {

  //! \brief Specialization of tk::grm::use for MeshConv's command line parser
  //! \author J. Bakosi
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords >;

  // MeshConv's CmdLine state

  // MeshConv's CmdLine grammar

  //! brief Match and set verbose switch (i.e., verbose or quiet output)
  struct verbose :
         tk::grm::process_cmd_switch< use< kw::verbose >, tag::verbose > {};

  //! brief Match and set reorder switch (i.e., reorder mesh nodes or not)
  struct reorder :
         tk::grm::process_cmd_switch< use< kw::reorder >, tag::reorder > {};

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
                     reorder,
                     help,
                     helpkw,
                     io< use< kw::input >, tag::input >,
                     io< use< kw::output >, tag::output > > {};

  //! \brief Grammar entry point: parse keywords until end of string
  struct read_string :
         tk::grm::read_string< keywords > {};

} // cmd::
} // meshconv::

#endif // MeshConvCmdLineGrammar_h
