// *****************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/CmdLine.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     MeshConv's command line definition
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the command-line parsing of the mesh
     converter, MeshConv.
*/
// *****************************************************************************
#ifndef MeshConvCmdLine_h
#define MeshConvCmdLine_h

#include <string>

#include <brigand/algorithms/for_each.hpp>

#include "Macro.hpp"
#include "Keywords.hpp"
#include "HelpFactory.hpp"
#include "MeshConv/Types.hpp"

namespace meshconv {
//! Mesh converter control facilitating user input to internal data transfer
namespace ctr {

//! Member data for tagged tuple
using CmdLineMembers = brigand::list<
    tag::io,         ios
  , tag::verbose,    bool
  , tag::chare,      bool
  , tag::reorder,    bool
  , tag::help,       bool
  , tag::quiescence, bool
  , tag::trace,      bool
  , tag::version,    bool
  , tag::license,    bool
  , tag::cmdinfo,    tk::ctr::HelpFactory
  , tag::ctrinfo,    tk::ctr::HelpFactory
  , tag::helpkw,     tk::ctr::HelpKw
  , tag::error,      std::vector< std::string >
>;

//! \brief CmdLine is a TaggedTuple specialized to MeshConv
//! \details The stack is a tagged tuple, a hierarchical heterogeneous data
//!    structure where all parsed information is stored.
//! \see Base/TaggedTuple.h
//! \see Control/MeshConv/Types.h
class CmdLine : public tk::TaggedTuple< CmdLineMembers > {

  public:
    //! \brief MeshConv command-line keywords
    //! \see tk::grm::use and its documentation
    using keywords = tk::cmd_keywords< kw::verbose
                                     , kw::charestate
                                     , kw::help
                                     , kw::helpkw
                                     , kw::input
                                     , kw::output
                                     , kw::reorder_cmd
                                     , kw::quiescence
                                     , kw::trace
                                     , kw::version
                                     , kw::license
                                     >;

    //! \brief Constructor: set defaults.
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type. While there is a
    //!   ctrinfo parameter, it is unused here, since meshconv does not have a
    //!   control file parser.
    //! \see walker::ctr::CmdLine
    CmdLine() {
      get< tag::verbose >() = false; // Use quiet output by default
      get< tag::chare >() = false; // No chare state output by default
      get< tag::reorder >() = false; // Do not reorder by default
      get< tag::trace >() = true; // Output call and stack trace by default
      get< tag::version >() = false; // Do not display version info by default
      get< tag::license >() = false; // Do not display license info by default
      // Initialize help: fill from own keywords
      brigand::for_each< keywords::set >( tk::ctr::Info(get<tag::cmdinfo>()) );
    }

    /** @name Pack/Unpack: Serialize CmdLine object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) { tk::TaggedTuple< CmdLineMembers >::pup(p); }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c CmdLine object reference
    friend void operator|( PUP::er& p, CmdLine& c ) { c.pup(p); }
    //@}
};

} // ctr::
} // meshconv::

#endif // MeshConvCmdLine_h
