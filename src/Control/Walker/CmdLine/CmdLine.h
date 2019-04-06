// *****************************************************************************
/*!
  \file      src/Control/Walker/CmdLine/CmdLine.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Walker's command line
  \details   Walker's command line
*/
// *****************************************************************************
#ifndef WalkerCmdLine_h
#define WalkerCmdLine_h

#include <string>

#include <brigand/algorithms/for_each.hpp>

#include "HelpFactory.h"
#include "Keywords.h"
#include "Walker/Types.h"

namespace walker {
//! Walker control facilitating user input to internal data transfer
namespace ctr {

//! Member data for tagged tuple
using CmdLineMembers = brigand::list<
    tag::io,             ios
  , tag::virtualization, kw::virtualization::info::expect::type
  , tag::verbose,        bool
  , tag::chare,          bool
  , tag::help,           bool
  , tag::helpctr,        bool
  , tag::quiescence,     bool
  , tag::trace,          bool
  , tag::cmdinfo,        tk::ctr::HelpFactory
  , tag::ctrinfo,        tk::ctr::HelpFactory
  , tag::helpkw,         tk::ctr::HelpKw
  , tag::error,          std::vector< std::string >
>;

//! CmdLine : Control< specialized to Walker >, see Types.h
class CmdLine : public tk::TaggedTuple< CmdLineMembers > {

  public:
    //! Walker command-line keywords
    using keywords = tk::cmd_keywords< kw::verbose
                                     , kw::charestate
                                     , kw::virtualization
                                     , kw::help
                                     , kw::helpctr
                                     , kw::helpkw
                                     , kw::control
                                     , kw::pdf
                                     , kw::stat
                                     , kw::trace
                                     , kw::quiescence
                                     >;

    //! \brief Constructor: set all defaults.
    //! \param[in] ctrinfo std::map of control file keywords and their info
    //!  \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type. The ctrinfo map
    //!   argument is optional. If not given, it is an empty std::map
    //!   constructed in-place and affects nothing. If given, it contains the
    //!   control file keywords, all of which are moved into the relevant slot
    //!   (tag::ctrinfo). This allows constructing, e.g., a CmdLine object both
    //!   with and without this information in place, which are both used at
    //!   different stages of the execution. For example, because the
    //!   command-line is parsed very early on during runtime while the input
    //!   deck is only parsed much later, the control-file keywords and their
    //!   information (owned by and generated by the input deck and its
    //!   constructor) is not yet available when the CmdLine object is
    //!   constructed. However, during command-line parsing it is still possible
    //!   to request information on a control file keyword, so it must be
    //!   available. The input deck is where all parsed information goes during
    //!   control file parsing and is stored at global scope, see, e.g.,
    //!   walker::g_inputdeck. This global-scope (still namespace-scope), input
    //!   deck object is thus created before command-line parsing. The input
    //!   deck object's constructor (working only on type information, available
    //!   at compile-time, of all the control file keywords) creates a run-time
    //!   map. This is a run-time map, but available before main() starts,
    //!   because it is const and it is initialized as a global-scope map. This
    //!   map is then passed in here as ctrinfo, and its contents inserted into
    //!   the CmdLine object, making the control-file keywords and their info
    //!   available during command-line parsing. Since the input deck stack
    //!   contains a copy of the command-line stack, the command-line stack must
    //!   be possible to be instantiated without passing the ctrinfo map,
    //!   otherwise it would be a mutual dependency.
    // cppcheck-suppress noExplicitConstructor
    CmdLine( tk::ctr::HelpFactory ctrinfo = tk::ctr::HelpFactory() ) {
      get< tag::io, tag::output >() = "out";
      get< tag::io, tag::pdf >() = "pdf";
      get< tag::io, tag::stat >() = "stat.txt";
      get< tag::virtualization >() = 0.0;
      get< tag::verbose >() = false; // Quiet output by default
      get< tag::chare >() = false; // No chare state output by default
      get< tag::trace >() = true; // Output call and stack trace by default
      // Initialize help: fill from own keywords + add map passed in
      brigand::for_each< keywords::set >( tk::ctr::Info(get<tag::cmdinfo>()) );
      get< tag::ctrinfo >() = std::move( ctrinfo );
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
} // walker::

#endif // WalkerCmdLine_h
