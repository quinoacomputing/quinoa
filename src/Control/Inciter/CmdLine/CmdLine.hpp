// *****************************************************************************
/*!
  \file      src/Control/Inciter/CmdLine/CmdLine.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's command line definition
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the command-line parsing of the
     computational shock hydrodynamics tool, Inciter.
*/
// *****************************************************************************
#ifndef InciterCmdLine_h
#define InciterCmdLine_h

#include <string>

#include <brigand/algorithms/for_each.hpp>

#include "QuinoaConfig.hpp"
#include "TaggedTuple.hpp"
#include "PrintUtil.hpp"
#include "Tags.hpp"
#include "HelpFactory.hpp"
#include "Keywords.hpp"
#include "Inciter/OutVar.hpp"

namespace inciter {
//! Inciter control facilitating user input to internal data transfer
namespace ctr {

//! IO parameters storage
using ios = tk::TaggedTuple< brigand::list<
    tag::nrestart,  int                             //!< Number of restarts
  , tag::control,   kw::control::info::expect::type //!< Control filename
  , tag::input,     kw::input::info::expect::type   //!< Input filename
  , tag::output,    kw::output::info::expect::type  //!< Output filename
    //! Refined output (output field data on a refined mesh)
  , tag::refined,   kw::refined::info::expect::type
  , tag::screen,    kw::screen::info::expect::type  //!< Screen output filename
    //! List of side sets to save as field output
  , tag::surface,   std::vector< kw::sideset::info::expect::type >
    //! Diagnostics filename
  , tag::diag,      kw::diagnostics_cmd::info::expect::type
  , tag::particles, std::string                     //!< Particles filename
  , tag::outvar,    std::vector< OutVar >           //!< Output variables
  , tag::restart,   kw::restart::info::expect::type //!< Restart dirname
> >;

//! Member data for tagged tuple
using CmdLineMembers = brigand::list<
    tag::io,             ios
  , tag::virtualization, kw::virtualization::info::expect::type
  , tag::verbose,        bool
  , tag::chare,          bool
  , tag::nonblocking,    bool
  , tag::benchmark,      bool
  , tag::feedback,       bool
  , tag::help,           bool
  , tag::helpctr,        bool
  , tag::quiescence,     bool
  , tag::trace,          bool
  , tag::version,        bool
  , tag::license,        bool
  , tag::cmdinfo,        tk::ctr::HelpFactory
  , tag::ctrinfo,        tk::ctr::HelpFactory
  , tag::helpkw,         tk::ctr::HelpKw
  , tag::error,          std::vector< std::string >
  , tag::lbfreq,         kw::lbfreq::info::expect::type
  , tag::rsfreq,         kw::rsfreq::info::expect::type
>;

//! \brief CmdLine : Control< specialized to Inciter >
//! \details The stack is a tagged tuple
//! \see Base/TaggedTuple.h
//! \see Control/Inciter/Types.h
class CmdLine : public tk::TaggedTuple< CmdLineMembers > {

  public:
    //! \brief Inciter command-line keywords
    //! \see tk::grm::use and its documentation
    using keywords = tk::cmd_keywords< kw::verbose
                                     , kw::charestate
                                     , kw::nonblocking
                                     , kw::benchmark
                                     , kw::feedback
                                     , kw::virtualization
                                     , kw::help
                                     , kw::helpctr
                                     , kw::helpkw
                                     , kw::control
                                     , kw::input
                                     , kw::output
                                     , kw::screen
                                     , kw::restart
                                     , kw::diagnostics_cmd
                                     , kw::quiescence
                                     , kw::lbfreq
                                     , kw::rsfreq
                                     , kw::trace
                                     , kw::version
                                     , kw::license
                                     >;

    //! Set of tags to ignore when printing this CmdLine
    using ignore =
      brigand::set< tag::cmdinfo
                  , tag::ctrinfo
                  , tag::helpkw >;

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
      get< tag::io, tag::nrestart >() = 0;
      get< tag::io, tag::output >() = "out";
      get< tag::io, tag::refined >() = false;
      get< tag::io, tag::screen >() =
        tk::baselogname( tk::inciter_executable() );
      get< tag::io, tag::diag >() = "diag";
      get< tag::io, tag::particles >() = "track.h5part";
      get< tag::io, tag::restart >() = "restart";
      get< tag::virtualization >() = 0.0;
      get< tag::verbose >() = false; // Quiet output by default
      get< tag::chare >() = false; // No chare state output by default
      get< tag::nonblocking>() = false; // Blocking migration by default
      get< tag::benchmark >() = false; // No benchmark mode by default
      get< tag::feedback >() = false; // No detailed feedback by default
      get< tag::lbfreq >() = 1; // Load balancing every time-step by default
      get< tag::rsfreq >() = 1000;// Chkpt/restart after this many time steps
      get< tag::trace >() = true; // Output call and stack trace by default
      get< tag::version >() = false; // Do not display version info by default
      get< tag::license >() = false; // Do not display license info by default
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

    //! Compute and return log file name
    //! \param[in] def Default log file name (so we don't mess with user's)
    //! \param[in] nrestart Number of times restarted
    //! \return Log file name
    std::string logname( const std::string& def, int nrestart ) const {
      if (get< tag::io, tag::screen >() != def)
        return get< tag::io, tag::screen >();
      else
        return tk::logname( tk::inciter_executable(), nrestart );
    }
};

} // ctr::
} // inciter::

#endif // InciterCmdLine_h
