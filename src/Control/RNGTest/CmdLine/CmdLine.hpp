// *****************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/CmdLine.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     RNGTest's command line
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the command-line parsing of the
     random number generator test suite, RNGTest.
*/
// *****************************************************************************
#ifndef RNGTestCmdLine_h
#define RNGTestCmdLine_h

#include <string>

#include <brigand/algorithms/for_each.hpp>

#include "QuinoaConfig.hpp"
#include "TaggedTuple.hpp"
#include "HelpFactory.hpp"
#include "Keywords.hpp"
#include "RNGTest/Types.hpp"

namespace rngtest {
//! RNGTest control facilitating user input to internal data transfer
namespace ctr {

//! Member data for tagged tuple
using CmdLineMembers = brigand::list<
    tag::io,         ios
  , tag::verbose,    bool
  , tag::chare,      bool
  , tag::help,       bool
  , tag::helpctr,    bool
  , tag::quiescence, bool
  , tag::trace,      bool
  , tag::version,    bool
  , tag::license,    bool
  , tag::cmdinfo,    tk::ctr::HelpFactory
  , tag::ctrinfo,    tk::ctr::HelpFactory
  , tag::helpkw,     tk::ctr::HelpKw
  , tag::error,      std::vector< std::string >
>;

//! CmdLine is a TaggedTuple specialized to RNGTest
class CmdLine : public tk::TaggedTuple< CmdLineMembers > {

  public:
    //! RNGTest command-line keywords
    //! \see tk::grm::use and its documentation
    using keywords = tk::cmd_keywords< kw::verbose
                                     , kw::charestate
                                     , kw::control
                                     , kw::help
                                     , kw::helpctr
                                     , kw::helpkw
                                     , kw::screen
                                     , kw::quiescence
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
      get< tag::io, tag::screen >() =
        tk::baselogname( tk::rngtest_executable() );
      get< tag::verbose >() = false; // Use quiet output by default
      get< tag::chare >() = false; // No chare state output by default
      get< tag::trace >() = true; // Output call and stack trace by default
      get< tag::version >() = false; // Do not display version info by default
      get< tag::license >() = false; // Do not display license info by default
      // Initialize help: fill from own keywords + add map passed in
      brigand::for_each< keywords::set >( tk::ctr::Info( get<tag::cmdinfo>()) );
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
        return tk::logname( tk::rngtest_executable(), nrestart );
    }
};

} // ctr::
} // rngtest::

#endif // RNGTestCmdLine_h
