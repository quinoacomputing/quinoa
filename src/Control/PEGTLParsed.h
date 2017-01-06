// *****************************************************************************
/*!
  \file      src/Control/PEGTLParsed.h
  \author    J. Bakosi
  \date      Fri 06 Jan 2017 01:45:04 PM MST
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Class to equip parsed classes with PEGTL instruments
  \details   Class to equip parsed classes with PEGTL instruments. This is used
    to track the parser's location so that we can detect the location of errors
    during parsing.
*/
// *****************************************************************************
#ifndef PEGTLParsed_h
#define PEGTLParsed_h

namespace tk {
namespace ctr {

struct unused {};

//! \brief PEGTLParsed to equip a PEGTL stack, used to store heterogeneous
//!   information during parsing, with location information
//! \details This is used to track the parser's location so that we can detect
//!   the location of errors
//! \author J. Bakosi
template< class Parsed, class cmdtag = unused, class Cmd = unused >
class PEGTLParsed : public Parsed {

  public:
    //! \brief Constructor
    //! \param[in] ctrinfo std::map of control file keywords and their info
    //! \details The ctrinfo map argument is optional. If not given, it is an
    //!    empty std::map constructed in-place and affects nothing. If given, it
    //!    contains the control file keywords, all of which are moved into the
    //!    relevant slot (tag::ctrinfo) into Parsed. This allows constructing,
    //!    e.g., a CmdLine object both with and without this information in
    //!    place, which are both used at different stages of the execution. For
    //!    example, because the command-line is parsed very early on during
    //!    runtime while the input deck is only parsed much later, the
    //!    control-file keywords and their information (owned by and generated
    //!    by the input deck and its constructor) is not yet available when the
    //!    CmdLine object is constructed. However, during command-line parsing
    //!    it is still possible to request information on a control file
    //!    keyword, so it must be available. The input deck is where all
    //!    parsed information goes during control file parsing and is stored at
    //!    global scope, see, e.g., walker::g_inputdeck. This global-scope
    //!    (still namespace-scope), input deck object is thus created before
    //!    command-line parsing. The input deck object's constructor (working
    //!    only on type information, available at compile-time, of all the
    //!    control file keywords) creates a run-time map. This is a run-time
    //!    map, but available before main() starts, because it is const and it
    //!    is initialized as a global-scope map. This map is then passed in here
    //!    as ctrinfo, and its contents inserted into the CmdLine object, making
    //!    the control-file keywords and their info available during
    //!    command-line parsing. Since the input deck stack contains a copy of
    //!    the command-line stack, the command-line stack must be possible to be
    //!    instantiated without passing the ctrinfo map, otherwise it would be a
    //!    mutual dependency.
    //! \author J. Bakosi
    explicit PEGTLParsed( HelpFactory ctrinfo = HelpFactory() ) :
      Parsed( ctrinfo ) {}

    //! \brief Constructor setting command line
    //! \details This constructor is used for control file parsers whose state
    //!   includes the previously parsed command line.
    //! \author J. Bakosi
    explicit PEGTLParsed( const Cmd& cl )
    { Parsed::template set< cmdtag >( cl ); }
};

} // ctr::
} // tk::

#endif // PEGTLParsed_h
