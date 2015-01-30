//******************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Sat 17 Jan 2015 06:53:25 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     UnitTest's command line
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the command-line parsing of the unit test
     suite, UnitTest.
*/
//******************************************************************************
#ifndef UnitTestCmdLine_h
#define UnitTestCmdLine_h

#include <string>

#include <boost/mpl/set.hpp>
#include <boost/mpl/for_each.hpp>

#include <Control.h>
#include <HelpFactory.h>
#include <Keywords.h>
#include <UnitTest/Types.h>

namespace unittest {
//! UnitTest control facilitating user input to internal data transfer
namespace ctr {

//! CmdLine : Control< specialized to UnitTest >
//! \details The stack is a tagged tuple
//! \see Base/TaggedTuple.h
//! \see Control/UnitTest/Types.h
//! \author J. Bakosi
class CmdLine : public tk::Control<
                  // tag        type
                  tag::verbose, bool,
                  tag::help,    bool,
                  tag::helpctr, bool,
                  tag::cmdinfo, tk::ctr::HelpFactory,
                  tag::ctrinfo, tk::ctr::HelpFactory,
                  tag::helpkw,  tk::ctr::HelpKw,
                  tag::error,   std::vector< std::string > > {
  public:
    //! \brief UnitTest command-line keywords
    //! \author J. Bakosi
    //! \see tk::grm::use and its documentation
    using keywords = boost::mpl::set< kw::verbose
                                    , kw::help
                                    , kw::helpkw
                                    >;

    //! \brief Constructor: set defaults.
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type. While there is a
    //!   ctrinfo parameter, it is unused here, since unittest does not have a
    //!   control file parser.
    //! \param[in] ctrinfo empty std::map, unused so that PEGTLParsed's
    //!    constructor can be kept generic that works with those executables
    //!    that do have both command-line argument and control file parsers
    //! \see walker::ctr::CmdLine
    CmdLine( tk::ctr::HelpFactory ctrnfo = tk::ctr::HelpFactory() ) {
      set< tag::verbose >( false ); // Use quiet output by default
      // Initialize help: fill from own keywords
      boost::mpl::for_each< keywords >( tk::ctr::Info( get< tag::cmdinfo >() ) );
    }

    /** @name Pack/Unpack: Serialize CmdLine object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[inout] p Charm++'s PUP::er serializer object reference
    //! \author J. Bakosi
    void pup( PUP::er& p ) {
      tk::Control< tag::verbose, bool,
                   tag::help,    bool,
                   tag::helpctr, bool,
                   tag::cmdinfo, tk::ctr::HelpFactory,
                   tag::ctrinfo, tk::ctr::HelpFactory,
                   tag::helpkw,  tk::ctr::HelpKw,
                   tag::error,   std::vector< std::string > >::pup(p);
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[inout] p Charm++'s PUP::er serializer object reference
    //! \param[inout] c CmdLine object reference
    //! \author J. Bakosi
    friend void operator|( PUP::er& p, CmdLine& c ) { c.pup(p); }
    //@}
};

} // ctr::
} // unittest::

#endif // UnitTestCmdLine_h
