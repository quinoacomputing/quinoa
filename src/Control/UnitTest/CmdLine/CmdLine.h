// *****************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/CmdLine.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     UnitTest's command line
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the command-line parsing of the unit test
     suite, UnitTest.
*/
// *****************************************************************************
#ifndef UnitTestCmdLine_h
#define UnitTestCmdLine_h

#include <string>

#include <brigand/algorithms/for_each.hpp>

#include "NoWarning/set.h"

#include "Macro.h"
#include "Control.h"
#include "HelpFactory.h"
#include "Keywords.h"
#include "UnitTest/Types.h"

namespace unittest {
//! UnitTest control facilitating user input to internal data transfer
namespace ctr {

//! CmdLine : Control< specialized to UnitTest >
//! \details The stack is a tagged tuple
//! \see Base/TaggedTuple.h
//! \see Control/UnitTest/Types.h
class CmdLine : public tk::Control<
                  // tag            type
                  tag::verbose,     bool,
                  tag::chare,       bool,
                  tag::help,        bool,
                  tag::helpctr,     bool,
                  tag::quiescence,  bool,
                  tag::cmdinfo,     tk::ctr::HelpFactory,
                  tag::ctrinfo,     tk::ctr::HelpFactory,
                  tag::helpkw,      tk::ctr::HelpKw,
                  tag::group,       std::string,
                  tag::error,       std::vector< std::string > > {
  public:
    //! \brief UnitTest command-line keywords
    //! \see tk::grm::use and its documentation
    using keywords = brigand::set< kw::verbose
                                 , kw::charestate
                                 , kw::help
                                 , kw::helpkw
                                 , kw::group
                                 , kw::quiescence
                                 >;

    //! \brief Constructor: set defaults.
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type. While there is a
    //!   ctrinfo parameter, it is unused here, since unittest does not have a
    //!   control file parser.
    //! \see walker::ctr::CmdLine
    CmdLine() {
      set< tag::verbose >( false ); // Use quiet output by default
      set< tag::chare >( false ); // No chare state output by default
      // Initialize help: fill from own keywords
      brigand::for_each< keywords >( tk::ctr::Info( get< tag::cmdinfo >() ) );
    }

    /** @name Pack/Unpack: Serialize CmdLine object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) {
      tk::Control< tag::verbose,    bool,
                   tag::chare,      bool,
                   tag::help,       bool,
                   tag::helpctr,    bool,
                   tag::quiescence, bool,
                   tag::cmdinfo,    tk::ctr::HelpFactory,
                   tag::ctrinfo,    tk::ctr::HelpFactory,
                   tag::helpkw,     tk::ctr::HelpKw,
                   tag::group,      std::string,
                   tag::error,      std::vector< std::string > >::pup(p);
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c CmdLine object reference
    friend void operator|( PUP::er& p, CmdLine& c ) { c.pup(p); }
    //@}
};

} // ctr::
} // unittest::

#endif // UnitTestCmdLine_h
