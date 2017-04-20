// *****************************************************************************
/*!
  \file      src/Control/HelpFactory.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Command-line and input deck help factory
  \details   This file contains some types that facilitate the generation of
     on-screen help.
*/
// *****************************************************************************
#ifndef HelpFactory_h
#define HelpFactory_h

#include <boost/optional.hpp>

#include "PUPUtil.h"
#include "Factory.h"
#include "Has.h"

namespace tk {
namespace ctr {

//! \brief Keyword information bundle
//! \details This bundle contains the information that is used to display
//!    on-screen help on all command-line arguments and control file keywords
//!    for an exectuable. This struct is stored in a container that associates
//!    keywords (used by a grammar and parser) to this struct. The container, an
//!    runtime, std::map, is filled by the CmdLine and InputDeck objects'
//!    constructors by one or more boost::mpl::for_each which loops through the
//!    set of all keywords used in a grammar. The maps are stored in the CmdLine
//!    and InputDeck objects (which are tagged tuples) and thus can be migrated
//!    through the network, thus the Charm++ parck/unpack routines are defined.
//! \see Info functor used to fill the std::maps
//! \author J. Bakosi
struct KeywordInfo {
  std::string shortDescription;           //!< Short description
  std::string longDescription;            //!< Long description
  boost::optional< std::string > alias;   //!< Keyword alias
  boost::optional< std::string > expt;    //!< Expected type description
  boost::optional< std::string > choices; //!< Expected choices descr.

  /** @name Pack/Unpack: Serialize KeywordInfo object for Charm++ */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \author J. Bakosi
  void pup( PUP::er& p ) {
    p | shortDescription;
    p | longDescription;
    p | alias;
    p | expt;
    p | choices;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] info KeywordInfo object reference
  //! \author J. Bakosi
  friend void operator|( PUP::er& p, KeywordInfo& info ) { info.pup(p); }
  ///@}
};

//! \brief A typedef for associating a keyword-string with its associated
//!   information stored in a KeywordInfo struct
//! \author J. Bakosi
using HelpFactory = std::map< std::string, KeywordInfo >;

//! \brief Help bundle on a single keyword
//! \details This is used for delivering help on a single keyword. This struct
//!    also differentiates between command-line arguments and control file
//!    keywords.
//! \author J. Bakosi
struct HelpKw {
  HelpFactory::key_type keyword;        //!< Keyword string
  HelpFactory::mapped_type info;        //!< Keyword information
  bool cmd;                             //!< True if command-line keyword

  /** @name Pack/Unpack: Serialize HelpKw object for Charm++ */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \author J. Bakosi
  void pup( PUP::er& p ) { p|keyword; p|info; p|cmd; }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] h HelpKw object reference
  //! \author J. Bakosi
  friend void operator|( PUP::er& p, HelpKw& h ) { h.pup(p); }
  ///@}
};

//! \brief Function object for filling a HelpFactory (std::map) with keywords
//!   and their associated information bundle
//! \details This struct is used as a functor to loop through a set of keywords
//!   at compile-time and generate code for filling up the std::map.
//! \author J. Bakosi
struct Info {
  //! Store reference to map we are filling
  tk::ctr::HelpFactory& m_factory;
  //! Constructor: store reference to map to fill
  Info( tk::ctr::HelpFactory& factory ) : m_factory( factory ) {}
  //! \brief Function call operator templated on the type that does the filling
  template< typename U > void operator()( U ) {
    m_factory[ U::string() ] = { U::shortDescription(),
                                 U::longDescription(),
                                 U::alias(),
                                 U::expt(),
                                 U::choices() };
  }
};

} // ctr::
} // tk::

#endif // HelpFactory_h
