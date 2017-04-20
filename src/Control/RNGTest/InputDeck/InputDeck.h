// *****************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/InputDeck.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Random number generator test suite input deck
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the control file parsing of the random
     number generator test suite, RNGTest.
*/
// *****************************************************************************
#ifndef RNGTestInputDeck_h
#define RNGTestInputDeck_h

#include "NoWarning/set.h"
#include "NoWarning/for_each.h"

#include "Control.h"
#include "QuinoaConfig.h"
#include "RNGTest/CmdLine/CmdLine.h"

namespace rngtest {
namespace ctr {

//! \brief InputDeck : Control< specialized to RNGTest >, see Types.h
//! \details The stack is a tagged tuple, a hierarchical heterogeneous data
//!    structure where all parsed information is stored.
//! \see Base/TaggedTuple.h
//! \see Control/RNGTest/Types.h
//! \author J. Bakosi
class InputDeck : public tk::Control<
                    // tag           type
                    tag::title,      kw::title::info::expect::type,
                    tag::selected,   selects,
                    tag::io,         ios,
                    tag::cmd,        CmdLine,
                    tag::param,      parameters,
                    tag::error,      std::vector< std::string > > {

  public:
    //! \brief RNGTest input deck keywords
    //! \details Since there are more than 20 and boost::mpl only allows maxium
    //!   20 items in a set by default (and I don't want to mess with
    //!   preprocessor-generated boost::mpl headers), the whole set is broken up
    //!   into several sets each containing 20 keywords.
    //! \see tk::grm::use and its documentation
    using keywords1 = boost::mpl::set< kw::title
                                     , kw::end
                                     , kw::smallcrush
                                     , kw::crush
                                     , kw::bigcrush
                                     , kw::cja
                                     , kw::cja_accurate
                                     #ifdef HAS_RNGSSE2
                                     , kw::rngsse_gm19
                                     , kw::rngsse_gm29
                                     , kw::rngsse_gm31
                                     , kw::rngsse_gm55
                                     , kw::rngsse_gm61
                                     , kw::rngsse_gq581
                                     , kw::rngsse_gq583
                                     , kw::rngsse_gq584
                                     , kw::rngsse_mt19937
                                     , kw::rngsse_lfsr113
                                     , kw::rngsse_mrg32k3a
                                     , kw::seqlen
                                     #endif
                                     >;
    using keywords2 = boost::mpl::set< kw::seed
                                     #ifdef HAS_MKL
                                     , kw::mkl_mcg31
                                     , kw::mkl_r250
                                     , kw::mkl_mrg32k3a
                                     , kw::mkl_mcg59
                                     , kw::mkl_wh
                                     , kw::mkl_mt19937
                                     , kw::mkl_mt2203
                                     , kw::mkl_sfmt19937
                                     , kw::mkl_sobol
                                     , kw::mkl_niederr
                                     , kw::mkl_nondeterm
                                     , kw::uniform_method
                                     , kw::gaussian_method
                                     , kw::beta_method
                                     , kw::standard
                                     , kw::accurate
                                     , kw::boxmuller
                                     , kw::boxmuller2
                                     , kw::icdf
                                     #endif
                                     >;
    using keywords3 = boost::mpl::set< kw::r123_threefry
                                     , kw::r123_philox
                                     >;


    //! \brief Constructor: set all defaults
    //! \param[in] cl Previously parsed and store command line
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type.
    //! \author J. Bakosi
    InputDeck( const CmdLine& cl = {} ) {
      // Set previously parsed command line
      set< tag::cmd >( cl );
      // Initialize help
      const auto& ctrinfoFill = tk::ctr::Info( get< tag::cmd, tag::ctrinfo >() );
      boost::mpl::for_each< keywords1 >( ctrinfoFill );
      boost::mpl::for_each< keywords2 >( ctrinfoFill );
      boost::mpl::for_each< keywords3 >( ctrinfoFill );
    }

    /** @name Pack/Unpack: Serialize InputDeck object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \author J. Bakosi
    void pup( PUP::er& p ) {
      tk::Control< tag::title,      kw::title::info::expect::type,
                   tag::selected,   selects,
                   tag::io,         ios,
                   tag::cmd,        CmdLine,
                   tag::param,      parameters,
                   tag::error,      std::vector< std::string > >::pup(p);
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i InputDeck object reference
    //! \author J. Bakosi
    friend void operator|( PUP::er& p, InputDeck& i ) { i.pup(p); }
    //@}
};

} // ctr::
} // rngtest::

#endif // RNGTestInputDeck_h
