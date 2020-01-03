// *****************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/InputDeck.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Random number generator test suite input deck
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the control file parsing of the random
     number generator test suite, RNGTest.
*/
// *****************************************************************************
#ifndef RNGTestInputDeck_h
#define RNGTestInputDeck_h

#include <brigand/algorithms/for_each.hpp>

#include "NoWarning/set.hpp"

#include "QuinoaConfig.hpp"
#include "TaggedTuple.hpp"
#include "RNGTest/CmdLine/CmdLine.hpp"

namespace rngtest {
namespace ctr {

//! Member data for tagged tuple
using InputDeckMembers = brigand::list<
    tag::cmd,        CmdLine
  , tag::title,      kw::title::info::expect::type
  , tag::selected,   selects
  , tag::param,      parameters
  , tag::error,      std::vector< std::string >
>;

//! \brief InputDeck is a TaggedTuple specialized to RNGTest
//! \details The stack is a tagged tuple, a hierarchical heterogeneous data
//!    structure where all parsed information is stored.
//! \see Base/TaggedTuple.h
//! \see Control/RNGTest/Types.h
class InputDeck : public tk::TaggedTuple< InputDeckMembers > {

  public:
    //! \brief RNGTest input deck keywords
    //! \see tk::grm::use and its documentation
    using keywords = brigand::set< kw::title
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
                                 , kw::gaussianmv_method
                                 , kw::beta_method
                                 , kw::standard
                                 , kw::accurate
                                 , kw::boxmuller
                                 , kw::boxmuller2
                                 , kw::icdf
                                 #endif
                                 , kw::seed
                                 , kw::r123_threefry
                                 , kw::r123_philox
                                 , kw::gamma_method
                                 , kw::gnorm
                                 , kw::gnorm_accurate
                                 >;

    //! Set of tags to ignore when printing this InputDeck
    using ignore = CmdLine::ignore;

    //! \brief Constructor: set all defaults
    //! \param[in] cl Previously parsed and store command line
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type.
    explicit InputDeck( const CmdLine& cl = {} ) {
      // Set previously parsed command line
      get< tag::cmd >() = cl;
      // Initialize help
      const auto& ctrinfoFill = tk::ctr::Info( get< tag::cmd, tag::ctrinfo >() );
      brigand::for_each< keywords >( ctrinfoFill );
    }

    /** @name Pack/Unpack: Serialize InputDeck object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) { tk::TaggedTuple< InputDeckMembers >::pup(p); }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i InputDeck object reference
    friend void operator|( PUP::er& p, InputDeck& i ) { i.pup(p); }
    //@}
};

} // ctr::
} // rngtest::

#endif // RNGTestInputDeck_h
