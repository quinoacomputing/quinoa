// *****************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/InputDeck.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Random number generator test suite input deck
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the control file parsing of the random
     number generator test suite, RNGTest.
*/
// *****************************************************************************
#ifndef RNGTestInputDeck_h
#define RNGTestInputDeck_h

#include <brigand/algorithms/for_each.hpp>

#include "NoWarning/set.h"

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


    //! \brief Constructor: set all defaults
    //! \param[in] cl Previously parsed and store command line
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type.
    InputDeck( const CmdLine& cl = {} ) {
      // Set previously parsed command line
      set< tag::cmd >( cl );
      // Initialize help
      const auto& ctrinfoFill = tk::ctr::Info( get< tag::cmd, tag::ctrinfo >() );
      brigand::for_each< keywords >( ctrinfoFill );
    }

    /** @name Pack/Unpack: Serialize InputDeck object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
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
    friend void operator|( PUP::er& p, InputDeck& i ) { i.pup(p); }
    //@}
};

} // ctr::
} // rngtest::

#endif // RNGTestInputDeck_h
