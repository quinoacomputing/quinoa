//******************************************************************************
/*!
  \file      src/Control/Breeze/InputDeck/InputDeck.h
  \author    J. Bakosi
  \date      Wed 08 Apr 2015 08:43:31 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Breeze's input deck definition
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the control file parsing of the
     computational fluid dynamics tool, Breeze.
*/
//******************************************************************************
#ifndef BreezeInputDeck_h
#define BreezeInputDeck_h

#include <limits>

#include <boost/mpl/set.hpp>
#include <boost/mpl/for_each.hpp>

#include <Control.h>
#include <Breeze/CmdLine/CmdLine.h>
#include <Breeze/Components.h>

namespace breeze {
namespace ctr {

//! \brief InputDeck : Control< specialized to Breeze >, see Types.h,
//! \details The stack is a tagged tuple, a hierarchical heterogeneous data
//!    structure where all parsed information is stored.
//! \see Base/TaggedTuple.h
//! \see Control/Breeze/Types.h
//! \author J. Bakosi
class InputDeck :
  public tk::Control< // tag           type
                      tag::title,      kw::title::info::expect::type,
                      tag::selected,   selects,
                      tag::discr,      discretization,
                      tag::component,  ncomps,
                      tag::interval,   intervals,
                      tag::cmd,        CmdLine,
                      tag::param,      parameters,
                      tag::stat,       std::vector< tk::ctr::Product >,
                      tag::pdf,        std::vector< tk::ctr::Probability >,
                      tag::error,      std::vector< std::string > > {

  public:
    //! \brief Breeze input deck keywords
    //! \details Since there are more than 20 and boost::mpl only allows maxium
    //!   20 items in a set by default (and I don't want to mess with
    //!   preprocessor-generated boost::mpl headers), the whole set is broken up
    //!   into several sets each containing 20 keywords.
    //! \see tk::grm::use and its documentation
    using keywords1 = boost::mpl::set< kw::precision
                                     , kw::end
                                     , kw::depvar
                                     , kw::title
                                     , kw::statistics
                                     , kw::interval
                                     , kw::pdfs
                                     , kw::pdf_filetype
                                     , kw::pdf_policy
                                     , kw::pdf_centering
                                     , kw::txt_float_format
                                     , kw::npar
                                     , kw::nstep
                                     , kw::term
                                     , kw::dt
                                     , kw::ttyi
                                     , kw::rngs
                                     , kw::rng
                                     , kw::hommix
                                     , kw::homhydro
                                     >;
    using keywords2 = boost::mpl::set< kw::homrt
                                     , kw::spinsflow
                                     , kw::hydro_slm
                                     , kw::freq_gamma
                                     , kw::hydro_slm
                                     , kw::hydro_glm
                                     , kw::beta
                                     , kw::mix_iem
                                     , kw::mix_iecm
                                     , kw::mix_dir
                                     , kw::mix_gendir
                                     , kw::mixrate_gamma
                                     , kw::pos_inviscid
                                     , kw::pos_viscous
                                     >;
    using keywords3 = boost::mpl::set< kw::rngsse_gm19
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
                                     >;
    using keywords4 = boost::mpl::set< kw::seed
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
                                     , kw::standard
                                     , kw::accurate
                                     , kw::boxmuller
                                     , kw::boxmuller2
                                     , kw::icdf
                                     #endif
                                     >;

    //! \brief Constructor: set defaults
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type.
    //! \author J. Bakosi
    InputDeck() {
      // Default discretization parameters
      set< tag::discr, tag::npar >( 1 );
      set< tag::discr, tag::nstep >
         ( std::numeric_limits< kw::nstep::info::expect::type >::max() );
      set< tag::discr, tag::term >( 1.0 );
      set< tag::discr, tag::dt >( 0.5 );
      // Default txt floating-point output precision in digits
      set< tag::prec, tag::stat >( std::cout.precision() );
      set< tag::prec, tag::pdf >( std::cout.precision() );
      // Default intervals
      set< tag::interval, tag::tty >( 1 );
      set< tag::interval, tag::dump >( 1 );
      set< tag::interval, tag::stat >( 1 );
      set< tag::interval, tag::pdf >( 1 );
      set< tag::interval, tag::glob >( 1 );
      // Default simplified Langevin hydro model parameters
      set< tag::param, tag::slm, tag::c0 >( 2.1 );
      // Default generalized Langevin hydro model parameters
      set< tag::param, tag::slm, tag::c0 >( 2.1 );
      // Default requested statistics
      set< tag::stat >( std::vector< tk::ctr::Product >() );
      // Initialize help: fill own keywords
      const auto& ctrinfoFill = tk::ctr::Info( get< tag::cmd, tag::ctrinfo >() );
      boost::mpl::for_each< keywords1 >( ctrinfoFill );
      boost::mpl::for_each< keywords2 >( ctrinfoFill );
      boost::mpl::for_each< keywords3 >( ctrinfoFill );
      boost::mpl::for_each< keywords4 >( ctrinfoFill );
    }

    /** @name Pack/Unpack: Serialize InputDeck object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[inout] p Charm++'s PUP::er serializer object reference
    //! \author J. Bakosi
    void pup( PUP::er& p ) {
      tk::Control< tag::title,      kw::title::info::expect::type,
                   tag::selected,   selects,
                   tag::discr,      discretization,
                   tag::component,  ncomps,
                   tag::interval,   intervals,
                   tag::cmd,        CmdLine,
                   tag::param,      parameters,
                   tag::stat,       std::vector< tk::ctr::Product >,
                   tag::pdf,        std::vector< tk::ctr::Probability >,
                   tag::error,      std::vector< std::string > >::pup(p);
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[inout] p Charm++'s PUP::er serializer object reference
    //! \param[inout] i InputDeck object reference
    //! \author J. Bakosi
    friend void operator|( PUP::er& p, InputDeck& i ) { i.pup(p); }
    //@}
};

} // ctr::
} // breeze::

#endif // BreezeInputDeck_h
