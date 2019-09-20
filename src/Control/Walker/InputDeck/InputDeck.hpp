// *****************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/InputDeck.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Walker's input deck
  \details   Walker's input deck
*/
// *****************************************************************************
#ifndef WalkerInputDeck_h
#define WalkerInputDeck_h

#include <limits>
#include <iostream>

#include <brigand/algorithms/for_each.hpp>

#include "NoWarning/set.hpp"

#include "QuinoaConfig.hpp"
#include "TaggedTuple.hpp"
#include "HelpFactory.hpp"
#include "Walker/CmdLine/CmdLine.hpp"
#include "Walker/Components.hpp"

namespace walker {
namespace ctr {

//! Member data for tagged tuple
using InputDeckMembers = brigand::list<
    tag::title,      kw::title::info::expect::type
  , tag::selected,   selects
  , tag::discr,      discretization
  , tag::prec,       precision
  , tag::flformat,   floatformat
  , tag::component,  ncomps
  , tag::interval,   intervals
  , tag::cmd,        CmdLine
  , tag::param,      parameters
  , tag::stat,       std::vector< tk::ctr::Product >
  , tag::pdf,        std::vector< tk::ctr::Probability >
  , tag::error,      std::vector< std::string >
>;

//! InputDeck : Control< specialized to Walker >, see Types.h
class InputDeck : public tk::TaggedTuple< InputDeckMembers > {

  public:
    //! \brief Walker input deck keywords
    //! \see See also tk::grm::use
    using keywords = brigand::set< kw::precision
                                 , kw::end
                                 , kw::depvar
                                 , kw::title
                                 , kw::statistics
                                 , kw::interval
                                 , kw::pdfs
                                 , kw::filetype
                                 , kw::pdf_policy
                                 , kw::pdf_centering
                                 , kw::txt_float_format
                                 , kw::npar
                                 , kw::nstep
                                 , kw::term
                                 , kw::dt
                                 , kw::ttyi
                                 , kw::rngs
                                 , kw::ncomp
                                 , kw::rng
                                 , kw::walker
                                 , kw::init
                                 , kw::coeff
                                 , kw::diag_ou
                                 , kw::ornstein_uhlenbeck
                                 , kw::skewnormal
                                 , kw::gamma
                                 , kw::dirichlet
                                 , kw::mixdirichlet
                                 , kw::gendir
                                 , kw::wrightfisher
                                 , kw::beta
                                 , kw::sde_sigmasq
                                 , kw::sde_theta
                                 , kw::sde_mu
                                 , kw::sde_mean
                                 , kw::sde_cov
                                 , kw::mean_gradient
                                 , kw::sde_T
                                 , kw::sde_lambda
                                 , kw::sde_b
                                 , kw::sde_S
                                 , kw::sde_c
                                 , kw::sde_kappa
                                 , kw::sde_omega
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
                                 , kw::r123_threefry
                                 , kw::r123_philox
                                 , kw::const_shear
                                 , kw::stationary
                                 , kw::position
                                 , kw::velocity
                                 , kw::inst_velocity
                                 , kw::seed
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
                                 , kw::constcoeff
                                 , kw::decay
                                 , kw::raw
                                 , kw::zero
                                 , kw::elem
                                 , kw::node
                                 , kw::txt
                                 , kw::gmshtxt
                                 , kw::gmshbin
                                 , kw::exodusii
                                 , kw::overwrite
                                 , kw::multiple
                                 , kw::evolution
                                 , kw::txt_float_default
                                 , kw::txt_float_fixed
                                 , kw::txt_float_scientific
                                 , kw::numfracbeta
                                 , kw::sde_rho2
                                 , kw::sde_rho
                                 , kw::sde_rcomma
                                 , kw::icdelta
                                 , kw::spike
                                 , kw::sde_bprime
                                 , kw::sde_kappaprime
                                 , kw::mixnumfracbeta
                                 , kw::mixmassfracbeta
                                 , kw::massfracbeta
                                 , kw::sde_r
                                 , kw::homogeneous
                                 , kw::homdecay
                                 , kw::montecarlo_homdecay
                                 , kw::hydrotimescale
                                 , kw::jointbeta
                                 , kw::jointdelta
                                 , kw::jointgaussian
                                 , kw::jointcorrgaussian
                                 , kw::jointdirichlet
                                 , kw::icbeta
                                 , kw::betapdf
                                 , kw::dirichletpdf
                                 , kw::sde_c0
                                 , kw::icgaussian
                                 , kw::icjointgaussian
                                 , kw::icdirichlet
                                 , kw::gaussian
                                 , kw::dissipation
                                 , kw::jointgamma
                                 , kw::hydrotimescales
                                 , kw::hydroproductions
                                 , kw::eq_A005H
                                 , kw::eq_A005S
                                 , kw::eq_A005L
                                 , kw::eq_A05H
                                 , kw::eq_A05S
                                 , kw::eq_A05L
                                 , kw::eq_A075H
                                 , kw::eq_A075S
                                 , kw::eq_A075L
                                 , kw::prod_A005H
                                 , kw::prod_A005S
                                 , kw::prod_A005L
                                 , kw::prod_A05H
                                 , kw::prod_A05S
                                 , kw::prod_A05L
                                 , kw::prod_A075H
                                 , kw::prod_A075S
                                 , kw::prod_A075L
                                 , kw::gnorm
                                 , kw::gnorm_accurate
                                 , kw::gamma_method
                                 , kw::icgamma
                                 , kw::gammapdf
                                 , kw::sde_c3
                                 , kw::sde_c4
                                 , kw::sde_com1
                                 , kw::sde_com2
                                 , kw::fullvar
                                 , kw::fluctuation
                                 , kw::product
                                 , kw::solve
                                 , kw::variant
                                 , kw::slm
                                 , kw::glm
                                 , kw::normalization
                                 , kw::light
                                 , kw::heavy
                                 >;

    //! \brief Constructor: set all defaults
    //! \param[in] cl Previously parsed and store command line
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type.
    explicit InputDeck( const CmdLine& cl = {} ) {
      // Set previously parsed command line
      get< tag::cmd >() = cl;
      // Default discretization parameters
      get< tag::discr, tag::npar >() = 1;
      get< tag::discr, tag::nstep >() =
        std::numeric_limits< kw::nstep::info::expect::type >::max();
      get< tag::discr, tag::term >() = 1.0;
      get< tag::discr, tag::dt >() = 0.5;
      // Default txt floating-point output precision in digits
      get< tag::prec, tag::stat >() = std::cout.precision();
      get< tag::prec, tag::pdf >() = std::cout.precision();
      // Default intervals
      get< tag::interval, tag::tty >() = 1;
      get< tag::interval, tag::stat >() = 1;
      get< tag::interval, tag::pdf >() = 1;
      // Default requested statistics
      get< tag::stat >() = std::vector< tk::ctr::Product >();
      // Initialize help
      const auto& ctrinfoFill = tk::ctr::Info( get< tag::cmd, tag::ctrinfo >() );
      brigand::for_each< keywords >( ctrinfoFill );
    }

    //! Extract moment names of requested statistics
    std::vector< std::string > momentNames( std::function<
      bool ( const std::vector< tk::ctr::Term >& ) > momentType ) const
    {
      std::vector< std::string > names;
      for (const auto& product : get< tag::stat >()) {
        if (momentType( product )) {
          names.emplace_back( std::string() );
          for (const auto& term : product)
            names.back() += tk::ctr::Term( term.var, term.field );
        }
      }
      return names;
    }

    //! Query if there are any statistics or PDFs to estimate
    //! \return True if there are any statistics or PDFs to estimate
    bool stat()
    { return !get< tag::stat >().empty() || !get< tag::pdf >().empty(); }

    //! Query if there are any PDFs to estimate
    //! \return True if there are any PDFs to estimate
    bool pdf() { return !get< tag::pdf >().empty(); }

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
} // Walker::

#endif // WalkerInputDeck_h
