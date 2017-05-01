// *****************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/InputDeck.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Walker's input deck
  \details   Walker's input deck
*/
// *****************************************************************************
#ifndef WalkerInputDeck_h
#define WalkerInputDeck_h

#include <limits>
#include <iostream>

#include "NoWarning/set.h"
#include "NoWarning/for_each.h"

#include "QuinoaConfig.h"
#include "Control.h"
#include "HelpFactory.h"
#include "Walker/CmdLine/CmdLine.h"
#include "Walker/Components.h"

namespace walker {
namespace ctr {

//! InputDeck : Control< specialized to Walker >, see Types.h
class InputDeck :
  public tk::Control< // tag           type
                      tag::title,      kw::title::info::expect::type,
                      tag::selected,   selects,
                      tag::discr,      discretization,
                      tag::prec,       precision,
                      tag::flformat,   floatformat,
                      tag::component,  ncomps,
                      tag::interval,   intervals,
                      tag::cmd,        CmdLine,
                      tag::param,      parameters,
                      tag::stat,       std::vector< tk::ctr::Product >,
                      tag::pdf,        std::vector< tk::ctr::Probability >,
                      tag::error,      std::vector< std::string > > {

  private:
    //! Function object for extracting dependent variable vectors from
    //! components
    struct depvar {
      //! Need to store reference to host class whose data we operate on
      const InputDeck* const m_host;
      //! Internal reference to the vector of vectors that will hold the
      //! extracted dependent variables for all equations
      std::vector< std::vector< char > >& m_vars;
      //! Constructor: store host object pointer and the reference to the
      //! container storing the dependent variables letters.
      depvar( const InputDeck* const host,
              std::vector< std::vector< char > >& vars ) :
        m_host( host ), m_vars( vars ) {}
      //! Functional call operator templated on the type that adds all dependent
      //! variables (stored in a vector) for type U. (There can be multiple
      //! systems of the same equation type, differentiated by a different
      //! dependen:t variable.)
      template< typename U > void operator()( U ) {
        m_vars.push_back( m_host->get< tag::param, U, tag::depvar >() );
      }
    };

  public:
    //! \brief Walker input deck keywords
    //! \details Since there are more than 20 and boost::mpl only allows maxium
    //!   20 items in a set by default (and I don't want to mess with
    //!   preprocessor-generated boost::mpl headers), the whole set is broken up
    //!   into several sets each containing 20 keywords.
    //! \see See also tk::grm::use
    //! \details If you add a new mpl::set, make sure you also add it to
    //!   waker::deck::use.
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
                                     , kw::ncomp
                                     , kw::rng
                                     , kw::walker
                                     >;
    using keywords2 = boost::mpl::set< kw::init
                                     , kw::coeff
                                     , kw::diag_ou
                                     , kw::ornstein_uhlenbeck
                                     , kw::skewnormal
                                     , kw::gamma
                                     , kw::dirichlet
                                     , kw::gendir
                                     , kw::wrightfisher
                                     , kw::beta
                                     , kw::sde_sigmasq
                                     , kw::sde_theta
                                     , kw::sde_mu
                                     , kw::sde_T
                                     , kw::sde_lambda
                                     , kw::sde_b
                                     , kw::sde_S
                                     , kw::sde_c
                                     , kw::sde_kappa
                                     , kw::sde_omega
                                     >;
    using keywords3 = boost::mpl::set< kw::cja
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
                                     , kw::beta_method
                                     , kw::standard
                                     , kw::accurate
                                     , kw::boxmuller
                                     , kw::boxmuller2
                                     , kw::icdf
                                     #endif
                                     >;
    using keywords5 = boost::mpl::set< kw::constant
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
                                     , kw::sde_rcomma
                                     , kw::icdelta
                                     >;
    using keywords6 = boost::mpl::set< kw::spike
                                     , kw::sde_bprime
                                     , kw::sde_kappaprime
                                     , kw::mixnumfracbeta
                                     , kw::mixmassfracbeta
                                     , kw::massfracbeta
                                     , kw::sde_r
                                     , kw::homdecay
                                     , kw::montecarlo_homdecay
                                     , kw::hydrotimescale
                                     , kw::jointbeta
                                     , kw::jointdelta
                                     , kw::jointbeta
                                     , kw::icbeta
                                     , kw::betapdf
                                     >;
    using keywords7 = boost::mpl::set< kw::hydrotimescales
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
                                     >;

    //! \brief Constructor: set all defaults
    //! \param[in] cl Previously parsed and store command line
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type.
    //! \author J. Bakosi
    explicit InputDeck( const CmdLine& cl = {} ) {
      // Set previously parsed command line
      set< tag::cmd >( cl );
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
      set< tag::interval, tag::stat >( 1 );
      set< tag::interval, tag::pdf >( 1 );
      // Default requested statistics
      set< tag::stat >( std::vector< tk::ctr::Product >() );
      // Initialize help
      const auto& ctrinfoFill = tk::ctr::Info( get< tag::cmd, tag::ctrinfo >() );
      boost::mpl::for_each< keywords1 >( ctrinfoFill );
      boost::mpl::for_each< keywords2 >( ctrinfoFill );
      boost::mpl::for_each< keywords3 >( ctrinfoFill );
      boost::mpl::for_each< keywords4 >( ctrinfoFill );
      boost::mpl::for_each< keywords5 >( ctrinfoFill );
      boost::mpl::for_each< keywords6 >( ctrinfoFill );
      boost::mpl::for_each< keywords7 >( ctrinfoFill );
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

    // Count number requested PDFs with given number of sample space dimensions
    template< std::size_t d >
    std::size_t npdf() const {
      std::size_t n = 0;
      for (const auto& bs : get< tag::discr, tag::binsize >())
        if (bs.size() == d) ++n;
      return n;
    }

    //! Extract vector of vector of dependent variables from components
    std::vector< std::vector< char > > depvars() const {
      std::vector< std::vector< char > > vars;
      boost::mpl::for_each< ncomps::tags >( depvar( this, vars ) );
      return vars;
    }

    //! Query if there are any statistics or PDFs to estimate
    //! \return True if there are any statistics or PDFs to estimate
    bool stat()
    { return !get< tag::stat >().empty() || !get< tag::pdf >().empty(); }

    //! Query if there are any PDFs to estimate
    //! \return True if there are any PDFs to estimate
    bool pdf() { return !get< tag::pdf >().empty(); }

    //! Pack/Unpack
    void pup( PUP::er& p ) {
      tk::Control< tag::title,      kw::title::info::expect::type,
                   tag::selected,   selects,
                   tag::discr,      discretization,
                   tag::prec,       precision,
                   tag::flformat,   floatformat,
                   tag::component,  ncomps,
                   tag::interval,   intervals,
                   tag::cmd,        CmdLine,
                   tag::param,      parameters,
                   tag::stat,       std::vector< tk::ctr::Product >,
                   tag::pdf,        std::vector< tk::ctr::Probability >,
                   tag::error,      std::vector< std::string > >::pup(p);
    }
    friend void operator|( PUP::er& p, InputDeck& c ) { c.pup(p); }
};

} // ctr::
} // Walker::

#endif // WalkerInputDeck_h
