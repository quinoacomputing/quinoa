// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/InputDeck.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Inciter's input deck definition
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the control file parsing of the
     computational shock hydrodynamics tool, Inciter.
*/
// *****************************************************************************
#ifndef InciterInputDeck_h
#define InciterInputDeck_h

#include <limits>
#include <iomanip>
#include <iostream>

#include "NoWarning/set.h"
#include "NoWarning/for_each.h"

#include "Control.h"
#include "Inciter/CmdLine/CmdLine.h"
#include "Inciter/Components.h"

namespace inciter {
namespace ctr {

//! \brief InputDeck : Control< specialized to Inciter >, see Types.h,
//! \details The stack is a tagged tuple, a hierarchical heterogeneous data
//!    structure where all parsed information is stored.
//! \see Base/TaggedTuple.h
//! \see Control/Inciter/Types.h
//! \author J. Bakosi
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
                      tag::error,      std::vector< std::string > > {

  public:
    //! \brief Inciter input deck keywords
    //! \details Since there are more than 20 and boost::mpl only allows maxium
    //!   20 items in a set by default (and I don't want to mess with
    //!   preprocessor-generated boost::mpl headers), the whole set is broken up
    //!   into several sets each containing 20 keywords.
    //! \see tk::grm::use and its documentation
    using keywords1 = boost::mpl::set< kw::title,
                                       kw::nstep,
                                       kw::term,
                                       kw::t0,
                                       kw::dt,
                                       kw::ttyi,
                                       kw::transport,
                                       kw::end,
                                       kw::shear_diff,
                                       kw::slot_cyl,
                                       kw::problem,
                                       kw::plotvar,
                                       kw::interval,
                                       kw::partitioning,
                                       kw::algorithm,
                                       kw::rcb,
                                       kw::rib,
                                       kw::hsfc,
                                       kw::phg,
                                       kw::inciter >;
    using keywords2 = boost::mpl::set< kw::ncomp,
                                       kw::pde_diffusivity,
                                       kw::pde_lambda,
                                       kw::pde_u0,
                                       kw::poisson,
                                       kw::dir_neu,
                                       kw::bc_dirichlet,
                                       kw::sideset,
                                       kw::compflow,
                                       kw::ic,
                                       kw::txt_float_format,
                                       kw::txt_float_default,
                                       kw::txt_float_fixed,
                                       kw::txt_float_scientific,
                                       kw::precision,
                                       kw::diagnostics,
                                       kw::material,
                                       kw::id,
                                       kw::mat_gamma,
                                       kw::mat_mu >;
    using keywords3 = boost::mpl::set< kw::mat_cv,
                                       kw::mat_k,
                                       kw::npar,
                                       kw::physics,
                                       kw::advection,
                                       kw::advdiff,
                                       kw::compflow_navierstokes,
                                       kw::compflow_euler,
                                       kw::user_defined,
                                       kw::vortical_flow,
                                       kw::pde_alpha,
                                       kw::pde_beta,
                                       kw::pde_p0,
                                       kw::laplace,
                                       kw::ctau,
                                       kw::cfl,
                                       kw::mj,
                                       kw::depvar >;
    using keywords4 = boost::mpl::set< kw::amr,
                                       kw::amr_initial,
                                       kw::amr_uniform >;
                                     
    //! \brief Constructor: set defaults
    //! \param[in] cl Previously parsed and store command line
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type.
    //! \author J. Bakosi
    InputDeck( const CmdLine& cl = {} ) {
      // Set previously parsed command line
      set< tag::cmd >( cl );
      // Default discretization parameters
      set< tag::discr, tag::nstep >
         ( std::numeric_limits< kw::nstep::info::expect::type >::max() );
      set< tag::discr, tag::term >
         ( std::numeric_limits< kw::term::info::expect::type >::max() );
      set< tag::discr, tag::t0 >( 0.0 );
      set< tag::discr, tag::dt >( 0.0 );
      set< tag::discr, tag::cfl >( 0.0 );
      set< tag::discr, tag::ctau >( 1.0 );
      // Default AMR settings
      set< tag::selected, tag::initialamr >( tk::ctr::InitialAMRType::NONE );
      // Default txt floating-point output precision in digits
      set< tag::prec, tag::diag >( std::cout.precision() );
      // Default intervals
      set< tag::interval, tag::tty >( 1 );
      set< tag::interval, tag::field >( 1 );
      set< tag::interval, tag::diag >( 1 );
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
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \author J. Bakosi
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
} // inciter::

#endif // InciterInputDeck_h
