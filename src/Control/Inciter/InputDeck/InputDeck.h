// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/InputDeck.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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

#include <brigand/algorithms/for_each.hpp>

#include "NoWarning/set.h"

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
class InputDeck :
  public tk::Control< // tag           type
                      tag::title,      kw::title::info::expect::type,
                      tag::selected,   selects,
                      tag::amr,        amr,
                      tag::discr,      discretization,
                      tag::prec,       precision,
                      tag::flformat,   floatformat,
                      tag::component,  ncomps,
                      tag::interval,   intervals,
                      tag::cmd,        CmdLine,
                      tag::param,      parameters,
                      tag::diag,       diagnostics,
                      tag::error,      std::vector< std::string > > {

  public:
    //! \brief Inciter input deck keywords
    //! \see tk::grm::use and its documentation
    using keywords = brigand::set< kw::title,
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
                                   kw::inciter,
                                   kw::ncomp,
                                   kw::pde_diffusivity,
                                   kw::pde_lambda,
                                   kw::pde_u0,
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
                                   kw::mat_mu,
                                   kw::mat_cv,
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
                                   kw::ctau,
                                   kw::cfl,
                                   kw::mj,
                                   kw::depvar,
                                   kw::nl_energy_growth,
                                   kw::pde_betax,
                                   kw::pde_betay,
                                   kw::pde_betaz,
                                   kw::pde_ce,
                                   kw::pde_kappa,
                                   kw::pde_r0,
                                   kw::rayleigh_taylor,
                                   kw::taylor_green,
                                   kw::filetype,
                                   kw::exodusii,
                                   kw::root,
                                   kw::error,
                                   kw::l2,
                                   kw::linf,
                                   kw::fct,
                                   kw::reorder,
                                   kw::amr,
                                   kw::amr_t0ref,
                                   kw::amr_dtref,
                                   kw::amr_dtfreq,
                                   kw::amr_initial,
                                   kw::amr_uniform,
                                   kw::amr_initial_conditions,
                                   kw::amr_coords,
                                   kw::amr_error,
                                   kw::amr_jump,
                                   kw::amr_hessian,
                                   kw::amr_refvar,
                                   kw::amr_initref,
                                   kw::amr_coordref,
                                   kw::amr_xminus,
                                   kw::amr_xplus,
                                   kw::amr_yminus,
                                   kw::amr_yplus,
                                   kw::amr_zminus,
                                   kw::amr_zplus,
                                   kw::scheme,
                                   kw::matcg,
                                   kw::diagcg,
                                   kw::dg,
                                   kw::dgp1,
                                   kw::flux,
                                   kw::laxfriedrichs,
                                   kw::hllc,
                                   kw::limiter,
                                   kw::cweight,
                                   kw::nolimiter,
                                   kw::wenop1,
                                   kw::bc_sym,
                                   kw::bc_inlet,
                                   kw::bc_outlet,
                                   kw::bc_extrapolate,
                                   kw::gauss_hump,
                                   kw::cyl_advect,
                                   kw::sod_shocktube >;

    //! \brief Constructor: set defaults
    //! \param[in] cl Previously parsed and store command line
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type.
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
      set< tag::discr, tag::fct >( true );
      set< tag::discr, tag::reorder >( false );
      set< tag::discr, tag::ctau >( 1.0 );
      set< tag::discr, tag::scheme >( SchemeType::MatCG );
      set< tag::discr, tag::flux >( FluxType::HLLC );
      set< tag::discr, tag::ndof >( 1 );
      set< tag::discr, tag::limiter >( LimiterType::NOLIMITER );
      set< tag::discr, tag::cweight >( 1.0 );
      // Default field output file type
      set< tag::selected, tag::filetype >( tk::ctr::FieldFileType::EXODUSII );
      // Default AMR settings
      set< tag::amr, tag::amr >( false );
      set< tag::amr, tag::t0ref >( false );
      set< tag::amr, tag::dtref >( false );
      set< tag::amr, tag::dtfreq >( 3 );
      set< tag::amr, tag::error >( AMRErrorType::JUMP );
      auto rmax =
        std::numeric_limits< kw::amr_xminus::info::expect::type >::max();
      set< tag::amr, tag::xminus >( rmax );
      set< tag::amr, tag::xplus >( rmax );
      set< tag::amr, tag::yminus >( rmax );
      set< tag::amr, tag::yplus >( rmax );
      set< tag::amr, tag::zminus >( rmax );
      set< tag::amr, tag::zplus >( rmax );
      // Default txt floating-point output precision in digits
      set< tag::prec, tag::diag >( std::cout.precision() );
      // Default intervals
      set< tag::interval, tag::tty >( 1 );
      set< tag::interval, tag::field >( 1 );
      set< tag::interval, tag::diag >( 1 );
      // Initialize help: fill own keywords
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
                   tag::amr,        amr,
                   tag::discr,      discretization,
                   tag::prec,       precision,
                   tag::flformat,   floatformat,
                   tag::component,  ncomps,
                   tag::interval,   intervals,
                   tag::cmd,        CmdLine,
                   tag::param,      parameters,
                   tag::diag,       diagnostics,
                   tag::error,      std::vector< std::string > >::pup(p);
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i InputDeck object reference
    friend void operator|( PUP::er& p, InputDeck& i ) { i.pup(p); }
    //@}
};

} // ctr::
} // inciter::

#endif // InciterInputDeck_h
