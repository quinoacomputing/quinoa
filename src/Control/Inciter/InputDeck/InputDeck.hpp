// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/InputDeck.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include "NoWarning/set.hpp"

#include "Inciter/CmdLine/CmdLine.hpp"
#include "Inciter/Components.hpp"

namespace inciter {

namespace ctr {

//! Member data for tagged tuple
using InputDeckMembers = brigand::list<
    tag::cmd,             CmdLine
  , tag::title,           kw::title::info::expect::type
  , tag::selected,        selects
  , tag::amr,             amr
  , tag::ale,             ale
  , tag::pref,            pref
  , tag::discr,           discretization
  , tag::prec,            precision
  , tag::flformat,        floatformat
  , tag::component,       ncomps
  , tag::sys,             std::map< tk::ctr::ncomp_t, tk::ctr::ncomp_t >
  , tag::output,          output_parameters
  , tag::param,           parameters
  , tag::couple,          couple
  , tag::diag,            diagnostics
  , tag::error,           std::vector< std::string >
  , tag::history,         history
>;

//! \brief InputDeck : Control< specialized to Inciter >, see Types.h,
//! \details The stack is a tagged tuple, a hierarchical heterogeneous data
//!    structure where all parsed information is stored.
//! \see Base/TaggedTuple.h
//! \see Control/Inciter/Types.h
class InputDeck : public tk::TaggedTuple< InputDeckMembers > {

  public:
    //! \brief Inciter input deck keywords
    //! \see tk::grm::use and its documentation
    using keywords = brigand::set< kw::title
                                 , kw::nstep
                                 , kw::term
                                 , kw::t0
                                 , kw::dt
                                 , kw::ttyi
                                 , kw::transport
                                 , kw::end
                                 , kw::shear_diff
                                 , kw::slot_cyl
                                 , kw::problem
                                 , kw::field_output
                                 , kw::refined
                                 , kw::interval_iter
                                 , kw::interval_time
                                 , kw::time_range
                                 , kw::partitioning
                                 , kw::algorithm
                                 , kw::rcb
                                 , kw::rib
                                 , kw::hsfc
                                 , kw::phg
                                 , kw::inciter
                                 , kw::ncomp
                                 , kw::nmat
                                 , kw::pde_diffusivity
                                 , kw::pde_lambda
                                 , kw::pde_u0
                                 , kw::bc_dirichlet
                                 , kw::sideset
                                 , kw::compflow
                                 , kw::multimat
                                 , kw::ic
                                 , kw::box
                                 , kw::meshblock
                                 , kw::lua
                                 , kw::materialid
                                 , kw::blockid
                                 , kw::volume
                                 , kw::mass
                                 , kw::density
                                 , kw::velocity
                                 , kw::position
                                 , kw::acceleration
                                 , kw::fntype
                                 , kw::fn
                                 , kw::move
                                 , kw::initiate
                                 , kw::impulse
                                 , kw::linear
                                 , kw::pressure
                                 , kw::energy
                                 , kw::energy_content
                                 , kw::temperature
                                 , kw::xmin
                                 , kw::xmax
                                 , kw::ymin
                                 , kw::ymax
                                 , kw::zmin
                                 , kw::zmax
                                 , kw::txt_float_format
                                 , kw::txt_float_default
                                 , kw::txt_float_fixed
                                 , kw::txt_float_scientific
                                 , kw::precision
                                 , kw::diagnostics
                                 , kw::history_output
                                 , kw::mesh
                                 , kw::filename
                                 , kw::location
                                 , kw::orientation
                                 , kw::reference
                                 , kw::couple
                                 , kw::material
                                 , kw::id
                                 , kw::eos
                                 , kw::stiffenedgas
                                 , kw::jwl
                                 , kw::mat_gamma
                                 , kw::mat_pstiff
                                 , kw::w_gru
                                 , kw::A_jwl
                                 , kw::B_jwl
                                 , kw::C_jwl
                                 , kw::R1_jwl
                                 , kw::R2_jwl
                                 , kw::rho0_jwl
                                 , kw::de_jwl
                                 , kw::rhor_jwl
                                 , kw::Pr_jwl
                                 , kw::mat_mu
                                 , kw::mat_cv
                                 , kw::mat_k
                                 , kw::physics
                                 , kw::advection
                                 , kw::advdiff
                                 , kw::navierstokes
                                 , kw::euler
                                 , kw::energy_pill
                                 , kw::user_defined
                                 , kw::vortical_flow
                                 , kw::pde_alpha
                                 , kw::pde_beta
                                 , kw::pde_p0
                                 , kw::ctau
                                 , kw::cfl
                                 , kw::dvcfl
                                 , kw::vortmult
                                 , kw::mj
                                 , kw::elem
                                 , kw::node
                                 , kw::depvar
                                 , kw::outvar
                                 , kw::outvar_density
                                 , kw::outvar_xmomentum
                                 , kw::outvar_ymomentum
                                 , kw::outvar_zmomentum
                                 , kw::outvar_specific_total_energy
                                 , kw::outvar_volumetric_total_energy
                                 , kw::outvar_xvelocity
                                 , kw::outvar_yvelocity
                                 , kw::outvar_zvelocity
                                 , kw::outvar_pressure
                                 , kw::outvar_material_indicator
                                 , kw::outvar_analytic
                                 , kw::nl_energy_growth
                                 , kw::pde_betax
                                 , kw::pde_betay
                                 , kw::pde_betaz
                                 , kw::pde_ce
                                 , kw::pde_kappa
                                 , kw::pde_r0
                                 , kw::rayleigh_taylor
                                 , kw::taylor_green
                                 , kw::filetype
                                 , kw::exodusii
                                 , kw::root
                                 , kw::error
                                 , kw::l2
                                 , kw::linf
                                 , kw::fct
                                 , kw::fctclip
                                 , kw::fcteps
                                 , kw::sysfct
                                 , kw::sysfctvar
                                 , kw::pelocal_reorder
                                 , kw::operator_reorder
                                 , kw::steady_state
                                 , kw::residual
                                 , kw::rescomp
                                 , kw::amr
                                 , kw::ale
                                 , kw::smoother
                                 , kw::laplace
                                 , kw::helmholtz
                                 , kw::meshvelocity
                                 , kw::meshvel_maxit
                                 , kw::meshvel_tolerance
                                 , kw::mesh_motion
                                 , kw::meshforce
                                 , kw::none
                                 , kw::sine
                                 , kw::fluid
                                 , kw::amr_t0ref
                                 , kw::amr_dtref
                                 , kw::amr_dtref_uniform
                                 , kw::amr_dtfreq
                                 , kw::amr_maxlevels
                                 , kw::amr_initial
                                 , kw::amr_uniform
                                 , kw::amr_uniform_derefine
                                 , kw::amr_initial_conditions
                                 , kw::amr_coords
                                 , kw::amr_error
                                 , kw::amr_jump
                                 , kw::amr_hessian
                                 , kw::amr_refvar
                                 , kw::amr_tolref
                                 , kw::amr_tolderef
                                 , kw::amr_edgelist
                                 , kw::amr_xminus
                                 , kw::amr_xplus
                                 , kw::amr_yminus
                                 , kw::amr_yplus
                                 , kw::amr_zminus
                                 , kw::amr_zplus
                                 , kw::pref
                                 , kw::pref_indicator
                                 , kw::pref_spectral_decay
                                 , kw::pref_non_conformity
                                 , kw::pref_ndofmax
                                 , kw::pref_tolref
                                 , kw::shock_detector_coeff
                                 , kw::scheme
                                 , kw::diagcg
                                 , kw::alecg
                                 , kw::dg
                                 , kw::p0p1
                                 , kw::dgp1
                                 , kw::dgp2
                                 , kw::pdg
                                 , kw::fv
                                 , kw::flux
                                 , kw::laxfriedrichs
                                 , kw::hllc
                                 , kw::upwind
                                 , kw::ausm
                                 , kw::hll
                                 , kw::limiter
                                 , kw::cweight
                                 , kw::nolimiter
                                 , kw::wenop1
                                 , kw::superbeep1
                                 , kw::vertexbasedp1
                                 , kw::accuracy_test
                                 , kw::limsol_projection
                                 , kw::prelax
                                 , kw::prelax_timescale
                                 , kw::intsharp
                                 , kw::intsharp_param
                                 , kw::bc_sym
                                 , kw::bc_inlet
                                 , kw::bc_outlet
                                 , kw::bc_farfield
                                 , kw::bc_extrapolate
                                 , kw::bc_timedep
                                 , kw::bc_stag
                                 , kw::bc_skip
                                 , kw::sponge
                                 , kw::point
                                 , kw::init_time
                                 , kw::front_width
                                 , kw::radius
                                 , kw::gauss_hump
                                 , kw::rotated_sod_shocktube
                                 , kw::cyl_advect
                                 , kw::cyl_vortex
                                 , kw::shedding_flow
                                 , kw::sod_shocktube
                                 , kw::sedov_blastwave
                                 , kw::interface_advection
                                 , kw::gauss_hump_compflow
                                 , kw::waterair_shocktube
                                 , kw::shock_hebubble
                                 , kw::underwater_ex
                                 , kw::shockdensity_wave
                                 , kw::equilinterface_advect
                                 , kw::sinewave_packet
                                 , kw::richtmyer_meshkov >;

    //! Set of tags to ignore when printing this InputDeck
    using ignore = CmdLine::ignore;

    //! \brief Constructor: set defaults
    //! \param[in] cl Previously parsed and store command line
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type.
    explicit InputDeck( const CmdLine& cl = {} ) {
      // Set previously parsed command line
      get< tag::cmd >() = cl;
      // Default discretization parameters
      get< tag::discr, tag::nstep >() =
         std::numeric_limits< kw::nstep::info::expect::type >::max();
      get< tag::discr, tag::term >() =
         std::numeric_limits< kw::term::info::expect::type >::max();
      get< tag::discr, tag::t0 >() = 0.0;
      get< tag::discr, tag::dt >() = 0.0;
      get< tag::discr, tag::cfl >() = 0.0;
      get< tag::discr, tag::fct >() = true;
      get< tag::discr, tag::fctclip >() = false;
      get< tag::discr, tag::ctau >() = 1.0;
      get< tag::discr, tag::fcteps >() =
        std::numeric_limits< tk::real >::epsilon();
      get< tag::discr, tag::pelocal_reorder >() = false;
      get< tag::discr, tag::operator_reorder >() = false;
      get< tag::discr, tag::steady_state >() = false;
      get< tag::discr, tag::residual >() = 1.0e-8;
      get< tag::discr, tag::rescomp >() = 1;
      get< tag::discr, tag::scheme >() = SchemeType::DiagCG;
      get< tag::discr, tag::ndof >() = 1;
      get< tag::discr, tag::limiter >() = LimiterType::NOLIMITER;
      get< tag::discr, tag::cweight >() = 1.0;
      get< tag::discr, tag::ndof >() = 1;
      get< tag::discr, tag::rdof >() = 1;
      get< tag::discr, tag::accuracy_test >() = false;
      get< tag::discr, tag::limsol_projection >() = true;
      get< tag::discr, tag::shock_detector_coeff >() = 1.0;
      // Default field output file type
      get< tag::selected, tag::filetype >() = tk::ctr::FieldFileType::EXODUSII;
      // Default AMR settings
      get< tag::amr, tag::amr >() = false;
      get< tag::amr, tag::t0ref >() = false;
      get< tag::amr, tag::dtref >() = false;
      get< tag::amr, tag::dtref_uniform >() = false;
      get< tag::amr, tag::dtfreq >() = 3;
      get< tag::amr, tag::maxlevels >() = 2;
      get< tag::amr, tag::error >() = AMRErrorType::JUMP;
      get< tag::amr, tag::tolref >() = 0.2;
      get< tag::amr, tag::tolderef >() = 0.05;
      auto rmax =
        std::numeric_limits< kw::amr_xminus::info::expect::type >::max() / 100;
      get< tag::amr, tag::xminus >() = rmax;
      get< tag::amr, tag::xplus >() = -rmax;
      get< tag::amr, tag::yminus >() = rmax;
      get< tag::amr, tag::yplus >() = -rmax;
      get< tag::amr, tag::zminus >() = rmax;
      get< tag::amr, tag::zplus >() = -rmax;
      // Default ALE settings
      get< tag::ale, tag::ale >() = false;
      get< tag::ale, tag::smoother >() = MeshVelocitySmootherType::NONE;
      get< tag::ale, tag::dvcfl >() = 0.0;
      get< tag::ale, tag::vortmult >() = 0.0;
      get< tag::ale, tag::maxit >() = 5;
      get< tag::ale, tag::tolerance >() = 1.0e-2;
      // Default p-refinement settings
      get< tag::pref, tag::pref >() = false;
      get< tag::pref, tag::indicator >() = PrefIndicatorType::SPECTRAL_DECAY;
      get< tag::pref, tag::ndofmax >() = 10;
      get< tag::pref, tag::tolref >() = 0.5;
      // Default txt floating-point output precision in digits
      get< tag::prec, tag::diag >() = std::cout.precision();
      get< tag::prec, tag::history >() = std::cout.precision();
      // Default intervals
      get< tag::output, tag::iter, tag::tty >() = 1;
      get< tag::output, tag::iter, tag::diag >() = 1;
      get< tag::output, tag::iter, tag::field >() =
        std::numeric_limits< kw::interval_iter::info::expect::type >::max();
      get< tag::output, tag::iter, tag::history >() =
        std::numeric_limits< kw::interval_iter::info::expect::type >::max();
      // Initialize help: fill own keywords
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

    //! Extract surface side set ids along which user wants to save solution
    //! \return Unique set of surface side set ids along which user wants to
    //!   save solution field variables
    //! \note This returns an ordered set so the order of the set ids are
    //!   always the same.
    std::set< int > outsets() const {
      std::set< int > ids;
      for (const auto& s : get< tag::cmd, tag::io, tag::surface >()) {
        std::stringstream conv( s );
        int num;
        conv >> num;
        ids.insert( num );
      }
      return ids;
    }

    //! Extract field output variable names the user wants to save
    //! \param[in] c Extract variables only with this centering
    //! \return Unique set of field output variable names user wants
    //! \note This returns an ordered set so the order of the variable names
    //!   are alphabetical and unique.
    std::set< OutVar > outvars( tk::Centering c ) const {
      std::set< OutVar > vars;
      for (const auto& v : get< tag::cmd, tag::io, tag::outvar >()) {
        if (v.centering == c) vars.insert( v );
      }
      return vars;
    }

    //! Extract field output variable names and aliases the user configured
    //! \return Map of field output variable names and alias for all
    //!    output variables the user configured
    std::map< std::string, std::string > outvar_aliases() const {
      std::map< std::string, std::string > aliases;
      for (const auto& v : get< tag::cmd, tag::io, tag::outvar >())
         if (!v.alias.empty()) {
           std::stringstream s;
           s << v;
           aliases[ s.str() ] = v.alias;
         }
      return aliases;
    }

    //! Extract list of mesh filenames (each assigned to a solver)
    std::vector< std::string > mesh() const {
      using PDETypes = parameters::Keys;
      std::vector< std::string > meshes;
      brigand::for_each< PDETypes >( Meshes( *this, meshes ) );
      return meshes;
    }

    //! Extract list of dependent variables (each configuring a solver)
    std::vector< char > depvar() const {
      using PDETypes = parameters::Keys;
      std::vector< char > depvar;
      brigand::for_each< PDETypes >( Depvar( *this, depvar ) );
      return depvar;
    }

    //! Query special point BC configuration
    //! \tparam eq PDE type to query
    //! \tparam sbc Special BC type to query, e.g., stagnation, skip
    //! \param[in] system Equation system id
    //! \return Vectors configuring the special points and their radii
    template< class eq, class sbc >
    std::tuple< std::vector< tk::real >, std::vector< tk::real > >
    specialBC( std::size_t system ) {
      const auto& bcspec = get< tag::param, eq, sbc >();
      const auto& point = bcspec.template get< tag::point >();
      const auto& radius = bcspec.template get< tag::radius >();
      std::vector< tk::real > pnt;
      std::vector< tk::real > rad;
      if (point.size() > system && radius.size() > system) {
        pnt = point[ system ];
        rad = radius[ system ];
      }
      Assert( pnt.size() == 3*rad.size(), "Size mismatch" );
      return { std::move(pnt), std::move(rad) };
    }

    //! Query scheme centering
    //! \return Scheme centering
    tk::Centering centering() const
    { return ctr::Scheme().centering( get< tag::discr, tag::scheme >() ); }

  private:
    //! Function object to extract the mesh filenames assigned to solvers
    //! \details This is instantiated for all PDE types at compile time. It goes
    //!   through all configured solvers (equation system configuration blocks)
    //!   and builds a list of all mesh filenames associated to all solvers in
    //!   the input file.
    struct Meshes {
      const InputDeck& inputdeck;
      std::vector< std::string >& filenames;
      explicit Meshes( const InputDeck& i, std::vector< std::string >& f )
        : inputdeck(i), filenames(f) {}
      template< typename eq > void operator()( brigand::type_<eq> ) {
        const auto& eq_mesh_filename =
           inputdeck.get< tag::param, eq, tag::mesh, tag::filename >();
        for (const auto& f : eq_mesh_filename) filenames.push_back( f );
      }
    };

    //! Function object to extract the dependent variables assigned to solvers
    //! \details This is instantiated for all PDE types at compile time. It goes
    //!   through all configured solvers (equation system configuration blocks)
    //!   and builds a list of all dependent variables associated to all solvers
    //!   in the input file.
    struct Depvar {
      const InputDeck& inputdeck;
      std::vector< char >& depvar;
      explicit Depvar( const InputDeck& i, std::vector< char >& d ) :
        inputdeck(i), depvar(d) {}
      template< typename eq > void operator()( brigand::type_<eq> ) {
        const auto& eq_depvar = inputdeck.get< tag::param, eq, tag::depvar >();
        for (const auto& d : eq_depvar) depvar.push_back( d );
      }
    };
};

} // ctr::
} // inciter::

#endif // InciterInputDeck_h
