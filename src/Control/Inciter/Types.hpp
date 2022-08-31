// *****************************************************************************
/*!
  \file      src/Control/Inciter/Types.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Types for Incitier's parsers
  \details   Types for Incitier's parsers. This file defines the components of
    the agged tuple that stores heteroegeneous objects in a hierarchical way.
    These components are therefore part of the grammar stack that is filled
    during parsing (both command-line argument parsing and control file
    parsing).
*/
// *****************************************************************************
#ifndef IncitierTypes_h
#define IncitierTypes_h

#include "Tags.hpp"
#include "Base/Types.hpp"
#include "Inciter/Options/PDE.hpp"
#include "Inciter/Options/Problem.hpp"
#include "Inciter/Options/Scheme.hpp"
#include "Inciter/Options/Limiter.hpp"
#include "Inciter/Options/Flux.hpp"
#include "Inciter/Options/Initiate.hpp"
#include "Inciter/Options/AMRInitial.hpp"
#include "Inciter/Options/AMRError.hpp"
#include "Inciter/Options/PrefIndicator.hpp"
#include "Inciter/Options/MeshVelocity.hpp"
#include "Inciter/Options/MeshVelocitySmoother.hpp"
#include "Inciter/Options/Material.hpp"
#include "Options/PartitioningAlgorithm.hpp"
#include "Options/TxtFloatFormat.hpp"
#include "Options/FieldFile.hpp"
#include "Options/Error.hpp"
#include "Options/UserTable.hpp"
#include "PUPUtil.hpp"
#include "OutVar.hpp"
#include "Transfer.hpp"

namespace inciter {
namespace ctr {

using namespace tao;

//! Storage of selected options
using selects = tk::TaggedTuple< brigand::list<
    tag::pde,         std::vector< ctr::PDEType >        //!< Partial diff eqs
  , tag::partitioner, tk::ctr::PartitioningAlgorithmType //!< Mesh partitioner
  , tag::filetype,    tk::ctr::FieldFileType       //!< Field output file type
> >;

//! Adaptive-mesh refinement options
using amr = tk::TaggedTuple< brigand::list<
    tag::amr,     bool                            //!< AMR on/off
  , tag::t0ref,   bool                            //!< AMR before t<0 on/off
  , tag::dtref,   bool                            //!< AMR during t>0 on/off
  , tag::dtref_uniform, bool                      //!< Force dtref uniform-only
  , tag::dtfreq,  kw::amr_dtfreq::info::expect::type //!< Refinement frequency
  , tag::maxlevels, kw::amr_maxlevels::info::expect::type //!< Max refine levels
  , tag::init,    std::vector< AMRInitialType >   //!< List of initial AMR types
  , tag::refvar,  std::vector< std::string >      //!< List of refinement vars
  , tag::id,      std::vector< std::size_t >      //!< List of refvar indices
  , tag::error,   AMRErrorType                    //!< Error estimator for AMR
  , tag::tolref,  tk::real                        //!< Refine tolerance
  , tag::tolderef, tk::real                       //!< De-refine tolerance
  //! List of edges-node pairs
  , tag::edge,    std::vector< kw::amr_edgelist::info::expect::type >
  //! Refinement tagging edges with end-point coordinates lower than x coord
  , tag::xminus,  kw::amr_xminus::info::expect::type
  //! Refinement tagging edges with end-point coordinates higher than x coord
  , tag::xplus,  kw::amr_xplus::info::expect::type
  //! Refinement tagging edges with end-point coordinates lower than y coord
  , tag::yminus,  kw::amr_yminus::info::expect::type
  //! Refinement tagging edges with end-point coordinates higher than y coord
  , tag::yplus,  kw::amr_yplus::info::expect::type
  //! Refinement tagging edges with end-point coordinates lower than z coord
  , tag::zminus,  kw::amr_zminus::info::expect::type
  //! Refinement tagging edges with end-point coordinates higher than z coord
  , tag::zplus,  kw::amr_zplus::info::expect::type
> >;

//! A list of side sets moving with a user-defined function in time
using moving_sides = tk::TaggedTuple< brigand::list<
  //! List of side sets to move
     tag::sideset, std::vector< kw::sideset::info::expect::type >
  //! User-defined table (function) type
  ,  tag::fntype,  tk::ctr::UserTableType
  //! Functions x(t), y(t), and z(t) to move the side sets with
  ,  tag::fn,      std::vector< tk::real >
> >;

//! A list of side sets along with a user-defined function for time dependent BC
using time_dependent_bc = tk::TaggedTuple< brigand::list<
  //! List of side sets on which to apply time dependent BC
     tag::sideset, std::vector< kw::sideset::info::expect::type >
  //! Functions p(t), rho(t), u(t), v(t) and w(t) to specify time dependent BC
  ,  tag::fn,      std::vector< tk::real >
> >;

//! ALE mesh motion options
using ale = tk::TaggedTuple< brigand::list<
  //! ALE on/off
    tag::ale,           bool
  //! Restrict mesh velocity dimensions (useful for d<3 dimensional problems)
  , tag::mesh_motion,   std::vector< kw::mesh_motion::info::expect::type >
  //! dvCFL (CFL mesh volume change) coefficient for ALE
  , tag::dvcfl,         kw::dvcfl::info::expect::type
  //!< Multiplier for vorticity in mesh velocity smoother
  , tag::vortmult,      kw::vortmult::info::expect::type
  //! Mesh velocity smoother linear solver max number of iterations
  , tag::maxit,         kw::meshvel_maxit::info::expect::type
  //! Mesh velocity smoother linear solver tolerance
  , tag::tolerance,     kw::meshvel_tolerance::info::expect::type
  //! Mesh velocity option
  , tag::meshvelocity,  MeshVelocityType
  //! Mesh velocity smoother option
  , tag::smoother    ,  MeshVelocitySmootherType
    //! Mesh velocity Dirichlet BC sidesets
  , tag::bcdir,         std::vector< kw::sideset::info::expect::type >
    //! Mesh velocity symmetry BC sidesets
  , tag::bcsym,         std::vector< kw::sideset::info::expect::type >
    //! Mesh force parameters
  , tag::meshforce,     std::vector< kw::meshforce::info::expect::type >
    //! List of side sets to move with a user-defined function
  , tag::move,          std::vector< moving_sides >
> >;

//! p-adaptive refinement options
using pref = tk::TaggedTuple< brigand::list<
    tag::pref,        bool                //!< p-refinement on/off
  , tag::indicator,   PrefIndicatorType   //!< Choice of adaptive indicator
  , tag::ndofmax,     std::size_t         //!< Max number of degree of freedom
  , tag::tolref,      tk::real            //!< Threshold of p-refinement
> >;

//! Discretization parameters storage
using discretization = tk::TaggedTuple< brigand::list<
    tag::nstep,  kw::nstep::info::expect::type  //!< Number of time steps
  , tag::term,   kw::term::info::expect::type   //!< Time to terminate
  , tag::t0,     kw::t0::info::expect::type     //!< Starting time
  , tag::dt,     kw::dt::info::expect::type     //!< Size of time step
  , tag::cfl,    kw::cfl::info::expect::type    //!< CFL coefficient
  , tag::pelocal_reorder, bool                  //!< PE-locality reordering
  , tag::operator_reorder, bool                 //!< Operator-access reordering
  , tag::steady_state, bool                     //!< March to steady state
  , tag::residual, kw::residual::info::expect::type //!< Convergence residual
  , tag::rescomp, kw::rescomp::info::expect::type //!< Convergence residual comp
  , tag::fct,    bool                           //!< FCT on/off
  , tag::fctclip,bool                           //!< FCT clipping limiter on/off
  , tag::fcteps, kw::fcteps::info::expect::type //!< FCT small number
  , tag::ctau,   kw::ctau::info::expect::type   //!< FCT mass diffisivity
  , tag::scheme, inciter::ctr::SchemeType       //!< Spatial discretization type
  , tag::limiter,inciter::ctr::LimiterType      //!< Limiter type
  , tag::cweight,kw::cweight::info::expect::type//!< WENO central stencil weight
  , tag::rdof,   std::size_t          //!< Number of reconstructed solution DOFs
  , tag::ndof,   std::size_t                    //!< Number of solution DOFs
  , tag::accuracy_test, bool                    //!< Accuracy test on/off
  , tag::limsol_projection, bool         //!< Limited solution projection on/off
  , tag::shock_detector_coeff, kw::shock_detector_coeff::info::expect::type
      //!< Shock detector coefficient for limiting on/off
> >;

//! ASCII output floating-point precision in digits
using precision = tk::TaggedTuple< brigand::list<
    //! Diagnostics output precision
    tag::diag, kw::precision::info::expect::type
    //! History output precision
  , tag::history, kw::precision::info::expect::type
> >;

//! ASCII output floating-point format
using floatformat = tk::TaggedTuple< brigand::list<
    tag::diag,    tk::ctr::TxtFloatFormatType  //!< Diagnostics output format
  , tag::history, tk::ctr::TxtFloatFormatType  //!< History output format
> >;

//! Output intervals in units of iteration count
using interval_iter = tk::TaggedTuple< brigand::list<
    //! TTY output interval
    tag::tty,     kw::ttyi::info::expect::type
    //! Field output interval
  , tag::field,   kw::interval_iter::info::expect::type
    //! History output interval
  , tag::history, kw::interval_iter::info::expect::type
    //! Diags output interval
  , tag::diag,    kw::interval_iter::info::expect::type
> >;

//! Output intervals in units of physics time
using interval_time = tk::TaggedTuple< brigand::list<
    //! Field output interval
    tag::field,   kw::interval_time::info::expect::type
    //! History output interval
  , tag::history, kw::interval_time::info::expect::type
> >;

//! Output time ranges in units of physics time
using time_range = tk::TaggedTuple< brigand::list<
    //! \brief Field output configuration: outer vector: multiple ranges, inner
    //!        vector: mintime, maxtime, dt
    tag::field,   std::vector<
                    std::vector< kw::time_range::info::expect::type > >
    //! \brief History output configuration: outer vector: multiple ranges,
    //!        inner vector: mintime, maxtime, dt
  , tag::history, std::vector<
                    std::vector< kw::time_range::info::expect::type > >
> >;

//! Output configuration parameters
using output_parameters = tk::TaggedTuple< brigand::list<
    //! Output intervals in units of iteration count
    tag::iter,  interval_iter
    //! Output intervals in units of physics time
  , tag::time,  interval_time
    //! Output time ranges in units of physics time
  , tag::range, time_range
> >;

//! History output parameters storage
using history = tk::TaggedTuple< brigand::list<
    tag::point,   std::vector< std::vector< kw::point::info::expect::type > >
  , tag::id,      std::vector< std::string >     //!< Point identifiers
> >;

//! IO parameters storage
using ios = tk::TaggedTuple< brigand::list<
    tag::nrestart,  int                             //!< Number of restarts
  , tag::control,   kw::control::info::expect::type //!< Control filename
  , tag::input,     kw::input::info::expect::type   //!< Input filename
  , tag::output,    kw::output::info::expect::type  //!< Output filename
    //! Refined output (output field data on a refined mesh)
  , tag::refined,   kw::refined::info::expect::type
  , tag::screen,    kw::screen::info::expect::type  //!< Screen output filename
    //! List of side sets to save as field output
  , tag::surface,   std::vector< kw::sideset::info::expect::type >
    //! Diagnostics filename
  , tag::diag,      kw::diagnostics_cmd::info::expect::type
  , tag::particles, std::string                     //!< Particles filename
  , tag::outvar,    std::vector< OutVar >           //!< Output variables
  , tag::restart,   kw::restart::info::expect::type //!< Restart dirname
> >;

//! Error/diagnostics output configuration
using diagnostics = tk::TaggedTuple< brigand::list<
  tag::error,       std::vector< tk::ctr::ErrorType > //!< Errors to compute
> >;

//! Initiation configuration for box IC
using InitiateParameters = tk::TaggedTuple< brigand::list<
    tag::init,          InitiateType
  , tag::point,         std::vector< kw::point::info::expect::type >
  , tag::front_width,   kw::front_width::info::expect::type
  , tag::velocity,      kw::velocity::info::expect::type
> >;

//! Box, given by coordinates, specifying physics variables
using box = tk::TaggedTuple< brigand::list<
    tag::xmin,          kw::xmin::info::expect::type
  , tag::xmax,          kw::xmax::info::expect::type
  , tag::ymin,          kw::ymin::info::expect::type
  , tag::ymax,          kw::ymax::info::expect::type
  , tag::zmin,          kw::zmin::info::expect::type
  , tag::zmax,          kw::zmax::info::expect::type
  , tag::orientation,   std::vector< kw::orientation::info::expect::type >
  , tag::materialid,    kw::materialid::info::expect::type
  , tag::mass,          kw::mass::info::expect::type
  , tag::density,       kw::density::info::expect::type
  , tag::velocity,      std::vector< kw::velocity::info::expect::type >
  , tag::pressure,      kw::pressure::info::expect::type
  , tag::energy,        kw::energy::info::expect::type
  , tag::energy_content,kw::energy_content::info::expect::type
  , tag::temperature,   kw::temperature::info::expect::type
  , tag::initiate,      InitiateParameters
> >;

//! Mesh block from mesh file specifying physics variables
using meshblock = tk::TaggedTuple< brigand::list<
    tag::blockid,       kw::blockid::info::expect::type
  , tag::materialid,    kw::materialid::info::expect::type
  , tag::volume,        kw::volume::info::expect::type
  , tag::mass,          kw::mass::info::expect::type
  , tag::density,       kw::density::info::expect::type
  , tag::velocity,      std::vector< kw::velocity::info::expect::type >
  , tag::pressure,      kw::pressure::info::expect::type
  , tag::energy,        kw::energy::info::expect::type
  , tag::energy_content,kw::energy_content::info::expect::type
  , tag::temperature,   kw::temperature::info::expect::type
  , tag::initiate,      InitiateParameters
> >;

//! Initial condition configuration
using ic = tk::TaggedTuple< brigand::list<
    tag::density,       std::vector<
                          std::vector< kw::density::info::expect::type > >
  , tag::materialid,    std::vector<
                          std::vector< kw::materialid::info::expect::type > >
  , tag::velocity,      std::vector<
                          std::vector< kw::velocity::info::expect::type > >
  , tag::pressure,      std::vector<
                          std::vector< kw::pressure::info::expect::type > >
  , tag::energy,        std::vector<
                          std::vector< kw::energy::info::expect::type > >
  , tag::temperature,   std::vector<
                          std::vector< kw::temperature::info::expect::type > >
  , tag::box,           std::vector< std::vector< box > >
  , tag::meshblock,     std::vector< std::vector< meshblock > >
> >;

//! Boundary conditions configuration (list of side sets for each eq system)
using bc = tk::TaggedTuple< brigand::list<
    tag::bcdir,             std::vector< std::vector<
                              kw::sideset::info::expect::type > >
  , tag::bcsym,             std::vector< std::vector<
                              kw::sideset::info::expect::type > >
  , tag::bcinlet,           std::vector< std::vector<
                              kw::sideset::info::expect::type > >
  , tag::bcoutlet,          std::vector< std::vector<
                              kw::sideset::info::expect::type > >
  , tag::bcfarfield,        std::vector< std::vector<
                              kw::sideset::info::expect::type > >
  , tag::bcextrapolate,     std::vector< std::vector<
                              kw::sideset::info::expect::type > >
> >;

//! Solver coupling
using couple = tk::TaggedTuple< brigand::list<
    tag::transfer,  std::vector< Transfer >     //!< List of mesh transfers
> >;

//! Mesh assignment and configuration
using mesh = tk::TaggedTuple< brigand::list<
    tag::id,          std::vector< std::size_t >
  , tag::filename,    std::vector< std::string >
  , tag::location,    std::vector<
                        std::vector< kw::location::info::expect::type > >
  , tag::orientation, std::vector<
                        std::vector< kw::orientation::info::expect::type > >
  , tag::reference,   std::vector< char >
> >;

//! Stagnation points parameters storage
using StagnationParameters = tk::TaggedTuple< brigand::list<
    tag::point,         std::vector<
                          std::vector< kw::point::info::expect::type > >
  , tag::radius,        std::vector<
                          std::vector< kw::radius::info::expect::type > >
> >;

//! Skip points parameters storage
using SkipParameters = tk::TaggedTuple< brigand::list<
    tag::point,         std::vector<
                          std::vector< kw::point::info::expect::type > >
  , tag::radius,        std::vector<
                          std::vector< kw::radius::info::expect::type > >
> >;

//! Sponge parameters storage
using SpongeParameters = tk::TaggedTuple< brigand::list<
    tag::sideset,       std::vector< std::vector<
                          kw::sideset::info::expect::type > >
  , tag::velocity,      std::vector< std::vector<
                          kw::velocity::info::expect::type > >
  , tag::pressure,      std::vector< std::vector<
                          kw::pressure::info::expect::type > >
> >;

//! Transport equation parameters storage
using TransportPDEParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::mesh,          mesh
  , tag::physics,       std::vector< PhysicsType >
  , tag::problem,       std::vector< ProblemType >
  , tag::diffusivity,   std::vector< std::vector<
                        kw::pde_diffusivity::info::expect::type > >
  , tag::lambda,        std::vector< std::vector<
                        kw::pde_lambda::info::expect::type > >
  , tag::u0,            std::vector< std::vector<
                        kw::pde_u0::info::expect::type > >
  , tag::bc,            bc
  , tag::bctimedep,     std::vector< std::vector< time_dependent_bc > >
  , tag::sponge,        SpongeParameters
  //! interface compression toggle
  , tag::intsharp,      std::vector< kw::intsharp::info::expect::type >
  //! interface compression parameter
  , tag::intsharp_param,
                      std::vector< kw::intsharp_param::info::expect::type >
> >;

//! Material configuration
using material = tk::TaggedTuple< brigand::list<
    //! material id
    tag::id,       std::vector< kw::id::info::expect::type >
    //! material EOS
  , tag::eos,      MaterialType
    //! Ratio of spec heats
  , tag::gamma,    std::vector< kw::mat_gamma::info::expect::type >
    //! EoS stiffness parameter
  , tag::pstiff,   std::vector< kw::mat_pstiff::info::expect::type >
    //! Gruneisen coefficient
  , tag::w_gru,    std::vector< kw::w_gru::info::expect::type >
    //! JWL EoS parameter A
  , tag::A_jwl,    std::vector< kw::A_jwl::info::expect::type >
    //! JWL EoS parameter B
  , tag::B_jwl,    std::vector< kw::B_jwl::info::expect::type >
    //! JWL EoS parameter C
  , tag::C_jwl,    std::vector< kw::C_jwl::info::expect::type >
    //! JWL EoS parameter R1
  , tag::R1_jwl,   std::vector< kw::R1_jwl::info::expect::type >
    //! JWL EoS parameter R2
  , tag::R2_jwl,   std::vector< kw::R2_jwl::info::expect::type >
    //! JWL EoS parameter rho0
  , tag::rho0_jwl, std::vector< kw::rho0_jwl::info::expect::type >
    //! JWL EoS parameter e0
  , tag::e0_jwl,   std::vector< kw::e0_jwl::info::expect::type >
    //! Dynamic viscosity
  , tag::mu,       std::vector< kw::mat_mu::info::expect::type >
    //! Spec. heat at const vol.
  , tag::cv,       std::vector< kw::mat_cv::info::expect::type >
    //! Heat conductivity
  , tag::k,        std::vector< kw::mat_k::info::expect::type >
> >;

//! Compressible flow equation parameters storage
using CompFlowPDEParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::mesh,          mesh
  , tag::physics,       std::vector< PhysicsType >
  , tag::problem,       std::vector< ProblemType >
  , tag::farfield_pressure, std::vector< kw::pressure::info::expect::type >
  , tag::farfield_density,  std::vector< kw::density::info::expect::type >
  , tag::farfield_velocity, std::vector< std::vector<
                              kw::velocity::info::expect::type > >
  , tag::bc,            bc
  , tag::bctimedep,     std::vector< std::vector< time_dependent_bc > >
  , tag::sponge,        SpongeParameters
  , tag::ic,            ic
  //! Stagnation boundary condition configuration storage
  , tag::stag,          StagnationParameters
  //! Skip boundary condition configuration storage
  , tag::skip,          SkipParameters
  //! System FCT character
  , tag::sysfct,        std::vector< int >
  //! Indices of system-FCT scalar components considered as a system
  , tag::sysfctvar,     std::vector<
                          std::vector< kw::sysfctvar::info::expect::type > >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::alpha,         std::vector< kw::pde_alpha::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::beta,          std::vector< kw::pde_beta::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::betax,         std::vector< kw::pde_betax::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::betay,         std::vector< kw::pde_betay::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::betaz,         std::vector< kw::pde_betaz::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::r0,            std::vector< kw::pde_r0::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::ce,            std::vector< kw::pde_ce::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::kappa,         std::vector< kw::pde_kappa::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::p0,            std::vector< kw::pde_p0::info::expect::type >
    //! Materials block
  , tag::material,      std::vector< std::vector< material > >
    //! Materials index/EoS map
  , tag::matidxmap,     tk::TaggedTuple< brigand::list<
      tag::eosidx,      std::vector< std::size_t >,
      tag::matidx,      std::vector< std::size_t > > >
    //! Flux function type
  , tag::flux,          std::vector< FluxType >
    //! Lua code (multiple blocks)
  , tag::lua,           std::vector< std::string >
> >;

//! Compressible flow equation parameters storage
using MultiMatPDEParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::mesh,          mesh
  , tag::physics,       std::vector< PhysicsType >
  , tag::problem,       std::vector< ProblemType >
  , tag::bc,            bc
  , tag::bctimedep,     std::vector< std::vector< time_dependent_bc > >
  , tag::ic,            ic
  , tag::farfield_pressure, std::vector< kw::pressure::info::expect::type >
  , tag::farfield_density,  std::vector< kw::density::info::expect::type >
  , tag::farfield_velocity, std::vector< std::vector<
                              kw::velocity::info::expect::type > >
  , tag::sponge,        SpongeParameters
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::alpha,         std::vector< kw::pde_alpha::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::beta,          std::vector< kw::pde_beta::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::betax,         std::vector< kw::pde_betax::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::betay,         std::vector< kw::pde_betay::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::betaz,         std::vector< kw::pde_betaz::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::r0,            std::vector< kw::pde_r0::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::ce,            std::vector< kw::pde_ce::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::kappa,         std::vector< kw::pde_kappa::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::p0,            std::vector< kw::pde_p0::info::expect::type >
    //! Materials block
  , tag::material,      std::vector< std::vector< material > >
    //! Materials index/EoS map
    //!   The following data structure is used to index into the correct
    //!   material vector entry. This is done using the following three maps:
    //!   1. eosidx: This vector provides the eos-index (value) in the
    //!   vector<tag::material> for the given user-spec material id (index).
    //!   2. matidx: This vector provides the material-index (value) inside the
    //!   vector<tag::material>[eosidx] block for the given user-specified
    //!   material id (index).
  , tag::matidxmap,     tk::TaggedTuple< brigand::list<
      tag::eosidx,      std::vector< std::size_t >,
      tag::matidx,      std::vector< std::size_t > > >
    //! number of materials
  , tag::nmat,          std::vector< kw::nmat::info::expect::type >
    //! pressure relaxation toggle
  , tag::prelax,        std::vector< kw::prelax::info::expect::type >
    //! pressure relaxation time scale
  , tag::prelax_timescale,
                      std::vector< kw::prelax_timescale::info::expect::type >
    //! interface compression toggle
  , tag::intsharp,      std::vector< kw::intsharp::info::expect::type >
    //! interface compression parameter
  , tag::intsharp_param,
                      std::vector< kw::intsharp_param::info::expect::type >
    //! Flux function type
  , tag::flux,          std::vector< FluxType >
> >;

//! Parameters storage
using parameters = tk::TaggedTuple< brigand::list<
    tag::transport, TransportPDEParameters
  , tag::compflow,  CompFlowPDEParameters
  , tag::multimat,  MultiMatPDEParameters
> >;

//! PEGTL location/position type to use throughout all of Inciter's parsers
using Location = pegtl::position;

} // ctr::
} // inciter::

#endif // IncitierTypes_h
