// *****************************************************************************
/*!
  \file      src/Control/Inciter/Types.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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
#include "Inciter/Options/AMRInitial.hpp"
#include "Inciter/Options/AMRError.hpp"
#include "Inciter/Options/PrefIndicator.hpp"
#include "Options/PartitioningAlgorithm.hpp"
#include "Options/TxtFloatFormat.hpp"
#include "Options/FieldFile.hpp"
#include "Options/Error.hpp"
#include "PUPUtil.hpp"

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
    tag::amr,     bool                             //!< AMR on/off
  , tag::t0ref,   bool                             //!< AMR before t<0 on/off
  , tag::dtref,   bool                             //!< AMR during t>0 on/off
  , tag::dtref_uniform, bool                       //!< Force dtref uniform-only
  , tag::dtfreq,  kw::amr_dtfreq::info::expect::type //!< Refinement frequency
  , tag::init,    std::vector< AMRInitialType >    //!< List of initial AMR types
  , tag::refvar,  std::vector< std::string >       //!< List of refinement vars
  , tag::id,      std::vector< std::size_t >       //!< List of refvar indices
  , tag::error,   AMRErrorType                     //!< Error estimator for AMR
  , tag::tolref,  tk::real                         //!< Refine tolerance
  , tag::tolderef, tk::real                        //!< De-refine tolerance
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
  , tag::fct,    bool                           //!< FCT on/off
  , tag::reorder,bool                           //!< reordering on/off
  , tag::ctau,   kw::ctau::info::expect::type   //!< FCT mass diffisivity
  , tag::scheme, inciter::ctr::SchemeType       //!< Spatial discretization type
  , tag::limiter,inciter::ctr::LimiterType      //!< Limiter type
  , tag::cweight,kw::cweight::info::expect::type//!< WENO central stencil weight
  , tag::flux,   inciter::ctr::FluxType         //!< Flux function type
  , tag::rdof,   std::size_t          //!< Number of reconstructed solution DOFs
  , tag::ndof,   std::size_t                   //!< Number of solution DOFs
> >;

//! ASCII output floating-point precision in digits
using precision = tk::TaggedTuple< brigand::list<
    tag::diag, kw::precision::info::expect::type //!< Diagnostics output precision
> >;

//! ASCII output floating-point format
using floatformat = tk::TaggedTuple< brigand::list<
    tag::diag, tk::ctr::TxtFloatFormatType  //!< Diagnostics output format
> >;

//! Output intervals storage
using intervals = tk::TaggedTuple< brigand::list<
    tag::tty,   kw::ttyi::info::expect::type      //!< TTY output interval
  , tag::field, kw::interval::info::expect::type  //!< Field output interval
  , tag::diag,  kw::interval::info::expect::type  //!< Diags output interval
> >;

//! IO parameters storage
using ios = tk::TaggedTuple< brigand::list<
    tag::control,     kw::control::info::expect::type //!< Control filename
  , tag::input,       std::string                     //!< Input filename
  , tag::output,      std::string                     //!< Output filename
  , tag::diag,        std::string                     //!< Diagnostics filename
  , tag::part,        std::string                     //!< Particles filename
  , tag::restart,     std::string                     //!< Restart dirname
> >;

//! Error/diagnostics output configuration
using diagnostics = tk::TaggedTuple< brigand::list<
  tag::error,       std::vector< tk::ctr::ErrorType > //!< Errors to compute
> >;

//! Transport equation parameters storage
using TransportPDEParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::physics,       std::vector< PhysicsType >
  , tag::problem,       std::vector< ProblemType >
  , tag::diffusivity,   std::vector< std::vector<
                        kw::pde_diffusivity::info::expect::type > >
  , tag::lambda,        std::vector< std::vector<
                        kw::pde_lambda::info::expect::type > >
  , tag::u0,            std::vector< std::vector<
                        kw::pde_u0::info::expect::type > >
  , tag::bcdir,         std::vector< std::vector<
                         kw::sideset::info::expect::type > >
  , tag::bcsym,         std::vector< std::vector<
                         kw::sideset::info::expect::type > >
  , tag::bcinlet,       std::vector< std::vector<
                         kw::sideset::info::expect::type > >
  , tag::bcoutlet,      std::vector< std::vector<
                         kw::sideset::info::expect::type > >
  , tag::bcextrapolate, std::vector< std::vector<
                         kw::sideset::info::expect::type > >
> >;

//! Compressible flow equation parameters storage
using CompFlowPDEParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::physics,       std::vector< PhysicsType >
  , tag::problem,       std::vector< ProblemType >
  , tag::bcdir,         std::vector< std::vector<
                          kw::sideset::info::expect::type > >
  , tag::bcsym,         std::vector< std::vector<
                          kw::sideset::info::expect::type > >
  , tag::bcinlet,       std::vector< std::vector<
                          kw::sideset::info::expect::type > >
  , tag::bcsubsonicoutlet,
                        std::vector< std::vector<
                          kw::sideset::info::expect::type > >
  , tag::farfield_pressure,
                        std::vector< kw::farfield_pressure::info::expect::type >
  , tag::bcextrapolate, std::vector< std::vector<
                         kw::sideset::info::expect::type > >
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
    //! Ratio of spec heats
  , tag::gamma,         std::vector<
                          std::vector< kw::mat_gamma::info::expect::type > >
    //! EoS stiffness parameter
  , tag::pstiff,        std::vector<
                          std::vector< kw::mat_pstiff::info::expect::type > >
    //! Dynamic viscosity
  , tag::mu,            std::vector<
                          std::vector< kw::mat_mu::info::expect::type > >
    //! Spec. heat at const vol.
  , tag::cv,            std::vector<
                          std::vector< kw::mat_cv::info::expect::type > >
    //! Heat conductivity
  , tag::k,             std::vector<
                          std::vector< kw::mat_k::info::expect::type > >
    //! total number of optional passive tracker particles for visualization
  , tag::npar,          std::vector< kw::npar::info::expect::type >
> >;

//! Compressible flow equation parameters storage
using MultiMatPDEParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::physics,       std::vector< PhysicsType >
  , tag::problem,       std::vector< ProblemType >
  , tag::bcdir,         std::vector< std::vector<
                       kw::sideset::info::expect::type > >
  , tag::bcsym,         std::vector< std::vector<
                       kw::sideset::info::expect::type > >
  , tag::bcinlet,       std::vector< std::vector<
                        kw::sideset::info::expect::type > >
  , tag::bcoutlet,      std::vector< std::vector<
                        kw::sideset::info::expect::type > >
  , tag::bcextrapolate, std::vector< std::vector<
                         kw::sideset::info::expect::type > >
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
    //! Ratio of spec heats
  , tag::gamma,         std::vector<
                          std::vector< kw::mat_gamma::info::expect::type > >
    //! EoS stiffness parameter
  , tag::pstiff,        std::vector<
                          std::vector< kw::mat_pstiff::info::expect::type > >
    //! Dynamic viscosity
  , tag::mu,            std::vector<
                          std::vector< kw::mat_mu::info::expect::type > >
    //! Spec. heat at const vol.
  , tag::cv,            std::vector<
                          std::vector< kw::mat_cv::info::expect::type > >
    //! Heat conductivity
  , tag::k,             std::vector<
                          std::vector< kw::mat_k::info::expect::type > >
  //! number of materials
  , tag::nmat,          std::vector< kw::nmat::info::expect::type >
> >;

//! Parameters storage
using parameters = tk::TaggedTuple< brigand::list<
    tag::transport,     TransportPDEParameters
  , tag::compflow,      CompFlowPDEParameters
  , tag::multimat,      MultiMatPDEParameters
> >;

//! PEGTL location/position type to use throughout all of Inciter's parsers
using Location = pegtl::position;

} // ctr::
} // inciter::

#endif // IncitierTypes_h
