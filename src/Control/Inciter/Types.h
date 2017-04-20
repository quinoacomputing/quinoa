// *****************************************************************************
/*!
  \file      src/Control/Inciter/Types.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
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

#include "Tags.h"
#include "Types.h"
#include "Inciter/Options/PDE.h"
#include "Inciter/Options/Problem.h"
#include "Inciter/Options/InitialAMR.h"
#include "Options/PartitioningAlgorithm.h"
#include "Options/TxtFloatFormat.h"
#include "PUPUtil.h"

namespace inciter {
namespace ctr {

//! Storage of selected options
using selects = tk::tuple::tagged_tuple<
  tag::pde,          std::vector< ctr::PDEType >,       //!< Partial diff eqs
  tag::partitioner,  tk::ctr::PartitioningAlgorithmType,//!< Mesh partitioner
  tag::initialamr,   tk::ctr::InitialAMRType            //!< Initial AMR type
>;

//! Discretization parameters storage
using discretization = tk::tuple::tagged_tuple<
  tag::nstep,     kw::nstep::info::expect::type, //!< Number of time steps
  tag::term,      kw::term::info::expect::type,  //!< Time to terminate
  tag::t0,        kw::t0::info::expect::type,    //!< Starting time
  tag::dt,        kw::dt::info::expect::type,    //!< Size of time step
  tag::cfl,       kw::cfl::info::expect::type,   //!< CFL coefficient
  tag::ctau,      kw::ctau::info::expect::type   //!< FCT mass diffisivity
>;

//! ASCII output floating-point precision in digits
using precision = tk::tuple::tagged_tuple<
  tag::diag, kw::precision::info::expect::type //!< Diagnostics output precision
>;

//! ASCII output floating-point format
using floatformat = tk::tuple::tagged_tuple<
  tag::diag, tk::ctr::TxtFloatFormatType  //!< Diagnostics output format
>;

//! Output intervals storage
using intervals = tk::tuple::tagged_tuple<
  tag::tty,   kw::ttyi::info::expect::type,       //!< TTY output interval
  tag::field, kw::interval::info::expect::type,   //!< Field output interval
  tag::diag,  kw::interval::info::expect::type    //!< Diags output interval
>;

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::control,     kw::control::info::expect::type,  //!< Control filename
  tag::input,       std::string,                      //!< Input filename
  tag::output,      std::string,                      //!< Output filename
  tag::diag,        std::string,                      //!< Diagnostics filename
  tag::part,        std::string                       //!< Particles filename
>;

//! Transport equation parameters storage
using TransportPDEParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::physics,      std::vector< PhysicsType >,
  tag::problem,     std::vector< ProblemType >,
  tag::diffusivity, std::vector< std::vector<
                      kw::pde_diffusivity::info::expect::type > >,
  tag::lambda,      std::vector< std::vector<
                      kw::pde_lambda::info::expect::type > >,
  tag::u0,          std::vector< std::vector<
                      kw::pde_u0::info::expect::type > >,
  tag::bcdir,       std::vector< std::vector<
                       kw::sideset::info::expect::type > >
>;

//! Poisson equation parameters storage
using PoissonPDEParameters = tk::tuple::tagged_tuple<
  tag::depvar,       std::vector< char >,
  tag::physics,      std::vector< PhysicsType >,
  tag::problem,      std::vector< ProblemType >,
  tag::bcdir,        std::vector< std::vector<
                       kw::sideset::info::expect::type > >
>;

//! Compressible flow equation parameters storage
using CompFlowPDEParameters = tk::tuple::tagged_tuple<
  tag::physics,      std::vector< PhysicsType >,
  tag::problem,      std::vector< ProblemType >,
  tag::bcdir,        std::vector< std::vector<
                       kw::sideset::info::expect::type > >,
  //! Parameter vector (for specific, e.g., verification, problems)
  tag::alpha, std::vector< kw::pde_alpha::info::expect::type >,
  //! Parameter vector (for specific, e.g., verification, problems)
  tag::beta, std::vector< kw::pde_beta::info::expect::type >,
  //! Parameter vector (for specific, e.g., verification, problems)
  tag::p0, std::vector< kw::pde_p0::info::expect::type >,
  //! Material ID
  tag::id,    std::vector< kw::id::info::expect::type >,
  //! Ratio of spec heats
  tag::gamma, std::vector< kw::mat_gamma::info::expect::type >,
  //! Dynamic viscosity
  tag::mu,    std::vector< kw::mat_mu::info::expect::type >,
  //! Spec. heat at const vol.
  tag::cv,    std::vector< kw::mat_cv::info::expect::type >,
  //! Heat conductivity
  tag::k,     std::vector< kw::mat_k::info::expect::type >,
  //! total number of optional passive tracker particles for visualization
  tag::npar,  std::vector< kw::npar::info::expect::type >
>;

//! Parameters storage
using parameters = tk::tuple::tagged_tuple<
  tag::transport,   TransportPDEParameters,
  tag::poisson,     PoissonPDEParameters,
  tag::compflow,    CompFlowPDEParameters
>;

//! PEGTL location/position type to use throughout all of Inciter's parsers
using Location = pegtl::position_info;

} // ctr::
} // inciter::

#endif // IncitierTypes_h
