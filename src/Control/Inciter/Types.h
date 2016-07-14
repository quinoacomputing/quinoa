// *****************************************************************************
/*!
  \file      src/Control/Inciter/Types.h
  \author    J. Bakosi
  \date      Mon 11 Jul 2016 10:54:08 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Types for Incitier's parsers
  \details   Types for Incitier's parsers. This file defines the components of the
    tagged tuple that stores heteroegeneous objects in a hierarchical way. These
    components are therefore part of the grammar stack that is filled during
    parsing (both command-line argument parsing and control file parsing).
*/
// *****************************************************************************
#ifndef IncitierTypes_h
#define IncitierTypes_h

#include "Tags.h"
#include "Types.h"
#include "Inciter/Options/PDE.h"
#include "Inciter/Options/Problem.h"
#include "Options/PartitioningAlgorithm.h"
#include "PUPUtil.h"

namespace inciter {
namespace ctr {

//! Storage of selected options
using selects = tk::tuple::tagged_tuple<
  tag::pde,          std::vector< ctr::PDEType >,       //!< Partial diff eqs
  tag::partitioner,  tk::ctr::PartitioningAlgorithmType //!< Mesh partitioner
>;

//! Discretization parameters storage
using discretization = tk::tuple::tagged_tuple<
  tag::nstep,     kw::nstep::info::expect::type, //!< Number of time steps
  tag::term,      kw::term::info::expect::type,  //!< Time to terminate
  tag::t0,        kw::t0::info::expect::type,    //!< Starting time
  tag::dt,        kw::dt::info::expect::type     //!< Size of time step
>;

//! Output intervals storage
using intervals = tk::tuple::tagged_tuple<
  tag::tty,   kw::ttyi::info::expect::type,       //!< TTY output interval
  tag::field, kw::interval::info::expect::type    //!< Field output interval
>;

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::control,     kw::control::info::expect::type,  //!< Control filename
  tag::input,       std::string,                      //!< Input filename
  tag::output,      std::string                       //!< Output filename
>;

//! Advection-diffusion transport equation parameters storage
using AdvDiffPDEParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::problem,     std::vector< ProblemType >,
  tag::diffusivity, std::vector< std::vector<
                      kw::pde_diffusivity::info::expect::type > >,
  tag::lambda,      std::vector< std::vector<
                      kw::pde_lambda::info::expect::type > >,
  tag::u0,          std::vector< std::vector<
                      kw::pde_u0::info::expect::type > >,
  tag::bc_dirichlet, std::vector< std::vector<
                       kw::sideset::info::expect::type > >
>;

//! Poisson equation parameters storage
using PoissonPDEParameters = tk::tuple::tagged_tuple<
  tag::depvar,       std::vector< char >,
  tag::problem,      std::vector< ProblemType >,
  tag::bc_dirichlet, std::vector< std::vector<
                       kw::sideset::info::expect::type > >
>;

//! Euler equation parameters storage
using EulerPDEParameters = tk::tuple::tagged_tuple<
  tag::problem,     std::vector< ProblemType >,
  tag::bc_dirichlet, std::vector< std::vector<
                       kw::sideset::info::expect::type > >
>;

//! Parameters storage
using parameters = tk::tuple::tagged_tuple<
  tag::advdiff,     AdvDiffPDEParameters,
  tag::poisson,     PoissonPDEParameters,
  tag::euler,       EulerPDEParameters
>;

//! PEGTL location type to use throughout Incitier's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // inciter::

#endif // IncitierTypes_h
