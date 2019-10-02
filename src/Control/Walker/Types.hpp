// *****************************************************************************
/*!
  \file      src/Control/Walker/Types.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Types for Walker's parsers
  \details   Types for Walker's parsers. This file defines the components of the
    tagged tuple that stores heteroegeneous objects in a hierarchical way. These
    components are therefore part of the grammar stack that is filled during
    parsing (both command-line argument parsing and control file parsing).
*/
// *****************************************************************************
#ifndef WalkerTypes_h
#define WalkerTypes_h

#include "Tags.hpp"
#include "Base/Types.hpp"
#include "RNGParam.hpp"
#include "Walker/Options/DiffEq.hpp"
#include "Walker/Options/InitPolicy.hpp"
#include "Walker/Options/CoeffPolicy.hpp"
#include "Walker/Options/HydroTimeScales.hpp"
#include "Walker/Options/HydroProductions.hpp"
#include "Options/PDFFile.hpp"
#include "Options/PDFPolicy.hpp"
#include "Options/PDFCentering.hpp"
#include "Options/TxtFloatFormat.hpp"
#include "Options/Depvar.hpp"
#include "Options/VelocityVariant.hpp"
#include "Options/Normalization.hpp"
#include "Options/RNG.hpp"
#include "QuinoaConfig.hpp"

namespace walker {
namespace ctr {

//! Storage of selected options
using selects = tk::TaggedTuple< brigand::list<
    tag::diffeq,       std::vector< ctr::DiffEqType >  //!< Differential eqs
  , tag::rng,          std::vector< tk::ctr::RNGType > //!< RNGs
  , tag::filetype,     tk::ctr::PDFFileType      //!< PDF output file type
  , tag::pdfpolicy,    tk::ctr::PDFPolicyType    //!< PDF output file policy
  , tag::pdfctr,       tk::ctr::PDFCenteringType //!< PDF output file centering
> >;

//! Discretization parameters storage
using discretization = tk::TaggedTuple< brigand::list<
    tag::npar,      kw::npar::info::expect::type  //!< Total number of particles
  , tag::nstep,     kw::nstep::info::expect::type   //!< Number of time steps
  , tag::term,      kw::term::info::expect::type    //!< Termination time
  , tag::dt,        kw::dt::info::expect::type      //!< Size of time step
  , tag::binsize,   std::vector< std::vector< tk::real > >  //!< PDF binsizes
  , tag::extent,    std::vector< std::vector< tk::real > >  //!< PDF extents
> >;

//! ASCII output floating-point precision in digits
using precision = tk::TaggedTuple< brigand::list<
    tag::stat, kw::precision::info::expect::type//!< Statistics output precision
  , tag::pdf,  kw::precision::info::expect::type  //!< PDF output precision
> >;

//! ASCII output floating-point format
using floatformat = tk::TaggedTuple< brigand::list<
    tag::stat, tk::ctr::TxtFloatFormatType  //!< Statistics output format
  , tag::pdf,  tk::ctr::TxtFloatFormatType  //!< PDF output format
> >;

//! Output intervals storage
using intervals = tk::TaggedTuple< brigand::list<
    tag::tty,  kw::ttyi::info::expect::type      //!< TTY output interval
  , tag::stat, kw::interval::info::expect::type  //!< Statistics output interval
  , tag::pdf,  kw::interval::info::expect::type  //!< PDF output interval
> >;

//! IO parameters storage
using ios = tk::TaggedTuple< brigand::list<
    tag::control,   kw::control::info::expect::type   //!< Control filename
  , tag::input,     std::string                       //!< Input filename
  , tag::output,    std::string                       //!< Output filename
  , tag::pdf,       kw::pdf::info::expect::type       //!< PDF filename
  , tag::stat,      kw::stat::info::expect::type      //!< Statistics filename
  , tag::pdfnames,  std::vector< std::string >        //!< PDF identifiers
> >;

//! Data for initialization (SDE initial conditions)
using Init = tk::TaggedTuple< brigand::list<
    tag::spike,         std::vector< std::vector< std::vector <
                          kw::spike::info::expect::type > > >
  , tag::betapdf,       std::vector< std::vector< std::vector <
                          kw::betapdf::info::expect::type > > >
  , tag::gamma,         std::vector< std::vector< std::vector <
                          kw::gammapdf::info::expect::type > > >
  , tag::dirichlet,     std::vector< std::vector<
                          kw::dirichletpdf::info::expect::type > >
  , tag::gaussian,      std::vector< std::vector< std::vector <
                          kw::gaussian::info::expect::type > > >
  , tag::mean,          std::vector< std::vector<
                          kw::sde_mu::info::expect::type > >
  , tag::cov,           std::vector< std::vector<
                          kw::sde_cov::info::expect::type > >
  , tag::jointgaussian, std::vector< std::vector< std::vector <
                          kw::gaussian::info::expect::type > > >
> >;

//! Dirichlet parameters storage
using DirichletParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::b,             std::vector< std::vector<
                          kw::sde_b::info::expect::type > >
  , tag::S,             std::vector< std::vector<
                          kw::sde_S::info::expect::type > >
  , tag::kappa,         std::vector< std::vector<
                          kw::sde_kappa::info::expect::type > >
  , tag::init,          Init
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
> >;

//! Mixture Dirichlet parameters storage
using MixDirichletParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::b,             std::vector< std::vector<
                          kw::sde_b::info::expect::type > >
  , tag::S,             std::vector< std::vector<
                          kw::sde_S::info::expect::type > >
  , tag::kappa,         std::vector< std::vector<
                          kw::sde_kappa::info::expect::type > >
  , tag::kappaprime,    std::vector< std::vector<
                          kw::sde_kappaprime::info::expect::type > >
  , tag::rho,           std::vector< std::vector<
                          kw::sde_rho::info::expect::type > >
  , tag::normalization, std::vector< ctr::NormalizationType >
  , tag::init,          Init
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
> >;

//! Generalized Dirichlet parameters storage
using GenDirichletParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::b,             std::vector< std::vector<
                        kw::sde_b::info::expect::type > >
  , tag::S,             std::vector< std::vector<
                        kw::sde_S::info::expect::type > >
  , tag::kappa,         std::vector< std::vector<
                        kw::sde_kappa::info::expect::type > >
  , tag::c,             std::vector< std::vector<
                        kw::sde_c::info::expect::type > >
  , tag::init,          Init
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
> >;

//! Wright-Fisher parameters storage
using WrightFisherParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::omega,         std::vector< std::vector<
                        kw::sde_omega::info::expect::type > >
  , tag::init,          Init
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
> >;

//! Ornstein-Uhlenbeck parameters storage
using OrnsteinUhlenbeckParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::sigmasq,       std::vector< std::vector<
                          kw::sde_sigmasq::info::expect::type > >
  , tag::theta,         std::vector< std::vector<
                          kw::sde_theta::info::expect::type > >
  , tag::mu,            std::vector< std::vector<
                          kw::sde_mu::info::expect::type > >
  , tag::init,          Init
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
> >;

//! Diagonal Ornstein-Uhlenbeck parameters storage
using DiagOrnsteinUhlenbeckParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::sigmasq,       std::vector< std::vector<
                          kw::sde_sigmasq::info::expect::type > >
  , tag::theta,         std::vector< std::vector<
                          kw::sde_theta::info::expect::type > >
  , tag::mu,            std::vector< std::vector<
                          kw::sde_mu::info::expect::type > >
  , tag::init,          Init
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
> >;

//! Skew-normal parameters storage
using SkewNormalParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::timescale,     std::vector< std::vector<
                          kw::sde_T::info::expect::type > >
  , tag::sigmasq,       std::vector< std::vector<
                          kw::sde_sigmasq::info::expect::type > >
  , tag::lambda,        std::vector< std::vector<
                          kw::sde_lambda::info::expect::type > >
  , tag::init,          Init
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
> >;

//! Gamma parameters storage
using GammaParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::b,             std::vector< std::vector<
                          kw::sde_b::info::expect::type > >
  , tag::S,             std::vector< std::vector<
                          kw::sde_S::info::expect::type > >
  , tag::kappa,         std::vector< std::vector<
                          kw::sde_kappa::info::expect::type > >
  , tag::init,          Init
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
> >;

//! Beta parameters storage
using BetaParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::b,             std::vector< std::vector<
                          kw::sde_b::info::expect::type > >
  , tag::S,             std::vector< std::vector<
                          kw::sde_S::info::expect::type > >
  , tag::kappa,         std::vector< std::vector<
                          kw::sde_kappa::info::expect::type > >
  , tag::init,          Init
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
> >;

//! Number-fraction beta parameters storage
using NumberFractionBetaParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::b,             std::vector< std::vector<
                          kw::sde_b::info::expect::type > >
  , tag::S,             std::vector< std::vector<
                          kw::sde_S::info::expect::type > >
  , tag::kappa,         std::vector< std::vector<
                          kw::sde_kappa::info::expect::type > >
  , tag::rho2,          std::vector< std::vector<
                          kw::sde_rho2::info::expect::type > >
  , tag::rcomma,        std::vector< std::vector<
                        kw::sde_rcomma::info::expect::type > >
  , tag::init,          Init
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
> >;

//! Mass-fraction beta parameters storage
using MassFractionBetaParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::b,             std::vector< std::vector<
                          kw::sde_b::info::expect::type > >
  , tag::S,             std::vector< std::vector<
                          kw::sde_S::info::expect::type > >
  , tag::kappa,         std::vector< std::vector<
                          kw::sde_kappa::info::expect::type > >
  , tag::rho2,          std::vector< std::vector<
                          kw::sde_rho2::info::expect::type > >
  , tag::r,             std::vector< std::vector<
                          kw::sde_rcomma::info::expect::type > >
  , tag::init,          Init
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
> >;

//! Mix number-fraction beta parameters storage
using MixNumberFractionBetaParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::bprime,        std::vector< std::vector<
                          kw::sde_bprime::info::expect::type > >
  , tag::S,             std::vector< std::vector<
                          kw::sde_S::info::expect::type > >
  , tag::kappaprime,    std::vector< std::vector<
                          kw::sde_kappaprime::info::expect::type > >
  , tag::rho2,          std::vector< std::vector<
                          kw::sde_rho2::info::expect::type > >
  , tag::rcomma,        std::vector< std::vector<
                          kw::sde_rcomma::info::expect::type > >
  , tag::init,          Init
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
> >;

//! Mix mass-fraction beta parameters storage
using MixMassFractionBetaParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::bprime,        std::vector< std::vector<
                          kw::sde_bprime::info::expect::type > >
  , tag::S,             std::vector< std::vector<
                          kw::sde_S::info::expect::type > >
  , tag::kappaprime,    std::vector< std::vector<
                           kw::sde_kappaprime::info::expect::type > >
  , tag::rho2,          std::vector< std::vector<
                          kw::sde_rho2::info::expect::type > >
  , tag::r,             std::vector< std::vector<
                          kw::sde_r::info::expect::type > >
  , tag::init,          Init
  , tag::mean_gradient, std::vector< std::vector<
                          kw::mean_gradient::info::expect::type > >
  , tag::hydrotimescales, std::vector< std::vector<
                             ctr::HydroTimeScalesType > >
  , tag::hydroproductions, std::vector< std::vector<
                             ctr::HydroProductionsType > >
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
  , tag::solve,         std::vector< ctr::DepvarType >
  , tag::dissipation,   std::vector< char >
  , tag::dissipation_id,std::vector< std::size_t >
  , tag::velocity,      std::vector< char >
  , tag::velocity_id,   std::vector< std::size_t >
> >;

//! Velocity parameters storage
using VelocityParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::c0,            std::vector< kw::sde_c0::info::expect::type >
  , tag::gravity,       std::vector< std::vector<
                          kw::gravity::info::expect::type > >
  , tag::position,      std::vector< char >
  , tag::position_id,   std::vector< std::size_t >
  , tag::dissipation,   std::vector< char >
  , tag::dissipation_id,std::vector< std::size_t >
  , tag::mixmassfracbeta, std::vector< char >
  , tag::mixmassfracbeta_id, std::vector< std::size_t >
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
  , tag::solve,         std::vector< ctr::DepvarType >
  , tag::variant,       std::vector< ctr::VelocityVariantType >
  , tag::init,          Init
  , tag::hydrotimescales, std::vector< std::vector< ctr::HydroTimeScalesType > >
  , tag::hydroproductions,std::vector< std::vector< ctr::HydroProductionsType > >
> >;

//! Position parameters storage
using PositionParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::velocity,      std::vector< char >
  , tag::velocity_id,   std::vector< std::size_t >
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
  , tag::solve,         std::vector< ctr::DepvarType >
  , tag::init,          Init
> >;

//! Dissipation parameters storage
using DissipationParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::velocity,      std::vector< char >
  , tag::velocity_id,   std::vector< std::size_t >
  , tag::c3,            std::vector< kw::sde_c3::info::expect::type >
  , tag::c4,            std::vector< kw::sde_c4::info::expect::type >
  , tag::com1,          std::vector< kw::sde_com1::info::expect::type >
  , tag::com2,          std::vector< kw::sde_com2::info::expect::type >
  , tag::rng,           std::vector< tk::ctr::RNGType >
  , tag::initpolicy,    std::vector< ctr::InitPolicyType >
  , tag::coeffpolicy,   std::vector< ctr::CoeffPolicyType >
  , tag::init,          Init
> >;

//! Parameters storage
using parameters = tk::TaggedTuple< brigand::list<
  #ifdef HAS_MKL
    tag::rngmkl,          tk::ctr::RNGMKLParameters,
  #endif
    tag::rngsse,          tk::ctr::RNGSSEParameters
  , tag::rng123,          tk::ctr::RNGRandom123Parameters
  , tag::dirichlet,       DirichletParameters
  , tag::mixdirichlet,    MixDirichletParameters
  , tag::gendir,          GenDirichletParameters
  , tag::wrightfisher,    WrightFisherParameters
  , tag::ou,              OrnsteinUhlenbeckParameters
  , tag::diagou,          DiagOrnsteinUhlenbeckParameters
  , tag::skewnormal,      SkewNormalParameters
  , tag::gamma,           GammaParameters
  , tag::beta,            BetaParameters
  , tag::numfracbeta,     NumberFractionBetaParameters
  , tag::massfracbeta,    MassFractionBetaParameters
  , tag::mixnumfracbeta,  MixNumberFractionBetaParameters
  , tag::mixmassfracbeta, MixMassFractionBetaParameters
  , tag::velocity,        VelocityParameters
  , tag::position,        PositionParameters
  , tag::dissipation,     DissipationParameters
> >;

} // ctr::
} // walker::

#endif // WalkerTypes_h
