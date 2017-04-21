// *****************************************************************************
/*!
  \file      src/Control/Walker/Types.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Types for Walker's parsers
  \details   Types for Walker's parsers. This file defines the components of the
    tagged tuple that stores heteroegeneous objects in a hierarchical way. These
    components are therefore part of the grammar stack that is filled during
    parsing (both command-line argument parsing and control file parsing).
*/
// *****************************************************************************
#ifndef WalkerTypes_h
#define WalkerTypes_h

#include "Tags.h"
#include "Types.h"
#include "RNGParam.h"
#include "Walker/Options/DiffEq.h"
#include "Walker/Options/InitPolicy.h"
#include "Walker/Options/CoeffPolicy.h"
#include "Walker/Options/HydroTimeScales.h"
#include "Walker/Options/HydroProductions.h"
#include "Options/PDFFile.h"
#include "Options/PDFPolicy.h"
#include "Options/PDFCentering.h"
#include "Options/TxtFloatFormat.h"
#include "Options/RNG.h"
#include "QuinoaConfig.h"

namespace walker {
namespace ctr {

//! Storage of selected options
using selects = tk::tuple::tagged_tuple<
  tag::diffeq,       std::vector< ctr::DiffEqType >,  //!< Differential eqs
  tag::rng,          std::vector< tk::ctr::RNGType >, //!< RNGs
  tag::pdffiletype,  tk::ctr::PDFFileType,      //!< PDF output file type
  tag::pdfpolicy,    tk::ctr::PDFPolicyType,    //!< PDF output file policy
  tag::pdfctr,       tk::ctr::PDFCenteringType  //!< PDF output file centering
>;

//! Discretization parameters storage
using discretization = tk::tuple::tagged_tuple<
  tag::npar,      kw::npar::info::expect::type,   //!< Total number of particles
  tag::nstep,     kw::nstep::info::expect::type,  //!< Number of time steps
  tag::term,      kw::term::info::expect::type,   //!< Termination time
  tag::dt,        kw::dt::info::expect::type,     //!< Size of time step
  tag::binsize,   std::vector< std::vector< tk::real > >, //!< PDF binsizes
  tag::extent,    std::vector< std::vector< tk::real > >  //!< PDF extents
>;

//! ASCII output floating-point precision in digits
using precision = tk::tuple::tagged_tuple<
  tag::stat, kw::precision::info::expect::type, //!< Statistics output precision
  tag::pdf,  kw::precision::info::expect::type  //!< PDF output precision
>;

//! ASCII output floating-point format
using floatformat = tk::tuple::tagged_tuple<
  tag::stat, tk::ctr::TxtFloatFormatType, //!< Statistics output format
  tag::pdf,  tk::ctr::TxtFloatFormatType  //!< PDF output format
>;

//! Output intervals storage
using intervals = tk::tuple::tagged_tuple<
  tag::tty,  kw::ttyi::info::expect::type,      //!< TTY output interval
  tag::stat, kw::interval::info::expect::type,  //!< Statistics output interval
  tag::pdf,  kw::interval::info::expect::type   //!< PDF output interval
>;

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::control,         kw::control::info::expect::type,  //!< Control filename
  tag::input,           std::string,                  //!< Input filename
  tag::output,          std::string,                  //!< Output filename
  tag::pdf,             kw::pdf::info::expect::type,  //!< PDF filename
  tag::stat,            kw::stat::info::expect::type, //!< Statistics filename
  tag::pdfnames,        std::vector< std::string >    //!< PDF identifiers
>;

//! Dirichlet parameters storage
using DirichletParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector<
                      kw::sde_b::info::expect::type > >,
  tag::S,           std::vector< std::vector<
                      kw::sde_S::info::expect::type > >,
  tag::kappa,       std::vector< std::vector<
                      kw::sde_kappa::info::expect::type > >,
  tag::spike,       std::vector< std::vector< std::vector <
                      kw::spike::info::expect::type > > >,
  tag::betapdf,     std::vector< std::vector< std::vector <
                      kw::betapdf::info::expect::type > > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Generalized Dirichlet parameters storage
using GenDirichletParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector<
                      kw::sde_b::info::expect::type > >,
  tag::S,           std::vector< std::vector<
                      kw::sde_S::info::expect::type > >,
  tag::kappa,       std::vector< std::vector<
                      kw::sde_kappa::info::expect::type > >,
  tag::c,           std::vector< std::vector<
                      kw::sde_c::info::expect::type > >,
  tag::spike,       std::vector< std::vector< std::vector <
                      kw::spike::info::expect::type > > >,
  tag::betapdf,     std::vector< std::vector< std::vector <
                      kw::betapdf::info::expect::type > > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Wright-Fisher parameters storage
using WrightFisherParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::omega,       std::vector< std::vector<
                      kw::sde_omega::info::expect::type > >,
  tag::spike,       std::vector< std::vector< std::vector <
                      kw::spike::info::expect::type > > >,
  tag::betapdf,     std::vector< std::vector< std::vector <
                      kw::betapdf::info::expect::type > > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Ornstein-Uhlenbeck parameters storage
using OrnsteinUhlenbeckParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::sigmasq,     std::vector< std::vector<
                      kw::sde_sigmasq::info::expect::type > >,
  tag::theta,       std::vector< std::vector<
                      kw::sde_theta::info::expect::type > >,
  tag::mu,          std::vector< std::vector<
                      kw::sde_mu::info::expect::type > >,
  tag::spike,       std::vector< std::vector< std::vector <
                      kw::spike::info::expect::type > > >,
  tag::betapdf,     std::vector< std::vector< std::vector <
                      kw::betapdf::info::expect::type > > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Diagonal Ornstein-Uhlenbeck parameters storage
using DiagOrnsteinUhlenbeckParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::sigmasq,     std::vector< std::vector<
                      kw::sde_sigmasq::info::expect::type > >,
  tag::theta,       std::vector< std::vector<
                      kw::sde_theta::info::expect::type > >,
  tag::mu,          std::vector< std::vector<
                      kw::sde_mu::info::expect::type > >,
  tag::spike,       std::vector< std::vector< std::vector <
                      kw::spike::info::expect::type > > >,
  tag::betapdf,     std::vector< std::vector< std::vector <
                      kw::betapdf::info::expect::type > > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Skew-normal parameters storage
using SkewNormalParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::timescale,   std::vector< std::vector<
                      kw::sde_T::info::expect::type > >,
  tag::sigmasq,     std::vector< std::vector<
                      kw::sde_sigmasq::info::expect::type > >,
  tag::lambda,      std::vector< std::vector<
                      kw::sde_lambda::info::expect::type > >,
  tag::spike,       std::vector< std::vector< std::vector <
                      kw::spike::info::expect::type > > >,
  tag::betapdf,     std::vector< std::vector< std::vector <
                      kw::betapdf::info::expect::type > > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Gamma parameters storage
using GammaParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector<
                      kw::sde_b::info::expect::type > >,
  tag::S,           std::vector< std::vector<
                      kw::sde_S::info::expect::type > >,
  tag::kappa,       std::vector< std::vector<
                      kw::sde_kappa::info::expect::type > >,
  tag::spike,       std::vector< std::vector< std::vector <
                      kw::spike::info::expect::type > > >,
  tag::betapdf,     std::vector< std::vector< std::vector <
                      kw::betapdf::info::expect::type > > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Beta parameters storage
using BetaParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector<
                      kw::sde_b::info::expect::type > >,
  tag::S,           std::vector< std::vector<
                      kw::sde_S::info::expect::type > >,
  tag::kappa,       std::vector< std::vector<
                      kw::sde_kappa::info::expect::type > >,
  tag::spike,       std::vector< std::vector< std::vector <
                      kw::spike::info::expect::type > > >,
  tag::betapdf,     std::vector< std::vector< std::vector <
                      kw::betapdf::info::expect::type > > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Number-fraction beta parameters storage
using NumberFractionBetaParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector<
                      kw::sde_b::info::expect::type > >,
  tag::S,           std::vector< std::vector<
                      kw::sde_S::info::expect::type > >,
  tag::kappa,       std::vector< std::vector<
                      kw::sde_kappa::info::expect::type > >,
  tag::rho2,        std::vector< std::vector<
                      kw::sde_rho2::info::expect::type > >,
  tag::rcomma,      std::vector< std::vector<
                      kw::sde_rcomma::info::expect::type > >,
  tag::spike,       std::vector< std::vector< std::vector <
                      kw::spike::info::expect::type > > >,
  tag::betapdf,     std::vector< std::vector< std::vector <
                      kw::betapdf::info::expect::type > > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Mass-fraction beta parameters storage
using MassFractionBetaParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector<
                      kw::sde_b::info::expect::type > >,
  tag::S,           std::vector< std::vector<
                      kw::sde_S::info::expect::type > >,
  tag::kappa,       std::vector< std::vector<
                      kw::sde_kappa::info::expect::type > >,
  tag::rho2,        std::vector< std::vector<
                      kw::sde_rho2::info::expect::type > >,
  tag::r,           std::vector< std::vector<
                      kw::sde_rcomma::info::expect::type > >,
  tag::spike,       std::vector< std::vector< std::vector <
                      kw::spike::info::expect::type > > >,
  tag::betapdf,     std::vector< std::vector< std::vector <
                      kw::betapdf::info::expect::type > > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Mix number-fraction beta parameters storage
using MixNumberFractionBetaParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::bprime,      std::vector< std::vector<
                      kw::sde_bprime::info::expect::type > >,
  tag::S,           std::vector< std::vector<
                      kw::sde_S::info::expect::type > >,
  tag::kappaprime, std::vector< std::vector<
                      kw::sde_kappaprime::info::expect::type > >,
  tag::rho2,       std::vector< std::vector<
                      kw::sde_rho2::info::expect::type > >,
  tag::rcomma,     std::vector< std::vector<
                      kw::sde_rcomma::info::expect::type > >,
  tag::spike,      std::vector< std::vector< std::vector <
                      kw::spike::info::expect::type > > >,
  tag::betapdf,     std::vector< std::vector< std::vector <
                      kw::betapdf::info::expect::type > > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Mix mass-fraction beta parameters storage
using MixMassFractionBetaParameters = tk::tuple::tagged_tuple<
  tag::depvar,          std::vector< char >,
  tag::bprime,          std::vector< std::vector<
                          kw::sde_bprime::info::expect::type > >,
  tag::S,               std::vector< std::vector<
                          kw::sde_S::info::expect::type > >,
  tag::kappaprime,     std::vector< std::vector<
                          kw::sde_kappaprime::info::expect::type > >,
  tag::rho2,           std::vector< std::vector<
                          kw::sde_rho2::info::expect::type > >,
  tag::r,              std::vector< std::vector<
                          kw::sde_r::info::expect::type > >,
  tag::spike,          std::vector< std::vector< std::vector <
                          kw::spike::info::expect::type > > >,
  tag::betapdf,         std::vector< std::vector< std::vector <
                          kw::betapdf::info::expect::type > > >,
  tag::hydrotimescales, std::vector< std::vector< ctr::HydroTimeScalesType > >,
  tag::hydroproductions,std::vector< std::vector< ctr::HydroProductionsType > >,
  tag::rng,             std::vector< tk::ctr::RNGType >,
  tag::initpolicy,      std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy,     std::vector< ctr::CoeffPolicyType >
>;

//! Parameters storage
using parameters = tk::tuple::tagged_tuple<
  #ifdef HAS_MKL
  tag::rngmkl,          tk::ctr::RNGMKLParameters,
  #endif
  tag::rngsse,          tk::ctr::RNGSSEParameters,
  tag::rng123,          tk::ctr::RNGRandom123Parameters,
  tag::dirichlet,       DirichletParameters,
  tag::gendir,          GenDirichletParameters,
  tag::wrightfisher,    WrightFisherParameters,
  tag::ou,              OrnsteinUhlenbeckParameters,
  tag::diagou,          DiagOrnsteinUhlenbeckParameters,
  tag::skewnormal,      SkewNormalParameters,
  tag::gamma,           GammaParameters,
  tag::beta,            BetaParameters,
  tag::numfracbeta,     NumberFractionBetaParameters,
  tag::massfracbeta,    MassFractionBetaParameters,
  tag::mixnumfracbeta,  MixNumberFractionBetaParameters,
  tag::mixmassfracbeta, MixMassFractionBetaParameters
>;

} // ctr::
} // walker::

#endif // WalkerTypes_h
