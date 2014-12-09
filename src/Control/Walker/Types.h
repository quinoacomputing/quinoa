//******************************************************************************
/*!
  \file      src/Control/Walker/Types.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 09:34:35 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Types for Walker's parsers
  \details   Types for Walker's parsers
*/
//******************************************************************************
#ifndef WalkerTypes_h
#define WalkerTypes_h

#include <Tags.h>
#include <Types.h>
#include <ControlTypes.h>
#include <Walker/Options/DiffEq.h>
#include <Options/InitPolicy.h>
#include <Options/CoeffPolicy.h>
#include <Options/PDFFile.h>
#include <Options/PDFPolicy.h>
#include <Options/PDFCentering.h>
#include <Options/TxtFloatFormat.h>
#include <Options/RNG.h>

namespace walker {
//! control and parsing
namespace ctr {

//! Storage of selected options
using selects = tk::tuple::tagged_tuple<
  tag::diffeq,       std::vector< ctr::DiffEqType >,  //!< Differential eqs
  tag::rng,          std::vector< tk::ctr::RNGType >, //!< RNGs
  tag::pdffiletype,  tk::ctr::PDFFileType,      //!< PDF output file type
  tag::pdfpolicy,    tk::ctr::PDFPolicyType,    //!< PDF output file policy
  tag::pdfctr,       tk::ctr::PDFCenteringType, //!< PDF output file centering
  tag::float_format, tk::ctr::TxtFloatFormatType//!< Text floating-point format
>;

//! Discretization parameters storage
using discretization = tk::tuple::tagged_tuple<
  tag::npar,      uint64_t,  //!< Total number of particles
  tag::nstep,     uint64_t,  //!< Number of time steps to take
  tag::term,      tk::real,  //!< Time to terminate time stepping
  tag::dt,        tk::real,  //!< Size of time step
  tag::binsize,   std::vector< std::vector< tk::real > >, //!< PDF binsizes
  tag::extent,    std::vector< std::vector< tk::real > >, //!< PDF extents
  tag::precision, std::streamsize  //!< Precision in digits
>;

//! Output intervals storage
using intervals = tk::tuple::tagged_tuple<
  tag::tty,  uint32_t,  //!< TTY output interval
  tag::stat, uint32_t,  //!< Statistics output interval
  tag::pdf,  uint32_t   //!< PDF output interval
>;

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::control,         std::string,  //!< Control filename
  tag::input,           std::string,  //!< Input filename
  tag::output,          std::string,  //!< Output filename
  tag::pdf,             std::string,  //!< PDF filename
  tag::stat,            std::string,  //!< Statistics filename
  tag::pdfnames,        std::vector< std::string >  //!< PDF identifiers
>;

//! Dirichlet parameters storage
using DirichletParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector< tk::real > >,
  tag::S,           std::vector< std::vector< tk::real > >,
  tag::kappa,       std::vector< std::vector< tk::real > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< tk::ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< tk::ctr::CoeffPolicyType >
>;

//! Generalized Dirichlet parameters storage
using GenDirichletParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector< tk::real > >,
  tag::S,           std::vector< std::vector< tk::real > >,
  tag::kappa,       std::vector< std::vector< tk::real > >,
  tag::c,           std::vector< std::vector< tk::real > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< tk::ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< tk::ctr::CoeffPolicyType >
>;

//! Wright-Fisher parameters storage
using WrightFisherParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::omega,       std::vector< std::vector< tk::real > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< tk::ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< tk::ctr::CoeffPolicyType >
>;

//! Ornstein-Uhlenbeck parameters storage
using OrnsteinUhlenbeckParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::sigma,       std::vector< std::vector< tk::real > >,
  tag::theta,       std::vector< std::vector< tk::real > >,
  tag::mu,          std::vector< std::vector< tk::real > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< tk::ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< tk::ctr::CoeffPolicyType >
>;

//! Diagonal Ornstein-Uhlenbeck parameters storage
using DiagOrnsteinUhlenbeckParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::sigma,       std::vector< std::vector< tk::real > >,
  tag::theta,       std::vector< std::vector< tk::real > >,
  tag::mu,          std::vector< std::vector< tk::real > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< tk::ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< tk::ctr::CoeffPolicyType >
>;

//! Skew-normal parameters storage
using SkewNormalParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::timescale,   std::vector< std::vector< tk::real > >,
  tag::sigma,       std::vector< std::vector< tk::real > >,
  tag::lambda,      std::vector< std::vector< tk::real > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< tk::ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< tk::ctr::CoeffPolicyType >
>;

//! Gamma parameters storage
using GammaParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector< tk::real > >,
  tag::S,           std::vector< std::vector< tk::real > >,
  tag::kappa,       std::vector< std::vector< tk::real > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< tk::ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< tk::ctr::CoeffPolicyType >
>;

//! Beta parameters storage
using BetaParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector< tk::real > >,
  tag::S,           std::vector< std::vector< tk::real > >,
  tag::kappa,       std::vector< std::vector< tk::real > >,
  tag::rng,         std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< tk::ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< tk::ctr::CoeffPolicyType >
>;

//! Parameters storage
using parameters = tk::tuple::tagged_tuple<
  #ifdef HAS_MKL
  tag::rngmkl,       tk::ctr::RNGMKLParameters,   //!< MKL RNG parameters
  #endif
  tag::rngsse,       tk::ctr::RNGSSEParameters,   //!< RNGSSE RNG parameters
  tag::dirichlet,    DirichletParameters,
  tag::gendir,       GenDirichletParameters,
  tag::wrightfisher, WrightFisherParameters,
  tag::ou,           OrnsteinUhlenbeckParameters,
  tag::diagou,       DiagOrnsteinUhlenbeckParameters,
  tag::skewnormal,   SkewNormalParameters,
  tag::gamma,        GammaParameters,
  tag::beta,         BetaParameters
>;

//! PEGTL location type to use throughout Walker's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // walker::

#endif // WalkerTypes_h
