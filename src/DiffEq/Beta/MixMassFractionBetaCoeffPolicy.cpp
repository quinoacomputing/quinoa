// *****************************************************************************
/*!
  \file      src/DiffEq/Beta/MixMassFractionBetaCoeffPolicy.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mass-fraction beta SDE coefficients policies
  \details   This file defines coefficients policy classes for the mass-fraction
             beta SDE, defined in DiffEq/Beta/MixMassFractionBeta.h. For general
             requirements on mixture mass-fraction beta SDE coefficients policy
             classes see the header file.
*/
// *****************************************************************************

#include <iostream>

#include "MixMassFractionBetaCoeffPolicy.hpp"
#include "TxtStatWriter.hpp"
#include "Walker/InputDeck/InputDeck.hpp"

namespace walker {

extern ctr::InputDeck g_inputdeck;

} // ::walker

walker::MixMassFracBetaCoeffDecay::MixMassFracBetaCoeffDecay(
  ncomp_t ncomp,
  const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
  const std::vector< kw::sde_S::info::expect::type >& S_,
  const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
  const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
  const std::vector< kw::sde_r::info::expect::type >& r_,
  std::vector< kw::sde_bprime::info::expect::type  >& bprime,
  std::vector< kw::sde_S::info::expect::type >& S,
  std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
  std::vector< kw::sde_rho2::info::expect::type >& rho2,
  std::vector< kw::sde_r::info::expect::type >& r,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_kappa::info::expect::type >& k )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] bprime_ Vector used to initialize coefficient vector bprime
//! \param[in] S_ Vector used to initialize coefficient vector S
//! \param[in] kprime_ Vector used to initialize coefficient vector kprime
//! \param[in] rho2_ Vector used to initialize coefficient vector rho2
//! \param[in] r_ Vector used to initialize coefficient vector r
//! \param[in,out] bprime Coefficient vector to be initialized
//! \param[in,out] S Coefficient vector to be initialized
//! \param[in,out] kprime Coefficient vector to be initialized
//! \param[in,out] rho2 Coefficient vector to be initialized
//! \param[in,out] r Coefficient vector to be initialized
//! \param[in,out] b Coefficient vector to be initialized
//! \param[in,out] k Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( bprime_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'b'");
  ErrChk( S_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'S'");
  ErrChk( kprime_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'kappa'");
  ErrChk( rho2_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'rho2'");
  ErrChk( r_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'r'");

  bprime = bprime_;
  S = S_;
  kprime = kprime_;
  rho2 = rho2_;
  r = r_;

  b.resize( bprime.size() );
  k.resize( kprime.size() );
}

void
walker::MixMassFracBetaCoeffDecay::update(
  char depvar,
  char,
  char,
  ctr::DepvarType,
  ctr::DepvarType /*scalar_solve*/,
  ncomp_t ncomp,
  const std::map< tk::ctr::Product, tk::real >& moments,
  const std::vector< kw::sde_bprime::info::expect::type  >& bprime,
  const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
  const std::vector< kw::sde_rho2::info::expect::type >&,
  const std::vector< kw::sde_r::info::expect::type >&,
  const std::vector< tk::Table >&,
  const std::vector< tk::Table >&,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_kappa::info::expect::type >& k,
  std::vector< kw::sde_S::info::expect::type >&,
  tk::real ) const
// *****************************************************************************
//  Update coefficients
//! \param[in] depvar Dependent variable
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] moments Map of statistical moments estimated
//! \param[in] bprime Coefficient vector b'
//! \param[in] kprime Coefficient vector kappa'
//! \param[in,out] b Coefficient vector to be updated
//! \param[in,out] k Coefficient vector to be updated
//! \details This where the mix mass-fraction beta SDE is made consistent
//!   with the no-mix and fully mixed limits by specifying the SDE
//!   coefficients, b and kappa as functions of b' and kappa'. We leave S
//!   unchanged.
// *****************************************************************************
{
  for (ncomp_t c=0; c<ncomp; ++c) {
    tk::real m = tk::ctr::lookup( tk::ctr::mean(depvar,c), moments );
    tk::real v = tk::ctr::lookup( tk::ctr::variance(depvar,c), moments );

    if (m<1.0e-8 || m>1.0-1.0e-8) m = 0.5;
    if (v<1.0e-8 || v>1.0-1.0e-8) v = 0.5;

    b[c] = bprime[c] * (1.0 - v / m / ( 1.0 - m ));
    k[c] = kprime[c] * v;
  }
}

walker::MixMassFracBetaCoeffHomDecay::MixMassFracBetaCoeffHomDecay(
  ncomp_t ncomp,
  const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
  const std::vector< kw::sde_S::info::expect::type >& S_,
  const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
  const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
  const std::vector< kw::sde_r::info::expect::type >& r_,
  std::vector< kw::sde_bprime::info::expect::type  >& bprime,
  std::vector< kw::sde_S::info::expect::type >& S,
  std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
  std::vector< kw::sde_rho2::info::expect::type >& rho2,
  std::vector< kw::sde_r::info::expect::type >& r,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_kappa::info::expect::type >& k )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] bprime_ Vector used to initialize coefficient vector bprime
//! \param[in] S_ Vector used to initialize coefficient vector S
//! \param[in] kprime_ Vector used to initialize coefficient vector kprime
//! \param[in] rho2_ Vector used to initialize coefficient vector rho2
//! \param[in] r_ Vector used to initialize coefficient vector r
//! \param[in,out] bprime Coefficient vector to be initialized
//! \param[in,out] S Coefficient vector to be initialized
//! \param[in,out] kprime Coefficient vector to be initialized
//! \param[in,out] rho2 Coefficient vector to be initialized
//! \param[in,out] r Coefficient vector to be initialized
//! \param[in,out] b Coefficient vector to be initialized
//! \param[in,out] k Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( bprime_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'b''");
  ErrChk( S_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'S'");
  ErrChk( kprime_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'kappa''");
  ErrChk( rho2_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'rho2'");
  ErrChk( r_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'r'");

  bprime = bprime_;
  S = S_;
  kprime = kprime_;
  rho2 = rho2_;
  r = r_;

  b.resize( bprime.size() );
  k.resize( kprime.size() );
}

void
walker::MixMassFracBetaCoeffHomDecay::update(
  char depvar,
  char,
  char,
  ctr::DepvarType,
  ctr::DepvarType /*scalar_solve*/,
  ncomp_t ncomp,
  const std::map< tk::ctr::Product, tk::real >& moments,
  const std::vector< kw::sde_bprime::info::expect::type  >& bprime,
  const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
  const std::vector< kw::sde_rho2::info::expect::type >& rho2,
  const std::vector< kw::sde_r::info::expect::type >& r,
  const std::vector< tk::Table >&,
  const std::vector< tk::Table >&,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_kappa::info::expect::type >& k,
  std::vector< kw::sde_S::info::expect::type >& S,
  tk::real ) const
// *****************************************************************************
//  Update coefficients
//! \param[in] depvar Dependent variable
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] moments Map of statistical moments estimated
//! \param[in] bprime Coefficient vector b'
//! \param[in] kprime Coefficient vector kappa'
//! \param[in] rho2 Coefficient vector rho2
//! \param[in] r Coefficient vector r
//! \param[in,out] b Coefficient vector to be updated
//! \param[in,out] k Coefficient vector to be updated
//! \param[in,out] S Coefficient vector to be updated
//! \details This where the mix mass-fraction beta SDE is made consistent
//!   with the no-mix and fully mixed limits by specifying the SDE
//!   coefficients, b and kappa as functions of b' and kappa'. We also
//!   specify S to force d\<rho\>/dt = 0, where \<rho\> = rho_2/(1+rY).
// *****************************************************************************
{
  using tk::ctr::lookup;
  using tk::ctr::mean;
  using tk::ctr::variance;
  using tk::ctr::cen3;

  // statistics nomenclature:
  //   Y = instantaneous mass fraction,
  //   R = instantaneous density,
  //   y = Y - <Y>, mass fraction fluctuation about its mean,
  //   r = R - <R>, density fluctuation about its mean,
  // <Y> = mean mass fraction,
  // <R> = mean density,

  for (ncomp_t c=0; c<ncomp; ++c) {
    tk::real m = lookup( mean(depvar,c), moments );            // <Y>
    tk::real v = lookup( variance(depvar,c), moments );        // <y^2>
    tk::real d = lookup( mean(depvar,c+ncomp), moments );      // <R>
    tk::real d2 = lookup( variance(depvar,c+ncomp), moments ); // <r^2>
    tk::real d3 = lookup( cen3(depvar,c+ncomp), moments );     // <r^3>

    if (m<1.0e-8 || m>1.0-1.0e-8) m = 0.5;
    if (v<1.0e-8 || v>1.0-1.0e-8) v = 0.5;
    b[c] = bprime[c] * (1.0 - v/m/(1.0-m));
    //b[c] = bprime[c] * (1.0 - v/M[c]/(1.0-M[c]));
    k[c] = kprime[c] * v;
    //b[c] = 1.0;
    //k[c] = 0.5*v/(m*(1.0-m));

    if (d < 1.0e-8) {
      std::cout << "d:" << d << " ";
      d = 0.5;
    }
    tk::real R = 1.0 + d2/d/d;
    tk::real B = -1.0/r[c]/r[c];
    tk::real C = (2.0+r[c])/r[c]/r[c];
    tk::real D = -(1.0+r[c])/r[c]/r[c];
    tk::real diff =
      B*d/rho2[c] +
      C*d*d*R/rho2[c]/rho2[c] +
      D*d*d*d*(1.0 + 3.0*d2/d/d + d3/d/d/d)/rho2[c]/rho2[c]/rho2[c];
    S[c] = (rho2[c]/d/R +
            2.0*k[c]/b[c]*rho2[c]*rho2[c]/d/d*r[c]*r[c]/R*diff - 1.0) / r[c];
    if (S[c] < 0.0 || S[c] > 1.0) {
      std::cout << S[c] << " ";
      S[c] = 0.5;
    }
  }
}

walker::MixMassFracBetaCoeffMonteCarloHomDecay::
MixMassFracBetaCoeffMonteCarloHomDecay(
   ncomp_t ncomp,
   const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
   const std::vector< kw::sde_S::info::expect::type >& S_,
   const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
   const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
   const std::vector< kw::sde_r::info::expect::type >& r_,
   std::vector< kw::sde_bprime::info::expect::type  >& bprime,
   std::vector< kw::sde_S::info::expect::type >& S,
   std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
   std::vector< kw::sde_rho2::info::expect::type >& rho2,
   std::vector< kw::sde_r::info::expect::type >& r,
   std::vector< kw::sde_b::info::expect::type  >& b,
   std::vector< kw::sde_kappa::info::expect::type >& k )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] bprime_ Vector used to initialize coefficient vector bprime
//! \param[in] S_ Vector used to initialize coefficient vector S
//! \param[in] kprime_ Vector used to initialize coefficient vector kprime
//! \param[in] rho2_ Vector used to initialize coefficient vector rho2
//! \param[in] r_ Vector used to initialize coefficient vector r
//! \param[in,out] bprime Coefficient vector to be initialized
//! \param[in,out] S Coefficient vector to be initialized
//! \param[in,out] kprime Coefficient vector to be initialized
//! \param[in,out] rho2 Coefficient vector to be initialized
//! \param[in,out] r Coefficient vector to be initialized
//! \param[in,out] b Coefficient vector to be initialized
//! \param[in,out] k Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( bprime_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'b''");
  ErrChk( S_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'S'");
  ErrChk( kprime_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'kappa''");
  ErrChk( rho2_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'rho2'");
  ErrChk( r_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'r'");

  bprime = bprime_;
  S = S_;
  kprime = kprime_;
  rho2 = rho2_;
  r = r_;

  b.resize( bprime.size() );
  k.resize( kprime.size() );
}

void
walker::MixMassFracBetaCoeffMonteCarloHomDecay::update(
  char depvar,
  char,
  char,
  ctr::DepvarType,
  ctr::DepvarType /*scalar_solve*/,
  ncomp_t ncomp,
  const std::map< tk::ctr::Product, tk::real >& moments,
  const std::vector< kw::sde_bprime::info::expect::type  >& bprime,
  const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
  const std::vector< kw::sde_rho2::info::expect::type >& rho2,
  const std::vector< kw::sde_r::info::expect::type >& r,
  const std::vector< tk::Table >&,
  const std::vector< tk::Table >&,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_kappa::info::expect::type >& k,
  std::vector< kw::sde_S::info::expect::type >& S,
  tk::real ) const
// *****************************************************************************
//  Update coefficients
//! \param[in] depvar Dependent variable
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] moments Map of statistical moments estimated
//! \param[in] bprime Coefficient vector b'
//! \param[in] kprime Coefficient vector kappa'
//! \param[in] rho2 Coefficient vector rho2
//! \param[in] r Coefficient vector r
//! \param[in,out] b Coefficient vector to be updated
//! \param[in,out] k Coefficient vector to be updated
//! \param[in,out] S Coefficient vector to be updated
//! \details This where the mix mass-fraction beta SDE is made consistent
//!   with the no-mix and fully mixed limits by specifying the SDE
//!   coefficients, b and kappa as functions of b' and kappa'. We also
//!   specify S to force d\<rho\>/dt = 0, where \<rho\> = rho_2/(1+rY).
// *****************************************************************************
{
  using tk::ctr::lookup;
  using tk::ctr::mean;
  using tk::ctr::variance;
  using tk::ctr::ord2;
  // statistics nomenclature:
  //   Y = instantaneous mass fraction,
  //   R = instantaneous density,
  //   y = Y - <Y>, mass fraction fluctuation about its mean,
  //   r = R - <R>, density fluctuation about its mean,
  // <Y> = mean mass fraction,
  // <R> = mean density,
  for (ncomp_t c=0; c<ncomp; ++c) {
    tk::real m = lookup( mean(depvar,c), moments );            // <Y>
    tk::real v = lookup( variance(depvar,c), moments );        // <y^2>
    tk::real r2 = lookup( ord2(depvar,c+ncomp), moments );     // <R^2>

    const tk::ctr::Term Y( static_cast<char>(std::toupper(depvar)),
                           c,
                           tk::ctr::Moment::ORDINARY );
    const tk::ctr::Term R( static_cast<char>(std::toupper(depvar)),
                           c+ncomp,
                           tk::ctr::Moment::ORDINARY );
    const tk::ctr::Term OneMinusY( static_cast<char>(std::toupper(depvar)),
                                   c+3*ncomp,
                                   tk::ctr::Moment::ORDINARY );

    const auto YR2 = tk::ctr::Product( { Y, R, R } );
    const auto Y1MYR3 = tk::ctr::Product( { Y, OneMinusY, R, R, R } );

    tk::real yr2 = lookup( YR2, moments );          // <RY^2>
    tk::real y1myr3 = lookup( Y1MYR3, moments );    // <Y(1-Y)R^3>

    if (m<1.0e-8 || m>1.0-1.0e-8) m = 0.5;
    if (v<1.0e-8 || v>1.0-1.0e-8) v = 0.5;
    b[c] = bprime[c] * (1.0 - v/m/(1.0-m));
    k[c] = kprime[c] * v;
    //b[c] = 1.0;
    //k[c] = 0.5*v/(m*(1.0-m));

    if (r2 < 1.0e-8) {
      std::cout << "r2:" << r2 << " ";
      r2 = 0.5;
    }
    S[c] = (yr2 + 2.0*k[c]/b[c]*r[c]/rho2[c]*y1myr3) / r2;
    if (S[c] < 0.0 || S[c] > 1.0) {
      std::cout << "S:" << S[c] << " ";
      S[c] = 0.5;
    }
  }
}

walker::MixMassFracBetaCoeffHydroTimeScale::
MixMassFracBetaCoeffHydroTimeScale(
  ncomp_t ncomp,
  const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
  const std::vector< kw::sde_S::info::expect::type >& S_,
  const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
  const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
  const std::vector< kw::sde_r::info::expect::type >& r_,
  std::vector< kw::sde_bprime::info::expect::type  >& bprime,
  std::vector< kw::sde_S::info::expect::type >& S,
  std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
  std::vector< kw::sde_rho2::info::expect::type >& rho2,
  std::vector< kw::sde_r::info::expect::type >& r,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_kappa::info::expect::type >& k )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] bprime_ Vector used to initialize coefficient vector bprime
//! \param[in] S_ Vector used to initialize coefficient vector S
//! \param[in] kprime_ Vector used to initialize coefficient vector kprime
//! \param[in] rho2_ Vector used to initialize coefficient vector rho2
//! \param[in] r_ Vector used to initialize coefficient vector r
//! \param[in,out] bprime Coefficient vector to be initialized
//! \param[in,out] S Coefficient vector to be initialized
//! \param[in,out] kprime Coefficient vector to be initialized
//! \param[in,out] rho2 Coefficient vector to be initialized
//! \param[in,out] r Coefficient vector to be initialized
//! \param[in,out] b Coefficient vector to be initialized
//! \param[in,out] k Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( bprime_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'b''");
  ErrChk( S_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'S'");
  ErrChk( kprime_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'kappa''");
  ErrChk( rho2_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'rho2'");
  ErrChk( r_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'r'");

  bprime = bprime_;
  S = S_;
  kprime = kprime_;
  rho2 = rho2_;
  r = r_;

  b.resize( bprime.size() );
  k.resize( kprime.size() );

  // Extra output besides normal statistics output
  m_extra_out_filename = "coeff";
  tk::TxtStatWriter sw( m_extra_out_filename );
  std::vector< std::string > names;
  for (ncomp_t c=0; c<ncomp; ++c) {
    names.push_back( "b" + std::to_string(c+1) );
    names.push_back( "S" + std::to_string(c+1) );
    names.push_back( "k" + std::to_string(c+1) );
  }
  sw.header( names, {}, {} );
}

void
walker::MixMassFracBetaCoeffHydroTimeScale::update(
  char depvar,
  char,
  char,
  ctr::DepvarType,
  ctr::DepvarType /*scalar_solve*/,
  ncomp_t ncomp,
  const std::map< tk::ctr::Product, tk::real >& moments,
  const std::vector< kw::sde_bprime::info::expect::type  >& bprime,
  const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
  const std::vector< kw::sde_rho2::info::expect::type >& rho2,
  const std::vector< kw::sde_r::info::expect::type >& r,
  const std::vector< tk::Table >& hts,
  const std::vector< tk::Table >& hp,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_kappa::info::expect::type >& k,
  std::vector< kw::sde_S::info::expect::type >& S,
  tk::real t ) const
// *****************************************************************************
//  Update coefficients
//! \param[in] depvar Dependent variable
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] moments Map of statistical moments estimated
//! \param[in] bprime Coefficient vector b'
//! \param[in] kprime Coefficient vector kappa'
//! \param[in] rho2 Coefficient vector rho2
//! \param[in] r Coefficient vector r
//! \param[in] hts (Inverse) hydrodynamics time scale (as input from DNS)
//! \param[in] hp Production/dissipation (as input from DNS)
//! \param[in,out] b Coefficient vector to be updated
//! \param[in,out] k Coefficient vector to be updated
//! \param[in,out] S Coefficient vector to be updated
//! \param[in] t Physical time at which to update coefficients
//! \details This where the mix mass-fraction beta SDE is made consistent
//!   with the no-mix and fully mixed limits by specifying the SDE
//!   coefficients, b and kappa as functions of b' and kappa'. Additionally,
//!   we pull in a hydrodynamic timescale from an external function. We also
//!   specify S to force d\<rho\>/dt = 0, where \<rho\> = rho_2/(1+rY).
// *****************************************************************************
{
  using tk::ctr::lookup;
  using tk::ctr::mean;
  using tk::ctr::variance;
  using tk::ctr::cen3;
  using tk::ctr::Product;

  if (m_it == 0)
    for (ncomp_t c=0; c<ncomp; ++c)
       m_s.push_back( S[c] );

  // Extra output besides normal statistics output
  tk::TxtStatWriter sw( m_extra_out_filename,
                        g_inputdeck.get< tag::flformat, tag::stat >(),
                        g_inputdeck.get< tag::prec, tag::stat >(),
                        std::ios_base::app );

  std::vector< tk::real > coeffs( ncomp * 3 );

  // statistics nomenclature:
  //   Y = instantaneous mass fraction,
  //   R = instantaneous density,
  //   y = Y - <Y>, mass fraction fluctuation about its mean,
  //   r = R - <R>, density fluctuation about its mean,
  // <Y> = mean mass fraction,
  // <R> = mean density,
  for (ncomp_t c=0; c<ncomp; ++c) {

    const tk::ctr::Term Y( static_cast<char>(std::toupper(depvar)),
                           c,
                           tk::ctr::Moment::ORDINARY );
    const tk::ctr::Term dens( static_cast<char>(std::toupper(depvar)),
                              c+ncomp,
                              tk::ctr::Moment::ORDINARY );
    const tk::ctr::Term s1( static_cast<char>(std::tolower(depvar)),
                            c+ncomp,
                            tk::ctr::Moment::CENTRAL );
    const tk::ctr::Term s2( static_cast<char>(std::tolower(depvar)),
                            c+ncomp*2,
                            tk::ctr::Moment::CENTRAL );

    const auto RY = tk::ctr::Product( { dens, Y } );
    tk::real ry = lookup( RY, moments );                       // <RY>
    const auto dscorr = tk::ctr::Product( { s1, s2 } );
    tk::real ds = -lookup( dscorr, moments );                  // b = -<rv>
    tk::real d = lookup( mean(depvar,c+ncomp), moments );      // <R>
    tk::real d2 = lookup( variance(depvar,c+ncomp), moments ); // <r^2>
    tk::real d3 = lookup( cen3(depvar,c+ncomp), moments );     // <r^3>
    tk::real yt = ry/d;

    // Sample hydrodynamics timescale and prod/diss at time t
    auto ts = hydrotimescale( t, hts[c] );  // eps/k
    auto pe = hydroproduction( t, hp[c] );  // P/eps = (dk/dt+eps)/eps

    tk::real a = r[c]/(1.0+r[c]*yt);
    tk::real bnm = a*a*yt*(1.0-yt);
    tk::real thetab = 1.0 - ds/bnm;
    tk::real f2 =
      1.0 / std::pow(1.0 + std::pow(pe-1.0,2.0)*std::pow(ds,0.25),0.5);
    tk::real b1 = m_s[0];
    tk::real b2 = m_s[1];
    tk::real b3 = m_s[2];
    tk::real eta = d2/d/d/ds;
    tk::real beta2 = b2*(1.0+eta*ds);
    tk::real Thetap = thetab*0.5*(1.0+eta/(1.0+eta*ds));
    tk::real beta3 = b3*(1.0+eta*ds);
    tk::real beta10 = b1 * (1.0+ds)/(1.0+eta*ds);
    tk::real beta1 = bprime[c] * 2.0/(1.0+eta+eta*ds) *
                  (beta10 + beta2*Thetap*f2 + beta3*Thetap*(1.0-Thetap)*f2);
    b[c] = beta1 * ts;
    k[c] = kprime[c] * beta1 * ts * ds * ds;
    //b[c] = bprime[c];
    //k[c] = kprime[c];
    //b[c] = bprime[c] + 0.25*std::sin(10.0*t);
    //k[c] = 1.0 + 0.25*std::sin(10.0*t);
    //k[c] = -(1.0 + std::sin(t)) * (S[c] - 1.0);

    tk::real R = 1.0 + d2/d/d;
    tk::real B = -1.0/r[c]/r[c];
    tk::real C = (2.0+r[c])/r[c]/r[c];
    tk::real D = -(1.0+r[c])/r[c]/r[c];
    tk::real diff =
      B*d/rho2[c] +
      C*d*d*R/rho2[c]/rho2[c] +
      D*d*d*d*(1.0 + 3.0*d2/d/d + d3/d/d/d)/rho2[c]/rho2[c]/rho2[c];
    S[c] = (rho2[c]/d/R +
           2.0*k[c]/b[c]*rho2[c]*rho2[c]/d/d*r[c]*r[c]/R*diff - 1.0) / r[c];
    //S[c] = 0.5 + 0.25*std::sin(10.0*t);

    // Implementation of a constraint for MixDirichlet for S_al to keep
    // d<rho>/dt = 0 to verify its behavior for beta (binary case). As input
    // file use mixmassfractbeta_mmS_A0.75.q.
    // auto R2 = lookup( Product({dens,dens}), moments ); // <R^2>
    // auto R2Yc = lookup( Product({dens,dens,Y}), moments ); // <R^2Yc>
    // auto R3Yc = lookup( Product({dens,dens,dens,Y}), moments ); // <R^3Yc>
    // auto R3Y2c = lookup( Product({dens,dens,dens,Y,Y}), moments ); // <R^3Y^2>
    // tk::real drYc = -r[c]/rho2[c]*R2;
    // tk::real drYcYc = -r[c]/rho2[c]*R2Yc;
    // tk::real drYc2YcYN = 2.0*std::pow(r[c]/rho2[c],2.0)*(R3Yc-R3Y2c);
    // S[c] = (drYcYc - k[c]/b[c]*drYc2YcYN) / drYc;

    coeffs[ 3*c+0 ] = b[c];
    coeffs[ 3*c+1 ] = S[c];
    coeffs[ 3*c+2 ] = k[c];
  }

  // Extra "stat" output of coefficients
  sw.stat( 0, t, coeffs, {}, {} );

  ++m_it;
}

walker::MixMassFracBetaCoeffInstVel::MixMassFracBetaCoeffInstVel(
  ncomp_t ncomp,
  const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
  const std::vector< kw::sde_S::info::expect::type >& S_,
  const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
  const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
  const std::vector< kw::sde_r::info::expect::type >& r_,
  std::vector< kw::sde_bprime::info::expect::type  >& bprime,
  std::vector< kw::sde_S::info::expect::type >& S,
  std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
  std::vector< kw::sde_rho2::info::expect::type >& rho2,
  std::vector< kw::sde_r::info::expect::type >& r,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_kappa::info::expect::type >& k )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] bprime_ Vector used to initialize coefficient vector bprime
//! \param[in] S_ Vector used to initialize coefficient vector S
//! \param[in] kprime_ Vector used to initialize coefficient vector kprime
//! \param[in] rho2_ Vector used to initialize coefficient vector rho2
//! \param[in] r_ Vector used to initialize coefficient vector r
//! \param[in,out] bprime Coefficient vector to be initialized
//! \param[in,out] S Coefficient vector to be initialized
//! \param[in,out] kprime Coefficient vector to be initialized
//! \param[in,out] rho2 Coefficient vector to be initialized
//! \param[in,out] r Coefficient vector to be initialized
//! \param[in,out] b Coefficient vector to be initialized
//! \param[in,out] k Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( bprime_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'b''");
  ErrChk( S_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'S'");
  ErrChk( kprime_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'kappa''");
  ErrChk( rho2_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'rho2'");
  ErrChk( r_.size() == ncomp,
    "Wrong number of mix mass-fraction beta SDE parameters 'r'");

  bprime = bprime_;
  S = S_;
  kprime = kprime_;
  rho2 = rho2_;
  r = r_;

  b.resize( bprime.size() );
  k.resize( kprime.size() );
}

void
walker::MixMassFracBetaCoeffInstVel::update(
  char depvar,
  char dissipation_depvar,
  char /*velocity_depvar*/,
  ctr::DepvarType /*velocity_solve*/,
  ctr::DepvarType solve,
  ncomp_t ncomp,
  const std::map< tk::ctr::Product, tk::real >& moments,
  const std::vector< kw::sde_bprime::info::expect::type  >& /*bprime*/,
  const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
  const std::vector< kw::sde_rho2::info::expect::type >& /*rho2*/,
  const std::vector< kw::sde_r::info::expect::type >& /*r*/,
  const std::vector< tk::Table >&,
  const std::vector< tk::Table >&,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_kappa::info::expect::type >& k,
  std::vector< kw::sde_S::info::expect::type >& S,
  tk::real ) const
// *****************************************************************************
//  Update coefficients
//! \param[in] depvar Dependent variable
//! \param[in] dissipation_depvar Dependent variable of coupled dissipation eq
//! \param[in] solve Enum selecting whether the full variable or its
//!   fluctuation is solved for
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] moments Map of statistical moments estimated
//! \param[in] kprime Coefficient vector kappa'
//! \param[in] rho2 Coefficient vector rho2
//! \param[in] r Coefficient vector r
//! \param[in,out] b Coefficient vector to be updated
//! \param[in,out] k Coefficient vector to be updated
//! \param[in,out] S Coefficient vector to be updated
//! \details This where the mix mass-fraction beta SDE is made consistent
//!   with the no-mix and fully mixed limits by specifying the SDE
//!   coefficients, b and kappa as functions of b' and kappa'. Additionally,
//!   we specify the hydrodynamic timescale by coupling to anothe SDE, computing
//!   the velocity field. We also specify S to force d\<rho\>/dt = 0, where
//!   \<rho\> = rho_2/(1+rY).
// *****************************************************************************
{
  using tk::ctr::lookup;
  using tk::ctr::mean;
  using tk::ctr::variance;
  using tk::ctr::cen3;

  if (m_it == 0) {
    for (ncomp_t c=0; c<ncomp; ++c) {
       m_s.push_back( S[c] );
    }
  }

  for (ncomp_t c=0; c<ncomp; ++c) {

    tk::real y2 = 0.0;
    tk::real ts = 1.0;

    if (solve == ctr::DepvarType::FULLVAR) {

      y2 = lookup( variance(depvar,c), moments );       // <y^2>

      // Access mean turbulence frequency from coupled dissipation model
      // hydroptimescale: eps/k = <O>
      if (dissipation_depvar != '-') {     // only if dissipation is coupled
        ts = lookup( mean(dissipation_depvar,0), moments );
      }

    } else if (solve == ctr::DepvarType::FLUCTUATION) {

      // Since we are solving for the fluctuating scalar, the "ordinary"
      // moments are really central moments
      auto d = static_cast< char >( std::toupper( depvar ) );
      tk::ctr::Term y( d, c, tk::ctr::Moment::ORDINARY );
      y2 = lookup( tk::ctr::Product( {y,y} ), moments );       // <y^2>

      // Access mean turbulence frequency from coupled dissipation model
      // hydroptimescale: eps/k = <O>
      if (dissipation_depvar != '-') {     // only if dissipation is coupled
        ts = lookup( mean(dissipation_depvar,0), moments );
      }

    } else Throw( "Depvar type not implemented" );

    tk::real beta1 = 2.0;
    b[c] = beta1 * ts;
    k[c] = kprime[c] * beta1 * ts * y2;

    S[c] = 0.5;
  }  

  ++m_it;
}
