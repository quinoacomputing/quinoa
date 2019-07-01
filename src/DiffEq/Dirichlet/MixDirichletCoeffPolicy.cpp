// *****************************************************************************
/*!
  \file      src/DiffEq/Dirichlet/MixDirichletCoeffPolicy.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mixture Dirichlet coefficients policies
  \details   This file defines coefficients policy classes for the mixture
             Dirichlet SDE, defined in DiffEq/Dirichlet/MixDirichlet.h.
             For general requirements on the mixture Dirichlet SDE coefficients
             policy classes see the header file.
*/
// *****************************************************************************

#include <iostream>

#include "MixDirichletCoeffPolicy.hpp"

std::vector< kw::sde_r::info::expect::type >
walker::MixDir_r( const std::vector< kw::sde_rho::info::expect::type >& rho,
                  ctr::NormalizationType norm )
// *****************************************************************************
//  Compute parameter vector r based on r_i = rho_N/rho_i - 1
//! \param[in] rho Parameter vector rho to MixDirichlet
//! \param[in] norm Normalization type (N=heavy or N=light)
//! \return Parameter vector r, determined by parameter vector rho
// *****************************************************************************
{
  Assert( rho.size() > 1, "Parameter vector rho must not be empty" );

  std::vector< kw::sde_r::info::expect::type > r( rho.size()-1 );

  for (std::size_t i=0; i<rho.size()-1; ++i) {
    if (norm == ctr::NormalizationType::LIGHT)
      r[i] = rho.back()/rho[i] + 1.0;
    else
      r[i] = rho.back()/rho[i] - 1.0;
  }

  //if (norm == ctr::NormalizationType::LIGHT)
  //  r.push_back( 2.0 );

  return r;
}

walker::MixDirichletCoeffConst::MixDirichletCoeffConst(
  tk::ctr::ncomp_type ncomp,
  ctr::NormalizationType norm,
  const std::vector< kw::sde_b::info::expect::type >& b_,
  const std::vector< kw::sde_S::info::expect::type >& S_,
  const std::vector< kw::sde_kappa::info::expect::type >& kprime_,
  const std::vector< kw::sde_rho::info::expect::type >& rho_,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_S::info::expect::type >& S,
  std::vector< kw::sde_kappa::info::expect::type >& kprime,
  std::vector< kw::sde_rho::info::expect::type >& rho,
  std::vector< kw::sde_r::info::expect::type >& r,
  std::vector< kw::sde_kappa::info::expect::type >& k )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] b_ Vector used to initialize coefficient vector b
//! \param[in] S_ Vector used to initialize coefficient vector S
//! \param[in] kprime_ Vector used to initialize coefficient vector kprime and k
//! \param[in] rho_ Vector used to initialize coefficient vector rho and r
//! \param[in,out] b Coefficient vector to be initialized
//! \param[in,out] S Coefficient vector to be initialized
//! \param[in,out] kprime Coefficient vector to be initialized
//! \param[in,out] rho Coefficient vector to be initialized
//! \param[in,out] r Coefficient vector to be initialized
//! \param[in,out] k Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( b_.size() == ncomp,
          "Wrong number of MixDirichlet SDE parameters 'b'");
  ErrChk( S_.size() == ncomp,
          "Wrong number of MixDirichlet SDE parameters 'S'");
  ErrChk( kprime_.size() == ncomp,
          "Wrong number of MixDirichlet SDE parameters 'kappaprime'");
  ErrChk( rho_.size() == ncomp+1,
          "Wrong number of MixDirichlet SDE parameters 'rho'");

  b = b_;
  S = S_;
  kprime = kprime_;
  rho = rho_;
  k.resize( kprime.size(), 0.0 );

  // Compute parameter vector r based on r_i = rho_N/rho_i - 1
  Assert( r.empty(), "Parameter vector r must be empty" );
  r = MixDir_r( rho, norm );
}

void
walker::MixDirichletCoeffConst::update(
  char /*depvar*/,
  ncomp_t ncomp,
  std::size_t /* density_offset */,
  std::size_t /* volume_offset */,
  const std::map< tk::ctr::Product, tk::real >& /*moments*/,
  const std::vector< kw::sde_rho::info::expect::type >& /*rho*/,
  const std::vector< kw::sde_r::info::expect::type >& /*r*/,
  const std::vector< kw::sde_kappa::info::expect::type >& kprime,
  const std::vector< kw::sde_b::info::expect::type >& /*b*/,
  std::vector< kw::sde_kappa::info::expect::type >& k,
  std::vector< kw::sde_kappa::info::expect::type >& S ) const
// *****************************************************************************
//  Update coefficients
//! \param[in] depvar Dependent variable
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] density_offset Offset of particle density in solution array
//!    relative to YN
//! \param[in] volume_offset Offset of particle specific volume in solution
//!    array relative to YN
//! \param[in] moments Map of statistical moments estimated
//! \param[in] rho Coefficient vector
//! \param[in] r Coefficient Vector
//! \param[in] kprime Coefficient vector
//! \param[in] b Coefficient vector
//! \param[in,out] k Coefficient vector to be updated
//! \param[in,out] S Coefficient vector to be updated
// *****************************************************************************
{
  using tk::ctr::lookup;
  using tk::ctr::mean;
  using tk::ctr::variance;
  using tk::ctr::Term;
  using tk::ctr::Moment;
  using tk::ctr::Product;

  for (ncomp_t c=0; c<ncomp; ++c) {
    k[c] = kprime[c];
  }

  for (ncomp_t c=0; c<ncomp; ++c) {
    if (S[c] < 0.0 || S[c] > 1.0) {
      std::cout << "S[" << c << "] bounds violated: " << S[c] << '\n';
    }
  }
}

walker::MixDirichletHomogeneous::MixDirichletHomogeneous(
  tk::ctr::ncomp_type ncomp,
  ctr::NormalizationType norm,
  const std::vector< kw::sde_b::info::expect::type >& b_,
  const std::vector< kw::sde_S::info::expect::type >& S_,
  const std::vector< kw::sde_kappa::info::expect::type >& kprime_,
  const std::vector< kw::sde_rho::info::expect::type >& rho_,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_S::info::expect::type >& S,
  std::vector< kw::sde_kappa::info::expect::type >& kprime,
  std::vector< kw::sde_rho::info::expect::type >& rho,
  std::vector< kw::sde_r::info::expect::type >& r,
  std::vector< kw::sde_kappa::info::expect::type >& k )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] b_ Vector used to initialize coefficient vector b
//! \param[in] S_ Vector used to initialize coefficient vector S
//! \param[in] kprime_ Vector used to initialize coefficient vector kprime and k
//! \param[in] rho_ Vector used to initialize coefficient vector rho and r
//! \param[in,out] b Coefficient vector to be initialized
//! \param[in,out] S Coefficient vector to be initialized
//! \param[in,out] kprime Coefficient vector to be initialized
//! \param[in,out] rho Coefficient vector to be initialized
//! \param[in,out] r Coefficient vector to be initialized
//! \param[in,out] k Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( b_.size() == ncomp,
          "Wrong number of MixDirichlet SDE parameters 'b'");
  ErrChk( S_.size() == ncomp,
          "Wrong number of MixDirichlet SDE parameters 'S'");
  ErrChk( kprime_.size() == ncomp,
          "Wrong number of MixDirichlet SDE parameters 'kappaprime'");
  ErrChk( rho_.size() == ncomp+1,
          "Wrong number of MixDirichlet SDE parameters 'rho'");

  b = b_;
  S = S_;
  kprime = kprime_;
  rho = rho_;
  k.resize( kprime.size(), 0.0 );

  // Compute parameter vector r based on r_i = rho_N/rho_i - 1
  Assert( r.empty(), "Parameter vector r must be empty" );
  r = MixDir_r( rho, norm );
}

void
walker::MixDirichletHomogeneous::update(
  char depvar,
  ncomp_t ncomp,
  std::size_t density_offset,
  std::size_t volume_offset,
  const std::map< tk::ctr::Product, tk::real >& moments,
  const std::vector< kw::sde_rho::info::expect::type >& rho,
  const std::vector< kw::sde_r::info::expect::type >& r,
  const std::vector< kw::sde_kappa::info::expect::type >& kprime,
  const std::vector< kw::sde_b::info::expect::type >& b,
  std::vector< kw::sde_kappa::info::expect::type >& k,
  std::vector< kw::sde_kappa::info::expect::type >& S ) const
// *****************************************************************************
//  Update coefficients
//! \param[in] depvar Dependent variable
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] density_offset Offset of particle density in solution array
//!    relative to YN
//! \param[in] volume_offset Offset of particle specific volume in solution
//!    array relative to YN
//! \param[in] moments Map of statistical moments estimated
//! \param[in] rho Coefficient vector
//! \param[in] r Coefficient Vector
//! \param[in] kprime Coefficient vector
//! \param[in] b Coefficient vector
//! \param[in,out] k Coefficient vector to be updated
//! \param[in,out] S Coefficient vector to be updated
// *****************************************************************************
{
  using tk::ctr::lookup;
  using tk::ctr::mean;
  using tk::ctr::variance;
  using tk::ctr::Term;
  using tk::ctr::Moment;
  using tk::ctr::Product;

  // note: ncomp = K = N-1

  // constraint: r_i = rho_N/rho_c - 1

  // statistics nomenclature:
  //   Y = instantaneous mass fraction
  //   R = instantaneous density
  //   y = Y - <Y>, mass fraction fluctuation about its mean
  // <Y> = mean mass fraction
  // <R> = mean density

  // <R>
  tk::real R = lookup( mean(depvar,ncomp+density_offset), moments );
  //if (R < 1.0e-8) R = 1.0;

  // b = -<rv>, density-specific-volume covariance
  Term rhoprime( static_cast<char>(std::tolower(depvar)),
                 ncomp+density_offset, Moment::CENTRAL );
  Term vprime( static_cast<char>(std::tolower(depvar)),
               ncomp+volume_offset, Moment::CENTRAL );
  //auto ds = -lookup( Product({rhoprime,vprime}), moments );

  // b. = -<ry.>/<R>
  std::vector< tk::real > bc( ncomp, 0.0 );
  for (ncomp_t c=0; c<ncomp; ++c) {
    Term tr( static_cast<char>(std::tolower(depvar)),
             ncomp+density_offset, Moment::CENTRAL );
    Term ty( static_cast<char>(std::tolower(depvar)),
             c, Moment::CENTRAL );
    bc[c] = -lookup( Product({tr,ty}), moments ) / R; // -<ryc>/<R>
  }

  Term tR( static_cast<char>(std::toupper(depvar)), ncomp+density_offset,
           Moment::ORDINARY );
  Term tYN( static_cast<char>(std::toupper(depvar)), ncomp, Moment::ORDINARY );

  auto R2YN = lookup( Product({tR,tR,tYN}), moments );

  std::vector< tk::real > y2( ncomp, 0.0 );
  std::vector< tk::real > RY( ncomp, 0.0 );
  std::vector< tk::real > R2Y( ncomp, 0.0 );
  std::vector< tk::real > R3YNY( ncomp, 0.0 );
  std::vector< tk::real > R3Y2( ncomp*ncomp, 0.0 );
  for (ncomp_t c=0; c<ncomp; ++c) {
    y2[c] = lookup(
      variance(static_cast<char>(std::tolower(depvar)),c), moments );
    Term tYc( static_cast<char>(std::toupper(depvar)), c, Moment::ORDINARY );
    RY[c] = lookup( Product({tR,tYc}), moments );    // <RYc>
    R2Y[c] = lookup( Product({tR,tR,tYc}), moments ); // <R^2Yc>
    R3YNY[c] = lookup( Product({tR,tR,tR,tYN,tYc}), moments ); // <R^3YNYc>
    for (ncomp_t d=0; d<ncomp; ++d) {
      Term tYd( static_cast<char>(std::toupper(depvar)), d, Moment::ORDINARY );
      // <R^3YcYd>
      R3Y2[c*ncomp+d] = lookup( Product({tR,tR,tR,tYc,tYd}), moments );
    }
    //std::cout << "R2Y: " << R2Y[c] << ' ';
  }
  //std::cout << std::endl;

  // Reynolds means

  // Reynolds means, Yc
  std::vector< tk::real > Y( ncomp, 0.0 );
  for (ncomp_t c=0; c<ncomp; ++c) {
    Y[c] = lookup( mean(depvar,c), moments );
    //std::cout << "Y: " << Y[c] << ' ';
  }
  //std::cout << std::endl;

  // sum of Yc
  tk::real sumY = 0.0;
  for (ncomp_t c=0; c<ncomp; ++c) sumY += Y[c];

//      // Y|Kc
//      std::vector< tk::real > YK( ncomp, 0.0 );
//      for (ncomp_t c=0; c<ncomp; ++c) {
//        YK[c] = sumY - lookup( mean(depvar,c), moments );
//        //std::cout << "YK: " << YK[c] << ' ';
//      }
//      //std::cout << std::endl;

  // Favre means

  // Ytc
  std::vector< tk::real > Yt( ncomp, 0.0 );
  for (ncomp_t c=0; c<ncomp; ++c) {
    Yt[c] = RY[c] / R;
    //std::cout << "Yt: " << Yt[c] << ' ';
  }
  //std::cout << std::endl;

  // sum of Ytc
  tk::real sumYt = 0.0;
  for (ncomp_t c=0; c<ncomp; ++c) sumYt += Yt[c];
  //std::cout << "sumYt: " << sumYt << '\n';

//      // Yt|Kc
//      std::vector< tk::real > YtK( ncomp, 0.0 );
//      for (ncomp_t c=0; c<ncomp; ++c) {
//        YtK[c] = sumYt - Yt[c];
//        //std::cout << "YtK: " << YtK[c] << ' ';
//      }
//      //std::cout << std::endl;

  // sum of <R^2Yc>
  tk::real sumR2Y = 0.0;
  std::vector< tk::real > sumR3Y2( ncomp, 0.0 );
  for (ncomp_t c=0; c<ncomp; ++c) {
    sumR2Y += R2Y[c];
    for (ncomp_t d=0; d<ncomp; ++d) sumR3Y2[c] += R3Y2[c*ncomp+d];
  }
  //std::cout << "sumR2Y: " << sumR2Y << std::endl;

  // <r^2>, density variance
  //auto rhovar = lookup(
  //  variance(static_cast<char>(std::tolower(depvar)),ncomp), moments );
  //std::cout << "<r^2>: " << rhovar << std::endl;

  // <R^2>
  //auto R2 = lookup( Product({tR,tR}), moments );

  for (ncomp_t c=0; c<ncomp; ++c) {

    //k[c] = kprime[c] * bc[c];
    //k[c] = kprime[c] * ds;
    k[c] = kprime[c];
    //k[c] = kprime[c] * y2[c];
    //if (k[c] < 0.0)
    // std::cout << "Positivity of k[" << c << "] violated: "
    //           << k[c] << '\n';

    //S[c] = 1.0/(1.0-YK[c]) - (1.0-Yt[c])/(1.0-YtK[c]);
    //S[c] = YK[c]/(1.0-YK[c]) - (1.0-Yt[c])*YtK[c]/(1.0-YtK[c]) + Yt[c];
    //S[c] = Yt[c] / ( 1.0 - sumYt + Yt[c] );
    //S[c] = ( -2.0*(r[c]/rho[ncomp]*R2Y[c])*(1.0-sumYt) +
    //         (r[c]/rho[ncomp]*(rhovar-sumR2Y))*Yt[c] ) /
    //       ( -2.0*(r[c]/rho[ncomp]*R2Y[c])*(1.0-sumYt) -
    //         (1.0-sumYt-Yt[c])*(r[c]/rho[ncomp]*(rhovar-sumR2Y)) );

    // correlation of density gradient wrt Y_alpha and Y_alpha (Y_alpha = Yc)
    //tk::real drYcYc = -r[c]/rho[ncomp]*R2Y[c];
    //tk::real drYcYN = -r[c]/rho[ncomp]*(R2-sumR2Y);
    //tk::real drYc2YcYN = 2.0*std::pow(r[c]/rho[ncomp],2.0)*(R3Y[c]-sumR3Y2[c]);
    //S[c] = (drYcYc - k[c]/b[c]*drYc2YcYN) / (drYcYN + drYcYc);
    //S[c] = Yt[c] / (1.0 - sumYt + Yt[c]) - k[c]/b[c]*drYc2YcYN / (drYcYN + drYcYc);

    S[c] = (R2Y[c] + 2.0*r[c]/rho[ncomp]*k[c]/b[c]*R3YNY[c]) / (R2Y[c] + R2YN);
    //S[c] = Yt[c] / (1.0 - sumYt + Yt[c]);

    //std::cout << "S[" << c << "] = " << S[c] << ", using " << b[c] << ", "
    //          << drYc2YcYN << ", " << drYcYN << ", " << drYcYc << '\n';

  }
  //std::cout << std::endl;

  for (ncomp_t c=0; c<ncomp; ++c) {
    if (S[c] < 0.0 || S[c] > 1.0) {
      std::cout << "S[" << c << "] bounds violated: " << S[c] << '\n';
      //S[c] = 0.5;
    }
  }
  //std::cout << std::endl;
}
