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

#include "MixDirichletCoeffPolicy.hpp"

std::vector< kw::sde_r::info::expect::type >
walker::MixDir_r( const std::vector< kw::sde_rho::info::expect::type >& rho )
// *****************************************************************************
//  Compute parameter vector r based on r_i = rho_N/rho_i - 1
//! \param[in] rho Parameter vector rho to MixDirichlet
//! \return Parameter vector r, determined by parameter vector rho
// *****************************************************************************
{
  Assert( rho.size() > 1, "Parameter vector rho must not be empty" );

  std::vector< kw::sde_r::info::expect::type > r( rho.size()-1 );

  for (std::size_t i=0; i<rho.size()-1; ++i) r[i] = rho.back()/rho[i] - 1.0;

  return r;
}

walker::MixDirichletHomCoeffConst::MixDirichletHomCoeffConst(
  tk::ctr::ncomp_t ncomp,
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
  k.resize( kprime.size() );

  // Compute parameter vector r based on r_i = rho_N/rho_i - 1
  Assert( r.empty(), "Parameter vector r must be empty" );
  r = MixDir_r( rho );
}

void
walker::MixDirichletHomCoeffConst::update(
  char depvar,
  ncomp_t ncomp,
  const std::map< tk::ctr::Product, tk::real >& moments,
  const std::vector< kw::sde_rho::info::expect::type >& rho,
  const std::vector< kw::sde_r::info::expect::type >& r,
  const std::vector< kw::sde_kappa::info::expect::type >& kprime,
  std::vector< kw::sde_kappa::info::expect::type >& k,
  std::vector< kw::sde_kappa::info::expect::type >& S ) const
// *****************************************************************************
//  Update coefficients
//! \param[in] depvar Dependent variable
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] moments Map of statistical moments estimated
//! \param[in] rho Coefficient vector
//! \param[in] r Coefficient Vector
//! \param[in] kprime Coefficient vector
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
  //   r = R - <R>, density fluctuation about its mean
  // <Y> = mean mass fraction
  // <R> = mean density

  // <R>
  tk::real R = lookup( mean(depvar,ncomp), moments );
  if (R < 1.0e-8) R = 1.0;

  // b = -<rv>, density-specific-volume covariance
  // Term rhoprime( static_cast<char>(std::tolower(depvar)),
  //                ncomp, Moment::CENTRAL );
  // Term vprime( static_cast<char>(std::tolower(depvar)),
  //              ncomp+1, Moment::CENTRAL );
  // auto ds = -lookup( Product({rhoprime,vprime}), moments );

  // b. = -<ry.>/<R>
  std::vector< tk::real > bc( ncomp, 0.0 );
  for (ncomp_t c=0; c<ncomp; ++c) {
    Term tr( static_cast<char>(std::tolower(depvar)),
             ncomp, Moment::CENTRAL );
    Term ty( static_cast<char>(std::tolower(depvar)),
             c, Moment::CENTRAL );
    bc[c] = -lookup( Product({tr,ty}), moments ) / R; // -<ryc>/<R>
    //std::cout << "RRY: " << RRY[c] << ' ';
  }

  std::vector< tk::real > RY( ncomp, 0.0 );
  std::vector< tk::real > RRY( ncomp, 0.0 );
  for (ncomp_t c=0; c<ncomp; ++c) {
    Term tR( static_cast<char>(std::toupper(depvar)),
             ncomp, Moment::ORDINARY );
    Term tY( static_cast<char>(std::toupper(depvar)),
             c, Moment::ORDINARY );
    RY[c] = lookup( Product({tR,tY}), moments );    // <RYc>
    RRY[c] = lookup( Product({tR,tR,tY}), moments ); // <R^2Yc>
    //std::cout << "RRY: " << RRY[c] << ' ';
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
  tk::real sumRRY = 0.0;
  for (ncomp_t c=0; c<ncomp; ++c) sumRRY += RRY[c];
  //std::cout << "sumRRY: " << sumRRY << std::endl;

  // <r^2>, density variance
  auto rhovar = lookup(
    variance(static_cast<char>(std::tolower(depvar)),ncomp), moments );
  //std::cout << "<r^2>: " << rhovar << std::endl;

  for (ncomp_t c=0; c<ncomp; ++c) {
    //S[c] = 1.0/(1.0-YK[c]) - (1.0-Yt[c])/(1.0-YtK[c]);
    //S[c] = YK[c]/(1.0-YK[c]) - (1.0-Yt[c])*YtK[c]/(1.0-YtK[c]) + Yt[c];
    //S[c] = Yt[c] / ( 1.0 - sumYt + Yt[c] );
    S[c] = ( -2.0*(r[c]/rho[ncomp]*RRY[c])*(1.0-sumYt) +
             (r[c]/rho[ncomp]*(rhovar-sumRRY))*Yt[c] ) /
           ( -2.0*(r[c]/rho[ncomp]*RRY[c])*(1.0-sumYt) -
             (1.0-sumYt-Yt[c])*(r[c]/rho[ncomp]*(rhovar-sumRRY)) );
    //std::cout << "S: " << S[c] << ", YKc: " << YK[c]
    //          << ", Ytc: " << Yt[c] << ", YtKc: " << YtK[c] << ' ';
    k[c] = kprime[c] * bc[c];
    //if (k[c] < 0.0)
    // std::cout << "Positivity of k[" << c << "] violated: "
    //           << k[c] << '\n';
  }
  //std::cout << std::endl;

  for (ncomp_t c=0; c<ncomp; ++c) {
    if (S[c] < 0.0 || S[c] > 1.0) {
      //std::cout << "S[" << c << "] bounds violated: " << S[c] << '\n';
      //S[c] = 0.5;
    }
  }
  //std::cout << std::endl;
}
