// *****************************************************************************
/*!
  \file      src/DiffEq/Langevin.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functionality implementing Langevin models for the velocity
  \details   Functionality implementing Langevin models for the velocity.
*/
// *****************************************************************************

#include <cctype>
#include <array>

#include "StatCtr.h"
#include "Langevin.h"

std::array< tk::real, 9 >
walker::slm( tk::real hts, tk::real C0 )
// *****************************************************************************
//  Calculate the 2nd order tensor Gij based on the simplified Langevin model
//! \param[in] hts Inverse hydrodynamics time scale, e.g., eps/k
//! \param[in] C0 Coefficient C0 in SLM
//! \return Tensor Gij computed based on the simplified Langevin model
// *****************************************************************************
{
  std::array< tk::real, 9 > G;
  G.fill( 0.0 );
  G[0] = G[4] = G[8] = -(0.5+0.75*C0) * hts;

  return G;
}

std::array< tk::real, 9 >
walker::glm( tk::real hts,
             tk::real C0,
             const std::array< tk::real, 6 >& rs,
             const std::array< tk::real, 9 >& dU )
// *****************************************************************************
//  Calculate the 2nd order tensor Gij based on the generalized Langevin model
//! \param[in] hts Inverse hydrodynamics time scale, e.g., eps/k
//! \param[in] C0 Coefficient C0 in SLM
//! \param[in] rs Reynolds stress
//! \param[in] dU Mean velocity gradient
//! \return Tensor Gij computed based on the simplified Langevin model
// *****************************************************************************
{
  // Generalized Langevion model coefficients
  std::array< tk::real, 2 > ALPHA{{ -(0.5 + 0.75*C0), 3.7 }};
  std::array< tk::real, 3 > BETA{{ -0.2, 0.8, -0.2 }};
  std::array< tk::real, 6 > GAMMA{{ -1.28, 3.01, -2.18, 0.0, 4.29, -3.09 }};

  // Compute Reynolds stress anisotropy
  tk::real tr = rs[0] + rs[1] + rs[2];
  std::array< tk::real, 9 > b{{ rs[0]/tr-1.0/3.0,
                                rs[3]/tr,
                                rs[4]/tr,
                                rs[3]/tr,
                                rs[1]/tr-1.0/3.0,
                                rs[5]/tr,
                                rs[4]/tr,
                                rs[5]/tr,
                                rs[2]/tr-1.0/3.0 }};
  // Compute Gij
  std::array< tk::real, 9 > G;
  G.fill( 0.0 );
  for (std::size_t i=0; i<3; ++i ) {
    // to main diagonal: hts * ALPHA1 * delta_ij + BETA1 * deltaij * d<Ul>/dxl
    G[i*3+i] += hts*ALPHA[0] + BETA[0]*(dU[0] + dU[4] + dU[8]);
    // to main diagonal: GAMMA1 * deltaij * bml * d<Um>/dxl
    tk::real dtmp = 0.0;
    for (std::size_t m=0; m<3; ++m)
      for (std::size_t l=0; l<3; ++l )
         dtmp += b[m*3+l]*dU[m*3+l];
    G[i*3+i] += GAMMA[0]*dtmp;
    // to main and off-diagonal
    for (std::size_t j=0; j<3; ++j) {
      G[i*3+j] += hts*ALPHA[1]*b[i*3+j] +  // eps/k * ALPHA2 * bij
                  BETA[1]*dU[i*3+j] +      // BETA2 * d<Ui>/dj
                  BETA[2]*dU[j*3+i] +      // BETA3 * d<Uj>/di
                  // GAMMA4 * bij * d<Ul>/dxl
                  GAMMA[3]*b[i*3+j]*(dU[0] + dU[4] + dU[8]);
      for (std::size_t l=0; l<3; ++l)
        G[i*3+j] += GAMMA[1]*b[j*3+l]*dU[i*3+l] + // GAMMA2 * bjl * d<Ui>/dxl
                    GAMMA[2]*b[j*3+l]*dU[l*3+i] + // GAMMA3 * bjl * d<Ul>/dxi
                    GAMMA[4]*b[i*3+l]*dU[l*3+j] + // GAMMA5 * bil * d<Ul>/dxj
                    GAMMA[5]*b[i*3+l]*dU[j*3+l];  // GAMMA6 * bil * d<Uj>/dxl
    }
  }

  return G;
}

std::array< tk::real, 6 >
walker::reynoldsStress( char depvar,
                        ctr::DepvarType solve,
                        const std::map< tk::ctr::Product, tk::real >& moments )
// *****************************************************************************
// Compute the Reynolds stress tensor
//! \param[in] depvar Dependent variable labeling a velocity eq
//! \param[in] solve Enum selecting what the velocity eq solved for
//!   (fluctuating velocity of full/instantaneous veelocity)
//! \param[in] moments Map of statistical moments
//! \return Symmetric part of the Reynolds stress tensor
// *****************************************************************************
{
  using tk::ctr::Product;

  // Extract diagonal of the Reynolds stress
  Product r11, r22, r33, r12, r13, r23;
  if (solve == ctr::DepvarType::FULLVAR) {

    using tk::ctr::variance;
    using tk::ctr::covariance;
    r11 = variance( depvar, 0 );
    r22 = variance( depvar, 1 );
    r33 = variance( depvar, 2 );
    r12 = covariance( depvar, 0, depvar, 1 );
    r13 = covariance( depvar, 0, depvar, 2 );
    r23 = covariance( depvar, 1, depvar, 2 );

  } else if (solve == ctr::DepvarType::FLUCTUATION) {

    // Since we are solving for the fluctuating velocity, the "ordinary"
    // moments, e.g., <U1U1>, are really central moments, i.e., <u1u1>.
    using tk::ctr::Term;
    using tk::ctr::Moment;
    auto d = static_cast< char >( std::toupper( depvar ) );
    Term u( d, 0, Moment::ORDINARY );
    Term v( d, 1, Moment::ORDINARY );
    Term w( d, 2, Moment::ORDINARY );
    r11 = Product( { u, u } );
    r22 = Product( { v, v } );
    r33 = Product( { w, w } );
    r12 = Product( { u, v } );
    r13 = Product( { u, w } );
    r23 = Product( { v, w } );

  } else Throw( "Depvar type not implemented" );

  // Compute nonzero components of the Reynolds stress
  return {{ lookup(r11,moments), lookup(r22,moments), lookup(r33,moments),
            lookup(r12,moments), lookup(r13,moments), lookup(r23,moments) }};
}

tk::real
walker::tke( char depvar,
             ctr::DepvarType solve,
             const std::map< tk::ctr::Product, tk::real >& moments )
// *****************************************************************************
// Compute the turbulent kinetic energy
//! \param[in] depvar Dependent variable labeling a velocity eq
//! \param[in] solve Enum selecting what the velocity eq solved for
//!   (fluctuating velocity of full/instantaneous veelocity)
//! \param[in] moments Map of statistical moments
//! \return Turbulent kinetic energy
// *****************************************************************************
{
  using tk::ctr::Product;

  // Extract diagonal of the Reynolds stress
  Product r11, r22, r33;
  if (solve == ctr::DepvarType::FULLVAR) {

    using tk::ctr::variance;
    using tk::ctr::covariance;
    r11 = variance( depvar, 0 );
    r22 = variance( depvar, 1 );
    r33 = variance( depvar, 2 );

  } else if (solve == ctr::DepvarType::FLUCTUATION) {

    // Since we are solving for the fluctuating velocity, the "ordinary"
    // moments, e.g., <U1U1>, are really central moments, i.e., <u1u1>.
    using tk::ctr::Term;
    using tk::ctr::Moment;
    auto d = static_cast< char >( std::toupper( depvar ) );
    Term u( d, 0, Moment::ORDINARY );
    Term v( d, 1, Moment::ORDINARY );
    Term w( d, 2, Moment::ORDINARY );
    r11 = Product( { u, u } );
    r22 = Product( { v, v } );
    r33 = Product( { w, w } );

  } else Throw( "Depvar type not implemented" );

  // Compute nonzero components of the Reynolds stress
  return (lookup(r11,moments) + lookup(r22,moments) + lookup(r33,moments))/2.0;
}
