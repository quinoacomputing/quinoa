// *****************************************************************************
/*!
  \file      src/DiffEq/VelocityCoeffPolicy.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Velocity equation coefficients policies
  \details   This file defines coefficients policy classes for the velocity
     equation for the fluctuating velocity in variable-density turbulence,
     defined in DiffEq/Velocity.h.

    General requirements on velocity equation coefficients policy classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficient, C0. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          kw::sde_c0::info::expect::type C0_,
          kw::sde_c0::info::expect::type& C0,
          std::array< tk::real, 9 >& dU )
      \endcode
      where
      - C0_ denotes a real value used to initialize the velocity system.
      - The reference C0 is to be initialized based on C0_.
      - _dU_ is an optionally prescribed mean velocity gradient.

    - Must define the function _update()_, called from Velocity::advance(),
      updating the model coefficients.
      Required signature:
      \code{.cpp}
        void update( char depvar,
                     char dissipation_depvar,
                     const std::map< tk::ctr::Product, tk::real >& moments,
                     const tk::Table& hts,
                     ctr::DepvarType solve,
                     ctr::VelocityVariantType variant,
                     kw::sde_c0::info::expect::type C0,
                     tk::real t,
                     tk::real& eps,
                     std::array< tk::real, 9 >& G ) const
      \endcode
      where _depvar_ is the dependent variable of the velocity equation,
      _dissipation_depvar_ is the dependent variable of the coupled dissipation
      equation, _moments_ if the map of computed statistical moments, _hts_ is
      a ctr::Table containing the inverse hydrodynamic timescale, _variant_ is
      the velocity model variant, _solve_ is the the lable of the dependent
      variable to solve for (full variable or fluctuation), _C0_ is the Langevin
      eq constat to use, _t_ is the physical time, _eps_ is the dissipation rate
      of turbulent kinetic energy to update, and _G_ is the G_{ij} tensor in the
      Langevin model to update.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::CoeffPolicyType type() noexcept {
          return ctr::CoeffPolicyType::HYDROTIMESCALE;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef VelocityCoeffPolicy_h
#define VelocityCoeffPolicy_h

#include <array>

#include <brigand/sequences/list.hpp>

#include "Types.h"
#include "SystemComponents.h"
#include "Walker/Options/CoeffPolicy.h"

namespace walker {

//! Calculate the 2nd order tensor Gij based on the simplified Langevin model
//! \param[in] hts Inverse hydrodynamics time scale, e.g., eps/k
//! \param[in] C0 Coefficient C0 in SLM
//! \return Tensor Gij computed based on the simplified Langevin model
static inline std::array< tk::real, 9 >
slm( tk::real hts, tk::real C0 ) {

  std::array< tk::real, 9 > G;
  G.fill( 0.0 );
  G[0] = G[4] = G[8] = -(0.5+0.75*C0) * hts;

  return G;
}

//! Calculate the 2nd order tensor Gij based on the generalized Langevin model
//! \param[in] hts Inverse hydrodynamics time scale, e.g., eps/k
//! \param[in] C0 Coefficient C0 in SLM
//! \param[in] rs Reynolds stress
//! \param[in] dU Mean velocity gradient
//! \return Tensor Gij computed based on the simplified Langevin model
static inline std::array< tk::real, 9 >
glm( tk::real hts,
     tk::real C0,
     const std::array< tk::real, 6 >& rs,
     const std::array< tk::real, 9 >& dU )
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

//! Velocity equation coefficients policy with prescribed mean shear
//! \details C0 is user-defined and we prescibe a hard-coded mean shear in the x
//!   direction
//! \see kw::const_shear_info
class Velocity_ConstShear {

  public:
    //! Constructor: initialize coefficients
    //! \param[in] C0_ Value of C0 parameter in the Langevin model
    //! \param[in,out] C0 Value of to set the C0 parameter in the Langevin model
    //! \param[in,out] dU Prescribed mean velocity gradient1
    Velocity_ConstShear( kw::sde_c0::info::expect::type C0_,
                         kw::sde_c0::info::expect::type& C0,
                         std::array< tk::real, 9 >& dU ) :
      m_dU( {{ 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0 }} )
    {
      C0 = C0_;
      dU = m_dU;
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONST_SHEAR; }

    //! Update the model coefficients (prescribing shear)
    //! \details Update the dissipation rate (eps) and G_{ij} based on the
    //!   turbulent kinetic energy (k) for a prescribed honmogeneous shear flow.
    void update( char depvar,
                 char dissipation_depvar,
                 const std::map< tk::ctr::Product, tk::real >& moments,
                 const tk::Table&,
                 ctr::DepvarType solve,
                 ctr::VelocityVariantType variant,
                 kw::sde_c0::info::expect::type C0,
                 tk::real,
                 tk::real& eps,
                 std::array< tk::real, 9 >& G ) const
    {
      using tk::ctr::lookup;
      using tk::ctr::mean;
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
        r11 = tk::ctr::Product( { u, u } );
        r22 = tk::ctr::Product( { v, v } );
        r33 = tk::ctr::Product( { w, w } );
        r12 = tk::ctr::Product( { u, v } );
        r13 = tk::ctr::Product( { u, w } );
        r23 = tk::ctr::Product( { v, w } );

      } else Throw( "Depvar type not implemented" );

      // Compute nonzero components of the Reynolds stress
      std::array< tk::real, 6 > rs{{ lookup(r11,moments),
                                     lookup(r22,moments),
                                     lookup(r33,moments),
                                     lookup(r12,moments),
                                     lookup(r13,moments),
                                     lookup(r23,moments) }};

      // Compute turbulent kinetic energy
      tk::real k = (rs[0] + rs[1] + rs[2]) / 2.0;

      // Access mean turbulence frequency
      tk::real O = lookup( mean(dissipation_depvar,0), moments );

      // compute turbulent kinetic energy dissipation rate
      eps = O*k;

      // update drift tensor based on the Langevin model variant configured
      if (variant == ctr::VelocityVariantType::SLM)     // simplified
        G = slm( O, C0 );
      else if (variant == ctr::VelocityVariantType::GLM)// generalized
        G = glm( O, C0, rs, m_dU );
      else Throw( "Velocity variant type not implemented" );
    }

  private:
    //! Mean velocity gradient prescribed for simpled 1D homogeneous shear
    std::array< tk::real, 9 > m_dU;
};

//! Velocity equation coefficients policy with prescribed constant dissipation
//! \details C0 is user-defined and we prescibe a hard-coded dissipation rate
//! \see kw::const_dissipation_info
class Velocity_ConstDissipation {

  public:
    //! Constructor: initialize coefficients
    //! \param[in] C0_ Value of C0 parameter in the Langevin model
    //! \param[in,out] C0 Value of to set the C0 parameter in the Langevin model
    //! \param[in,out] dU Prescribed mean velocity gradient1
    Velocity_ConstDissipation( kw::sde_c0::info::expect::type C0_,
                               kw::sde_c0::info::expect::type& C0,
                               std::array< tk::real, 9 >& dU ) :
      m_dU( {{ 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0 }} )
    {
      C0 = C0_;
      dU = m_dU;
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONST_DISSIPATION; }

    //! Update the model coefficients (prescribing shear)
    //! \details Update the dissipation rate (eps) and G_{ij} based on the
    //!   turbulent kinetic energy (k) for a prescribed honmogeneous shear flow.
    void update( char depvar,
                 char,
                 const std::map< tk::ctr::Product, tk::real >& moments,
                 const tk::Table&,
                 ctr::DepvarType solve,
                 ctr::VelocityVariantType variant,
                 kw::sde_c0::info::expect::type C0,
                 tk::real,
                 tk::real& eps,
                 std::array< tk::real, 9 >& G ) const
    {
      using tk::ctr::lookup;
      using tk::ctr::mean;
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
        r11 = tk::ctr::Product( { u, u } );
        r22 = tk::ctr::Product( { v, v } );
        r33 = tk::ctr::Product( { w, w } );
        r12 = tk::ctr::Product( { u, v } );
        r13 = tk::ctr::Product( { u, w } );
        r23 = tk::ctr::Product( { v, w } );

      } else Throw( "Depvar type not implemented" );

      // Compute nonzero components of the Reynolds stress
      std::array< tk::real, 6 > rs{{ lookup(r11,moments),
                                     lookup(r22,moments),
                                     lookup(r33,moments),
                                     lookup(r12,moments),
                                     lookup(r13,moments),
                                     lookup(r23,moments) }};

      // Assign constant mean turbulence frequency
      tk::real O = 1.0;

      // Assign constant turbulent kinetic energy dissipation rate
      eps = 1.0;

      // update drift tensor based on the Langevin model variant configured
      if (variant == ctr::VelocityVariantType::SLM)     // simplified
        G = slm( O, C0 );
      else if (variant == ctr::VelocityVariantType::GLM)// generalized
        G = glm( O, C0, rs, m_dU );
      else Throw( "Velocity variant type not implemented" );
    }

  private:
    //! Mean velocity gradient prescribed for simpled 1D homogeneous shear
    std::array< tk::real, 9 > m_dU;
};

//! Velocity equation coefficients policy with DNS hydrodynamics time scale
//! \details C0 is user-defined and we pull in a hydrodynamic timescale from an
//!   external function (from DNS).
//! \see kw::hydrotimescale_info
class Velocity_HydroTimeScale {

  public:
    //! Constructor: initialize coefficients
    //! \param[in] C0_ Value of C0 parameter in the Langevin model
    //! \param[in,out] C0 Value of to set the C0 parameter in the Langevin model
    Velocity_HydroTimeScale( kw::sde_c0::info::expect::type C0_,
                             kw::sde_c0::info::expect::type& C0,
                             std::array< tk::real, 9 >& )
    {
      C0 = C0_;
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::HYDROTIMESCALE; }

    //! \brief Update the model coefficients
    //! \details Update the dissipation rate (eps) based on eps/k (from DNS) and
    //!    the turbulent kinetic energy (k) (from the SDE)
    void update( char depvar,
                 char,
                 const std::map< tk::ctr::Product, tk::real >& moments,
                 const tk::Table& hts,
                 ctr::DepvarType,
                 ctr::VelocityVariantType,
                 kw::sde_c0::info::expect::type C0,
                 tk::real t,
                 tk::real& eps,
                 std::array< tk::real, 9 >& G ) const
    {
      using tk::ctr::lookup;
      using tk::ctr::variance;

      // Extract diagonal of the Reynolds stress
      const auto R11 = variance( depvar, 0 );
      const auto R22 = variance( depvar, 1 );
      const auto R33 = variance( depvar, 2 );
      // compute turbulent kinetic energy
      tk::real k = ( lookup(R11,moments) +
                     lookup(R22,moments) +
                     lookup(R33,moments) ) / 2.0;

      // Sample hydrodynamics timescale and prod/diss at time t
      auto ts = hydrotimescale( t, hts );  // eps/k

      // compute turbulent kinetic energy dissipation rate
      eps = ts * k;

      // update drift tensor based on the simplified Langevin model
      G.fill( 0.0 );
      G[0] = G[4] = G[8] = -(0.5+0.75*C0) * ts;
    }

    //! Sample the inverse hydrodynamics time scale at time t
    //! \param[in] t Time at which to sample inverse hydrodynamics time scale
    //! \param[in] ts Hydro time scale table to sample
    //! \return Sampled value from discrete table of inverse hydro time scale
    tk::real hydrotimescale( tk::real t, const tk::Table& ts ) const
    { return tk::sample( t, ts ); }
};

//! List of all Velocity's coefficients policies
using VelocityCoeffPolicies = brigand::list< Velocity_HydroTimeScale
                                           , Velocity_ConstShear
                                           , Velocity_ConstDissipation
                                           >;

} // walker::

#endif // VelocityCoeffPolicy_h
