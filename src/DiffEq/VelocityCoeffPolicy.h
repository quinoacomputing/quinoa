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
                     kw::sde_c0::info::expect::type C0,
                     tk::real t,
                     tk::real& eps,
                     std::array< tk::real, 9 >& G ) const
      \endcode
      where _depvar_ is the dependent variable of the velocity equation,
      _dissipation_depvar_ is the dependent variable of the coupled dissipation
      equation, _moments_ if the map of computed statistical moments, _hts_ is
      a ctr::Table containing the inverse hydrodynamic timescale, _solve_ is the
      the lable of the dependent variable to solve for (full variable or
      fluctuation), _C0_ is the Langevin eq constat to use, _t_ is the physical
      time, _eps_ is the dissipation rate of turbulent kinetic energy to update,
      and _G_ is the G_{ij} tensor in the Langevin model to update.

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

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "SystemComponents.h"
#include "Walker/Options/CoeffPolicy.h"

namespace walker {

//! Velocity equation coefficients policy with prescribed mean shear
//! \details C0 is user-defined and we prescibe a hard-coded mean shear in the x
//!   direction
//! \see kw::hydrotimescale_info
class Velocity_ConstShear {

  public:
    //! Constructor: initialize coefficients
    //! \param[in] C0_ Value of C0 parameter in the Langevin model
    //! \param[in,out] C0 Value of to set the C0 parameter in the Langevin model
    //! \param[in,out] dU Prescribed mean velocity gradient1
    Velocity_ConstShear( kw::sde_c0::info::expect::type C0_,
                         kw::sde_c0::info::expect::type& C0,
                         std::array< tk::real, 9 >& dU )
    {
      C0 = C0_;
      dU = {{ 0.0, 1.0, 0.0,
              0.0, 0.0, 0.0,
              0.0, 0.0, 0.0 }};
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
                 kw::sde_c0::info::expect::type C0,
                 tk::real,
                 tk::real& eps,
                 std::array< tk::real, 9 >& G ) const
    {
      using tk::ctr::lookup;
      using tk::ctr::mean;
      using tk::ctr::Product;

      // Extract diagonal of the Reynolds stress
      Product r11, r22, r33;
      if (solve == ctr::DepvarType::FULLVAR) {

        using tk::ctr::variance;
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
        r11 = tk::ctr::Product( { u, u } );
        r22 = tk::ctr::Product( { v, v } );
        r33 = tk::ctr::Product( { w, w } );

      } else Throw( "Depvar type not implemented" );

      // compute turbulent kinetic energy
      tk::real k = ( lookup(r11,moments) +
                     lookup(r22,moments) +
                     lookup(r33,moments) ) / 2.0;

      // Access mean turbulence frequency
      tk::real O = lookup( mean(dissipation_depvar,0), moments );

      // compute turbulent kinetic energy dissipation rate
      eps = O*k;

      // update drift tensor based on the simplified Langevin model
      G.fill( 0.0 );
      G[0] = G[4] = G[8] = -(0.5+0.75*C0) * O;
    }
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
using VelocityCoeffPolicies = boost::mpl::vector< Velocity_HydroTimeScale
                                                , Velocity_ConstShear
                                                >;

} // walker::

#endif // VelocityCoeffPolicy_h
