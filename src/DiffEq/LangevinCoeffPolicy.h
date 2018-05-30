// *****************************************************************************
/*!
  \file      src/DiffEq/LangevinCoeffPolicy.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Langevin equation coefficients policies
  \details   This file defines coefficients policy classes for the Langevin
     equation for the fluctuating velocity in variable-density turbulence,
     defined in DiffEq/Langevin.h.

    General requirements on Langevin equation coefficients policy classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficient, C0. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          kw::sde_c0::info::expect::type C0_,
          kw::sde_c0::info::expect::type& C0 )
      \endcode
      where
      - C0_ denote a real value used to initialize the Langevin system.
      - The reference C0 is to be initialized based on C0_.

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
#ifndef LangevinCoeffPolicy_h
#define LangevinCoeffPolicy_h

#include <array>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "SystemComponents.h"
#include "Walker/Options/CoeffPolicy.h"

namespace walker {

//! Langevin equation coefficients policy with DNS hydrodynamics time scale
//! \details C0 is user-defined and we pull in a hydrodynamic timescale from an
//!   external function (from DNS).
//! \see kw::hydrotimescale_info
class LangevinCoeffHydroTimeScale {

  public:
    //! Constructor: initialize coefficients
    LangevinCoeffHydroTimeScale( kw::sde_c0::info::expect::type C0_,
                                 kw::sde_c0::info::expect::type& C0 )
    {
      C0 = C0_;
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::HYDROTIMESCALE; }

    //! \brief Update the dissipation rate (eps)
    //! \details Update the dissipation rate (eps) based on eps/k (from DNS) and
    //!    the turbulent kinetic energy (k) (from the SDE)
    void update( char depvar,
                 const std::map< tk::ctr::Product, tk::real >& moments,
                 const tk::Table& hts,
                 kw::sde_c0::info::expect::type C0,
                 tk::real t,
                 tk::real& eps,
                 std::array< tk::real, 9 >& G ) const
    {
      using tk::ctr::lookup;
      using tk::ctr::mean;

      const tk::ctr::Term u1( static_cast<char>(std::toupper(depvar)),
                              0, tk::ctr::Moment::CENTRAL );
      const tk::ctr::Term u2( static_cast<char>(std::toupper(depvar)),
                              1, tk::ctr::Moment::CENTRAL );
      const tk::ctr::Term u3( static_cast<char>(std::toupper(depvar)),
                              2, tk::ctr::Moment::CENTRAL );

      // Extract diagonal of the Reynolds stress
      const auto R11 = tk::ctr::Product( {u1,u1} );
      const auto R22 = tk::ctr::Product( {u2,u2} );
      const auto R33 = tk::ctr::Product( {u3,u3} );
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
      G[0] = G[4] = G[8] = -(1.0/2.0 + 3.0/4.0*C0) * ts;
    }

    //! Sample the inverse hydrodynamics time scale at time t
    //! \param[in] t Time at which to sample inverse hydrodynamics time scale
    //! \param[in] ts Hydro time scale table to sample
    //! \return Sampled value from discrete table of inverse hydro time scale
    tk::real hydrotimescale( tk::real t, const tk::Table& ts ) const
    { return tk::sample( t, ts ); }
};

//! List of all beta's coefficients policies
using LangevinCoeffPolicies =
  boost::mpl::vector< LangevinCoeffHydroTimeScale >;

} // walker::

#endif // LangevinCoeffPolicy_h
