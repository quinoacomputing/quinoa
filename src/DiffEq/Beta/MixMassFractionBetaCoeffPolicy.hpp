// *****************************************************************************
/*!
  \file      src/DiffEq/Beta/MixMassFractionBetaCoeffPolicy.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mix mass-fraction beta SDE coefficients policies
  \details   This file declares coefficients policy classes for the mix
    mass-fraction beta SDE, defined in DiffEq/Beta/MixMassFractionBeta.h.

    General requirements on mix mass-fraction beta SDE coefficients policy
    classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, b, S, kappa, rho2, and r. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_t ncomp,
          const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
          const std::vector< kw::sde_S::info::expect::type >& S_,
          const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
          const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
          const std::vector< kw::sde_r::info::expect::type >& r_,
          std::vector< kw::sde_bprime::info::expect::type  >& bprime,
          std::vector< kw::sde_S::info::expect::type >& S,
          std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
          std::vector< kw::sde_rho2::info::expect::type >& rho2_,
          std::vector< kw::sde_r::info::expect::type >& r_ );
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of
        mix mass-fraction beta SDEs.
      - Constant references to bprime_, S_, kprime_, rho2_, and r_, which
        denote vectors of real values used to initialize the parameter
        vectors of the system of mix mass-fraction beta SDEs. The length of
        the vectors must be equal to the number of components given by ncomp.
      - References to bprime, S, kprime, rho2_, and r, which denote the
        parameter vectors to be initialized based on bprime_, S_, kprime_,
        rho2_, and r_.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::CoeffPolicyType type() noexcept {
          return ctr::CoeffPolicyType::DECAY;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.

    - Must define the function _update()_, called from
      MixMassFractionBeta::advance(), updating the model coefficients.
      Required signature:
      \code{.cpp}
        void update(
          char depvar,
          char dissipation_depvar,
          char velocity_depvar,
          ctr::DepvarType velocity_solve,
          ctr::DepvarType solve,
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
          std::vector< kw::sde_S::info::expect::type >& S ) const {}
      \endcode
      where _depvar_ is the dependent variable associated with the mix
      mass-fraction beta SDE, specified in the control file by the user,
      _dissipation_depvar_ is a character labeling the coupled dissipation
      equation, _velocity_depvar_ is an is a character labeling the coupled
      velocity equation, _velocity_solve_ is an enum selecting whether the
      coupled velocity equation solves for full variable or its fluctuation,
      _solve_ is an enum selecting whether the mixmassfracbeta (scalar) equation
      solves for full variable or its fluctuation, _ncomp_ is the number of
      components in the system, _moments_ is the map associating moment IDs
      (tk::ctr::vector< tk::ctr::Term >) to values of statistical moments,
      _bprime_, _kprime_, rho2, r, are user-defined parameters, and _b_, _k_,
      _S_, are the SDE parameters computed, see
      DiffEq/Beta/MixMassFractionBeta.h.

      The constant reference to hts, denotes a vector of y=f(x) functions (see
      src/DiffEq/Beta/HydroTimeScales.h and
      src/Control/Walker/Options/HydroTimescales.h) used to configure the
      inverse hydrodynamics time scales (extracted from direct numerical
      simulations) of the system of mix mass-fraction beta SDEs if the
      MixMassFracBetaCoeffHydroTimeScaleHomDecay coefficients policy is
      selected. The length of this vector must be equal to the number of
      components given by ncomp. Note that hts is only used by
      MixMassFracBetaCoeffHydroTimeScaleHomDecay.

      The constant reference to hp, denotes a vector of y=f(x) functions (see
      src/DiffEq/Beta/HydroProductions.h and
      src/Control/Walker/Options/HydroProductions.h) used to configure the
      turbulent kinetic energy production divided by the dissipation rate, P/e,
      a measure of the non-eqilibrium nature of the turbulent flow (extracted
      from direct numerical simulations) of the system of mix mass-fraction beta
      SDEs if the MixMassFracBetaCoeffHydroTimeScaleHomDecay
      coefficients policy is selected. The length of this vector must be equal
      to the number of components given by ncomp. Note that hp is only used by
      MixMassFracBetaCoeffHydroTimeScaleHomDecay.
*/
// *****************************************************************************
#ifndef MixMassFractionBetaCoeffPolicy_h
#define MixMassFractionBetaCoeffPolicy_h

#include <brigand/sequences/list.hpp>

#include "Types.hpp"
#include "Table.hpp"
#include "Walker/Options/CoeffPolicy.hpp"
#include "Langevin.hpp"

namespace walker {

//! \brief Mix mass-fraction beta SDE decay coefficients policy
//! \details User-defined parameters b' and kappa' are constants in time and
//!   ensure decay in the evolution of <y^2>.
class MixMassFracBetaCoeffDecay {

    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor: initialize coefficients
    MixMassFracBetaCoeffDecay(
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
      std::vector< kw::sde_kappa::info::expect::type >& k );

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::DECAY; }

    //! Update coefficients
    void update(
      char depvar,
      char,
      char,
      ctr::DepvarType /*velocity_solve*/,
      ctr::DepvarType /*solve*/,
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
      tk::real ) const;
};

//! \brief Mix mass-fraction beta SDE homogneous decay coefficients policy
//! \details User-defined parameters b' and kappa' are constants in time and
//!   ensure decay in the evolution of \<y^2\>. Additionally, S is constrained
//!   to make d\<rho\>/dt = 0, where \<rho\> = rho_2/(1+rY).
class MixMassFracBetaCoeffHomDecay {

    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor: initialize coefficients
    MixMassFracBetaCoeffHomDecay(
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
      std::vector< kw::sde_kappa::info::expect::type >& k );

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::HOMOGENEOUS_DECAY; }

    //! Update coefficients
    void update(
      char depvar,
      char,
      char,
      ctr::DepvarType /*velocity_solve*/,
      ctr::DepvarType /*solve*/,
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
      tk::real ) const;
};

//! \brief Mix mass-fraction beta SDE Monte Carlo homogenous decay coefficients
//!   policy
//! \details User-defined parameters b' and kappa' are constants in time and
//!   ensure decay in the evolution of \<y^2\>. Additionally, S is constrained
//!   to make d\<rho\>/dt = 0, where \<rho\> = rho_2/(1+rY). This is the same as
//!   the specification in MixMassFracBetaCoeffHomDecay, but uses more advanced
//!   statistics, available from the Monte Carlo simulation, which yield a
//!   simpler formula for the coefficient S.
class MixMassFracBetaCoeffMonteCarloHomDecay {

    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor: initialize coefficients
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
      std::vector< kw::sde_kappa::info::expect::type >& k );

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::MONTE_CARLO_HOMOGENEOUS_DECAY; }

    //! Update coefficients
    void update(
      char depvar,
      char,
      char,
      ctr::DepvarType /*velocity_solve*/,
      ctr::DepvarType /*solve*/,
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
      tk::real ) const;
};

//! \brief Mix mass-fraction beta SDE coefficients policy with DNS hydrodynamics
//!   time scale
//! \details User-defined parameters b' and kappa' are functions of an
//!   externally, e.g., DNS-, provided hydrodynamics time scale ensuring decay
//!   in the evolution of \<y^2\>. Additionally, S is constrained to
//!   make d\<rho\>/dt = 0, where \<rho\> = rho_2/(1+rY). Additionally,
//!   we pull in a hydrodynamic timescale from an external function.
//! \see kw::hydrotimescale_info
class MixMassFracBetaCoeffHydroTimeScale {

    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor: initialize coefficients
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
      std::vector< kw::sde_kappa::info::expect::type >& k );

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::HYDROTIMESCALE; }

    //! Update coefficients b', kappa', and S
    void update(
      char depvar,
      char,
      char,
      ctr::DepvarType /*velocity_solve*/,
      ctr::DepvarType /*solve*/,
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
      tk::real t ) const;

    //! Sample the inverse hydrodynamics time scale at time t
    //! \param[in] t Time at which to sample inverse hydrodynamics time scale
    //! \param[in] ts Hydro time scale table to sample
    //! \return Sampled value from discrete table of inverse hydro time scale
    tk::real hydrotimescale( tk::real t, const tk::Table& ts ) const
    { return tk::sample( t, ts ); }

    //! Sample the hydrodynamics production/dissipation rate (P/e) at time t
    //! \param[in] t Time at which to sample hydrodynamics P/e
    //! \param[in] p P/e table to sample
    //! \return Sampled value from discrete table of P/e
    tk::real hydroproduction( tk::real t, const tk::Table& p ) const
    { return tk::sample( t, p ); }

    mutable std::size_t m_it = 0;
    mutable std::vector< tk::real > m_s;
    mutable std::string m_extra_out_filename;
};

//! \brief Mix mass-fraction beta SDE coefficients policy coupled to velocity
//! \details User-defined parameters b' and kappa' are functions of P/eps and
//!    k/eps from a coupled velocity model. Additionally, S is constrained to
//!   make d\<rho\>/dt = 0, where \<rho\> = rho_2/(1+rY).
class MixMassFracBetaCoeffInstVel {

    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor: initialize coefficients
    MixMassFracBetaCoeffInstVel(
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
      std::vector< kw::sde_kappa::info::expect::type >& k );

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::INSTANTANEOUS_VELOCITY; }

    //! Update coefficients
    void update(
      char depvar,
      char dissipation_depvar,
      char /*velocity_depvar*/,
      ctr::DepvarType /*velocity_solve*/,
      ctr::DepvarType solve,
      ncomp_t ncomp,
      const std::map< tk::ctr::Product, tk::real >& moments,
      const std::vector< kw::sde_bprime::info::expect::type  >& /*bprime*/,
      const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2,
      const std::vector< kw::sde_r::info::expect::type >& r,
      const std::vector< tk::Table >&,
      const std::vector< tk::Table >&,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_kappa::info::expect::type >& k,
      std::vector< kw::sde_S::info::expect::type >& S,
      tk::real ) const;

    mutable std::size_t m_it = 0;
    mutable std::vector< tk::real > m_s;
};

//! List of all mix mass-fraction beta's coefficients policies
using MixMassFracBetaCoeffPolicies =
  brigand::list< MixMassFracBetaCoeffDecay
               , MixMassFracBetaCoeffHomDecay
               , MixMassFracBetaCoeffMonteCarloHomDecay
               , MixMassFracBetaCoeffHydroTimeScale
               , MixMassFracBetaCoeffInstVel
               >;

} // walker::

#endif // MixMassFractionBetaCoeffPolicy_h
