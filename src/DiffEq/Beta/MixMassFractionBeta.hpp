// *****************************************************************************
/*!
  \file      src/DiffEq/Beta/MixMassFractionBeta.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     System of mix mass-fraction beta SDEs
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs) with linear drift and quadratic diagonal
    diffusion, whose invariant is the joint [beta
    distribution](http://en.wikipedia.org/wiki/Beta_distribution). There are two
    differences compared to the plain beta SDE (see DiffEq/Beta.h):

    - First, the parameters, b, and kappa are specified via functions that
    constrain the beta SDE to be consistent with the turbulent mixing process.
    In particular, the SDE is made consistent with the no-mix and fully mixed
    limits. See, e.g., MixMassFractionBetaCoeffConst::update().

    - Second, there two additional random variables computed, the same as also
    computed by the mass-fraction beta equation, see also
    DiffEq/MassFractionBeta.h.

    In a nutshell, the equation integrated governs a set of scalars,
    \f$0\!\le\!Y_\alpha\f$, \f$\alpha\!=\!1,\dots,N\f$, as

    @m_class{m-show-m}

    \f[
       \mathrm{d}Y_\alpha(t) = \frac{b_\alpha}{2}\left(S_\alpha - Y_\alpha\right)
       \mathrm{d}t + \sqrt{\kappa_\alpha Y_\alpha(1-Y_\alpha)}
       \mathrm{d}W_\alpha(t), \qquad \alpha=1,\dots,N
    \f]

    @m_class{m-hide-m}

    \f[ \begin{split}
       \mathrm{d}Y_\alpha(t) = \frac{b_\alpha}{2}\left(S_\alpha - Y_\alpha\right)
       \mathrm{d}t + \sqrt{\kappa_\alpha Y_\alpha(1-Y_\alpha)}
       \mathrm{d}W_\alpha(t), \\ \alpha=1,\dots,N
    \end{split} \f]

    with parameter vectors \f$b_\alpha = \Theta b'_\alpha > 0\f$, \f$
    \newcommand{\irv}[1]{\langle{#1^2}\rangle} \kappa_\alpha = \kappa' \irv{x} >
    0\f$, and \f$0 < S_\alpha < 1\f$. This is similar to DiffEq/Beta.h, but the
    parameters, \f$b\f$ and \f$\kappa\f$ constrained. Here \f$
    \newcommand{\irv}[1]{\langle{#1^2}\rangle}
    \newcommand{\irmean}[1]{{\langle{#1}\rangle}} \Theta = 1 - \irv{x} /
    [ \irmean{Y} (1-\irmean{Y}) ]\f$. The fluctuation about the mean, \f$
    \newcommand{\irmean}[1]{{\langle{#1}\rangle}} \irmean{Y} \f$, is defined as
    usual: \f$ \newcommand{\irmean}[1]{{\langle{#1}\rangle}} x = Y - \irmean{Y}
    \f$, and \f$b'\f$ and \f$ \kappa'\f$ are user-specified constants. Also,
    \f$\mathrm{d}W_\alpha(t)\f$ is an isotropic vector-valued
    [Wiener process](http://en.wikipedia.org/wiki/Wiener_process) with
    independent increments. The invariant distribution is the joint beta
    distribution. This system of SDEs consists of N independent equations. For
    more on the beta SDE, see https://doi.org/10.1080/14685248.2010.510843.

    In addition to integrating the above SDE, there are two additional functions
    of \f$ Y_\alpha \f$ are computed as
    \f[ \begin{aligned}
      \rho(Y_\alpha) & = \frac{ \rho_{2\alpha} }{ 1 + r_\alpha Y_\alpha } \\
      V(Y_\alpha) & = \frac{1}{ \rho(Y_\alpha) }
    \end{aligned} \f]
    These equations compute the instantaneous mixture density, \f$ \rho \f$, and
    instantaneous specific volume, \f$ V_\alpha \f$, for equation \f$ \alpha \f$
    in the system. These quantities are used in binary mixing of
    variable-density turbulence between two fluids with constant densities, \f$
    \rho_1, \f$ and \f$ \rho_2 \f$. The additional parameters, \f$ \rho_2 \f$
    and \f$ r \f$ are user input parameters and kept constant during
    integration. Since we compute the above variables, \f$\rho,\f$ and \f$V\f$,
    and call them mixture density and specific volume, respectively, \f$Y\f$,
    governed by the beta SDE is a mass fraction.

    _All of this is unpublished, but will be linked in here once published_.
*/
// *****************************************************************************
#ifndef MixMassFractionBeta_h
#define MixMassFractionBeta_h

#include <vector>
#include <tuple>

#include "InitPolicy.hpp"
#include "MixMassFractionBetaCoeffPolicy.hpp"
#include "RNG.hpp"
#include "Particles.hpp"
#include "Table.hpp"
#include "CoupledEq.hpp"
#include "HydroTimeScales.hpp"
#include "HydroProductions.hpp"
#include "Walker/Options/HydroTimeScales.hpp"
#include "Walker/Options/HydroProductions.hpp"

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief MixMassFractionBeta SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see
//!     DiffEq/MixMassFractionBetaCoeffPolicy.h
template< class Init, class Coefficients >
class MixMassFractionBeta {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::mixmassfracbeta;

  public:
    //! Number of derived variables computed by the SDE
    //! \details Derived variables: density, specific volume, 1 - mass fraction
    //! \warning If you change this, you must also change the constant in
    //!   Velocity::m_mixmassfracbeta_ncomp in DiffEq/Velocity/Velocity.hpp.
    //!   To see where, search for NUMDERIVED there.
    static const std::size_t NUMDERIVED = 3;

    //! \brief Constructor
    //! \param[in] c Index specifying which system of mix mass-fraction beta
    //!   SDEs to construct. There can be multiple mixmassfracbeta ... end blocks
    //!   in a control file. This index specifies which mix mass-fraction beta
    //!   SDE system to instantiate. The index corresponds to the order in which
    //!   the mixmassfracbeta ... end blocks are given the control file.
    explicit MixMassFractionBeta( ncomp_t c ) :
      m_c( c ),
      m_depvar( g_inputdeck.get< tag::param, eq, tag::depvar >().at(c) ),
      // divide by the number of derived variables computed, see derived()
      m_ncomp( g_inputdeck.get< tag::component >().get< eq >().at(c) /
               (NUMDERIVED + 1) ),
      m_offset( g_inputdeck.get< tag::component >().offset< eq >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, eq, tag::rng >().at(c) ) ) ),
      m_solve( g_inputdeck.get< tag::param, eq, tag::solve >().at(c) ),
      m_velocity_coupled( coupled< eq, tag::velocity >( c ) ),
      m_velocity_depvar( depvar< eq, tag::velocity >( c ) ),
      m_velocity_offset( offset< eq, tag::velocity, tag::velocity_id >( c ) ),
      m_velocity_solve(
        m_velocity_coupled ?
          g_inputdeck.get< tag::param, tag::velocity, tag::solve >().at(
            system_id< eq, tag::velocity, tag::velocity_id >( c ) ) :
        ctr::DepvarType::FULLVAR ),
      m_dissipation_coupled( coupled< eq, tag::dissipation >( c ) ),
      m_dissipation_depvar( depvar< eq, tag::dissipation >( c ) ),
      m_dissipation_offset(
        offset< eq, tag::dissipation, tag::dissipation_id >( c ) ),
      m_dY( initScalarGradient() ),
      m_bprime(),
      m_S(),
      m_kprime(),
      m_rho2(),
      m_r(),
      m_b(),
      m_k(),
      coeff(
        m_ncomp,
        g_inputdeck.get< tag::param, eq, tag::bprime >().at(c),
        g_inputdeck.get< tag::param, eq, tag::S >().at(c),
        g_inputdeck.get< tag::param, eq, tag::kappaprime >().at(c),
        g_inputdeck.get< tag::param, eq, tag::rho2 >().at(c),
        g_inputdeck.get< tag::param, eq, tag::r >().at(c),
        m_bprime, m_S, m_kprime, m_rho2, m_r, m_b, m_k )
    {
      // Zero prescribed scalar gradient if full variable is solved for
      if (m_solve == ctr::DepvarType::FULLVAR) m_dY.fill( 0.0 );
      // Populate inverse hydrodynamics time scales and hydrodyanmics
      // production/dissipation extracted from DNS
      if ( Coefficients::type() == ctr::CoeffPolicyType::HYDROTIMESCALE ) {
        // Configure inverse hydrodyanmics time scale from DNS
        const auto& hts =
          g_inputdeck.get< tag::param, eq, tag::hydrotimescales >().at(c);
        ctr::HydroTimeScales ot;
        // cppcheck-suppress useStlAlgorithm
        for (auto t : hts) m_hts.push_back( ot.table(t) );
        Assert( m_hts.size() == m_ncomp, "Number of inverse hydro time scale "
          "tables associated does not match the components integrated" );

        // Configure hydrodyanmics production/dissipation from DNS
        const auto& hp =
          g_inputdeck.get< tag::param, eq, tag::hydroproductions >().at(c);
        ctr::HydroProductions op;
        // cppcheck-suppress useStlAlgorithm
        for (auto t : hp) m_hp.push_back( op.table(t) );
        Assert( m_hp.size() == m_ncomp, "Number of hydro "
          "production/dissipation tables associated does not match the "
          "components integrated" );
      }
    }

    //! Initalize SDE, prepare for time integration
    //! \param[in] stream Thread (or more precisely stream) ID 
    //! \param[in,out] particles Array of particle properties 
    void initialize( int stream, tk::Particles& particles ) {
      //! Set initial conditions using initialization policy
      Init::template init< eq >
        ( g_inputdeck, m_rng, stream, particles, m_c, m_ncomp, m_offset );
      // Initialize values derived from primary prognostic variable
      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p)
        for (ncomp_t i=0; i<m_ncomp; ++i)
          derived( particles, p, i );
    }

    //! \brief Advance particles according to the system of mix mass-fraction
    //!   beta SDEs
    //! \param[in,out] particles Array of particle properties
    //! \param[in] stream Thread (or more precisely stream) ID
    //! \param[in] dt Time step size
    //! \param[in] t Physical time of the simulation
    //! \param[in] moments Map of statistical moments
    void advance( tk::Particles& particles,
                  int stream,
                  tk::real dt,
                  tk::real t,
                  const std::map< tk::ctr::Product, tk::real >& moments )
    {
      // Update SDE coefficients
      coeff.update( m_depvar, m_dissipation_depvar, m_velocity_depvar,
                    m_velocity_solve, m_solve, m_ncomp, moments, m_bprime,
                    m_kprime, m_rho2, m_r, m_hts, m_hp, m_b, m_k, m_S, t );

      const auto eps = std::numeric_limits< tk::real >::epsilon();

      // Advance particles
      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        std::vector< tk::real > dW( m_ncomp );
        m_rng.gaussian( stream, m_ncomp, dW.data() );

        // Access coupled particle velocity
        tk::real u = 0.0, v = 0.0, w = 0.0;
        using std::abs;
        if (abs(m_dY[0]) > eps || abs(m_dY[1]) > eps || abs(m_dY[2]) > eps) {
          u = particles( p, 0, m_velocity_offset );
          v = particles( p, 1, m_velocity_offset );
          w = particles( p, 2, m_velocity_offset );
        }

        // Advance all m_ncomp scalars
        for (ncomp_t i=0; i<m_ncomp; ++i) {
          tk::real& Y = particles( p, i, m_offset );
          tk::real d = m_k[i] * Y * (1.0 - Y) * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          Y += 0.5*m_b[i]*(m_S[i] - Y)*dt + d*dW[i]
             - (m_dY[0]*u - m_dY[1]*v - m_dY[2]*w)*dt;
          // Compute instantaneous values derived from updated Y
          derived( particles, p, i );
        }
      }
    }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const char m_depvar;                //!< Dependent variable
    const ncomp_t m_ncomp;              //!< Number of components
    const ncomp_t m_offset;             //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator
    const ctr::DepvarType m_solve;      //!< Depndent variable to solve for

    const bool m_velocity_coupled;      //!< True if coupled to velocity
    const char m_velocity_depvar;       //!< Coupled velocity dependent variable
    const ncomp_t m_velocity_offset;    //!< Offset of coupled velocity eq
    //! Quantity the coupled velocity eq solves for
    const ctr::DepvarType m_velocity_solve;

    const bool m_dissipation_coupled;   //!< True if coupled to dissipation
    const char m_dissipation_depvar;    //!< Depvar of coupled dissipation eq
    const ncomp_t m_dissipation_offset; //!< Offset of coupled dissipation eq

    std::array< tk::real, 3 > m_dY;     //! Prescribed mean scalar gradient

    //! Coefficients
    std::vector< kw::sde_bprime::info::expect::type > m_bprime;
    std::vector< kw::sde_S::info::expect::type > m_S;
    std::vector< kw::sde_kappaprime::info::expect::type > m_kprime;
    std::vector< kw::sde_rho2::info::expect::type > m_rho2;
    std::vector< kw::sde_r::info::expect::type > m_r;
    std::vector< kw::sde_b::info::expect::type > m_b;
    std::vector< kw::sde_kappa::info::expect::type > m_k;

    //! Coefficients policy
    Coefficients coeff;

    //! Selected inverse hydrodynamics time scales (if used) for each component
    //! \details This is only used if the coefficients policy is
    //!   MixMassFracBetaCoeffHydroTimeScale. See constructor.
    std::vector< tk::Table > m_hts;

    //! Selected hydrodynamics production/dissipation (if used) for each comp.
    //! \details This is only used if the coefficients policy is
    //!   MixMassFracBetaCoeffHydroTimeScale. See constructor.
    std::vector< tk::Table > m_hp;

    //! \brief Return density for mass fraction
    //! \details Functional wrapper around the dependent variable of the beta
    //!   SDE. This function returns the instantaneous density, rho,
    //!   based on the mass fraction, Y, and parameters rho2 and r'.
    //! \param[in] Y Instantaneous value of the mass fraction, Y
    //! \param[in] i Index specifying which (of multiple) parameters to use
    //! \return Instantaneous value of the density, rho
    tk::real rho( tk::real Y, ncomp_t i ) const {
      return m_rho2[i] / ( 1.0 + m_r[i] * Y );
    }

    //! \brief Return specific volume for mass fraction
    //! \details Functional wrapper around the dependent variable of the beta
    //!   SDE. This function returns the instantaneous specific volume, V,
    //!   based on the mass fraction, Y, and parameters rho2 and r'.
    //! \param[in] Y Instantaneous value of the mass fraction, Y
    //! \param[in] i Index specifying which (of multiple) parameters to use
    //! \return Instantaneous value of the specific volume, V
    tk::real vol( tk::real Y, ncomp_t i ) const {
      return ( 1.0 + m_r[i] * Y ) / m_rho2[i];
    }

    //! Compute instantaneous values derived from updated Y
    //! \param[in,out] particles Particle properties array
    //! \param[in] p Particle index
    //! \param[in] i Component index
    void derived( tk::Particles& particles, ncomp_t p, ncomp_t i ) const {
      tk::real& Y = particles( p, i, m_offset );
      particles( p, m_ncomp+i, m_offset ) = rho( Y, i );
      particles( p, m_ncomp*2+i, m_offset ) = vol( Y, i );
      particles( p, m_ncomp*3+i, m_offset ) = 1.0 - Y;
    }

    //! Initialize imposed mean scalar gradient from user input
    std::array< tk::real, 3 > initScalarGradient() {
      const auto& mg = g_inputdeck.get< tag::param, eq, tag::mean_gradient >();
      std::array< tk::real, 3 > dY{{ 0.0, 0.0, 0.0 }};
      if (mg.size() > m_c) {
        const auto& g = mg[ m_c ];
        dY = {{ g[0], g[1], g[2] }};
      }
      Assert( dY.size() == 3, "Mean scalar gradient vector size must be 3" );
      return dY;
    }
};

} // walker::

#endif // MixMassFractionBeta_h
