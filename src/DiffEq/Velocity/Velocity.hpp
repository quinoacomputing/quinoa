// *****************************************************************************
/*!
  \file      src/DiffEq/Velocity/Velocity.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     A model for velocity in variable-density turbulence
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs) to model the fluctuating velocity components
    in homogeneous variable-density turbulence. This model is an extension of
    the generalized Langevin (GLM) model for constant-density flows by Haworth &
    Pope (https://doi.org/10.1063/1.865723). The extension is roughly along the
    lines of https://doi.org/10.1080/14685248.2011.554419.
*/
// *****************************************************************************
#ifndef Velocity_h
#define Velocity_h

#include <array>
#include <vector>
#include <cmath>

#include "InitPolicy.hpp"
#include "VelocityCoeffPolicy.hpp"
#include "RNG.hpp"
#include "Particles.hpp"
#include "CoupledEq.hpp"

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief Velocity SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see DiffEq/VelocityCoeffPolicy.h
template< class Init, class Coefficients >
class Velocity {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::velocity;

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of velocity SDEs to construct
    //!   There can be multiple velocity ... end blocks in a control file. This
    //!   index specifies which velocity SDE system to instantiate. The index
    //!   corresponds to the order in which the velocity ... end blocks are
    //!   given the control file.
    explicit Velocity( ncomp_t c ) :
      m_c( c ),
      m_depvar( g_inputdeck.get< tag::param, eq, tag::depvar >().at(c) ),
      m_solve( g_inputdeck.get< tag::param, eq, tag::solve >().at(c) ),
      m_mixmassfracbeta_coupled( coupled< eq, tag::mixmassfracbeta >( c ) ),
      m_mixmassfracbeta_depvar( depvar< eq, tag::mixmassfracbeta >( c ) ),
      m_mixmassfracbeta_offset(
        offset< eq, tag::mixmassfracbeta, tag::mixmassfracbeta_id >( c ) ),
      m_mixmassfracbeta_ncomp(
        // The magic number, 4, below is MixMassFractionBeta::NUMDERIVED + 1,
        // but cannot be given as such, because that would lead to circular
        // dependencies of Velocity depending on MixMassfractionBeta, and vice
        // versa.
        ncomp< eq, tag::mixmassfracbeta, tag::mixmassfracbeta_id >( c ) / 4 ),
      m_numderived( numderived() ),
      m_ncomp( g_inputdeck.get< tag::component, eq >().at(c) - m_numderived ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< eq >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, eq, tag::rng >().at(c) ) ) ),
      m_position_coupled( coupled< eq, tag::position >( c ) ),
      m_position_depvar( depvar< eq, tag::position >( c ) ),
      m_position_offset( offset< eq, tag::position, tag::position_id >( c ) ),
      m_dissipation_coupled( coupled< eq, tag::dissipation >( c ) ),
      m_dissipation_depvar( depvar< eq, tag::dissipation >( c ) ),
      m_dissipation_offset(
        offset< eq, tag::dissipation, tag::dissipation_id >( c ) ),
      m_U( {{ tk::ctr::mean( m_depvar, 0 ),
              tk::ctr::mean( m_depvar, 1 ),
              tk::ctr::mean( m_depvar, 2 ) }} ),
      m_variant( g_inputdeck.get< tag::param, eq, tag::variant >().at(c) ),
      m_c0(),
      m_G(),
      m_coeff( g_inputdeck.get< tag::param, eq, tag::c0 >().at(c), m_c0, m_dU ),
      m_gravity( { 0.0, 0.0, 0.0 } )
    {
      Assert( m_ncomp == 3, "Velocity eq number of components must be 3" );
      // Zero prescribed mean velocity gradient if full variable is solved for
      if (m_solve == ctr::DepvarType::FULLVAR) m_dU.fill( 0.0 );
      // Populate inverse hydrodynamics time scales extracted from DNS
      if ( Coefficients::type() == ctr::CoeffPolicyType::HYDROTIMESCALE ) {
        // Configure inverse hydrodyanmics time scale from DNS
        const auto& hts =
          g_inputdeck.get< tag::param, eq, tag::hydrotimescales >().at(c);
        Assert( hts.size() == 1,
                "Velocity eq Hydrotimescales vector size must be 1" );
        m_hts = ctr::HydroTimeScales().table( hts[0] );
      }
      // Initialize gravity body force if configured
      const auto& gravity =
        g_inputdeck.get< tag::param, eq, tag::gravity >().at(c);
      if (!gravity.empty()) {
        m_gravity[0] = gravity[0];
        m_gravity[1] = gravity[1];
        m_gravity[2] = gravity[2];
      }
    }

    //! Compute number of derived variables
    //! \return Number of derived variables computed
    std::size_t numderived() const {
      if (m_solve == ctr::DepvarType::PRODUCT)  // solve for momentum
        // 3 velocity components for each coupled mass fraction
        return m_mixmassfracbeta_ncomp * 3;
      else
        return 0;
    }

    //! Initalize SDE, prepare for time integration
    //! \param[in] stream Thread (or more precisely stream) ID
    //! \param[in,out] particles Array of particle properties
    void initialize( int stream, tk::Particles& particles ) {
      // Set initial conditions using initialization policy
      Init::template init< eq >
        ( g_inputdeck, m_rng, stream, particles, m_c, m_ncomp, m_offset );
    }

    //! \brief Advance particles according to the system of velocity SDEs
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
      using ctr::DepvarType;
      const auto epsilon = std::numeric_limits< tk::real >::epsilon();

      // Update coefficients
      tk::real eps = 0.0;
      m_coeff.update( m_depvar, m_dissipation_depvar, moments, m_hts, m_solve,
                      m_variant, m_c0, t, eps, m_G );

      // Access mean velocity (if needed)
      std::array< tk::real, 3 > U{{ 0.0, 0.0, 0.0 }};
      if (m_solve == DepvarType::FULLVAR || m_solve == DepvarType::PRODUCT) {
        using tk::ctr::lookup;
        U[0] = lookup( m_U[0], moments );
        U[1] = lookup( m_U[1], moments );
        U[2] = lookup( m_U[2], moments );
      }

      // Modify G with the mean velocity gradient
      for (std::size_t i=0; i<9; ++i) m_G[i] -= m_dU[i];

      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        std::vector< tk::real > dW( m_ncomp );
        m_rng.gaussian( stream, m_ncomp, dW.data() );
        // Access particle velocity
        tk::real& Up = particles( p, 0, m_offset );
        tk::real& Vp = particles( p, 1, m_offset );
        tk::real& Wp = particles( p, 2, m_offset );
        // Compute diffusion
        tk::real d = m_c0 * eps * dt;
        d = (d > 0.0 ? std::sqrt(d) : 0.0);
        // Compute velocity fluctuation
        tk::real u = Up - U[0];
        tk::real v = Vp - U[1];
        tk::real w = Wp - U[2];
        // Update particle velocity
        Up += (m_G[0]*u + m_G[1]*v + m_G[2]*w)*dt + d*dW[0];
        Vp += (m_G[3]*u + m_G[4]*v + m_G[5]*w)*dt + d*dW[1];
        Wp += (m_G[6]*u + m_G[7]*v + m_G[8]*w)*dt + d*dW[2];
        // Optionally compute particle velocities derived from particle momentum
        if (m_solve == ctr::DepvarType::PRODUCT) {      // if solve for momentum
          for (ncomp_t i=0; i<m_numderived/3; ++i) {
            auto rhoi = particles( p, m_mixmassfracbeta_ncomp+i,
                                  m_mixmassfracbeta_offset );
            if (std::abs(rhoi) > epsilon) {
              particles( p, m_ncomp+(i*3)+0, m_offset ) = Up/rhoi;
              particles( p, m_ncomp+(i*3)+1, m_offset ) = Vp/rhoi;
              particles( p, m_ncomp+(i*3)+2, m_offset ) = Wp/rhoi;
              Up += rhoi * m_gravity[0] * dt;
              Vp += rhoi * m_gravity[1] * dt;
              Wp += rhoi * m_gravity[2] * dt;
            }
          }
        } else {
          Up += m_gravity[0] * dt;
          Vp += m_gravity[1] * dt;
          Wp += m_gravity[2] * dt;
        }
      }
    }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const char m_depvar;                //!< Dependent variable
    const ctr::DepvarType m_solve;      //!< Dependent variable to solve for

    //! True if coupled to mixmassfracbeta
    const bool m_mixmassfracbeta_coupled;
    //! Depvar of coupled mixmassfracbeta eq
    const char m_mixmassfracbeta_depvar;
    //! Offset of coupled mixmassfracbeta eq
    const ncomp_t m_mixmassfracbeta_offset;
    //! Number of scalar components in coupled mixmassfracbeta eq
    const ncomp_t m_mixmassfracbeta_ncomp;

    const ncomp_t m_numderived;         //!< Number of derived variables
    const ncomp_t m_ncomp;              //!< Number of components
    const ncomp_t m_offset;             //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator

    const bool m_position_coupled;      //!< True if coupled to position
    const char m_position_depvar;       //!< Coupled position dependent variable
    const ncomp_t m_position_offset;    //!< Offset of coupled position eq

    const bool m_dissipation_coupled;   //!< True if coupled to dissipation
    const char m_dissipation_depvar;    //!< Coupled dissipation dependent var
    const ncomp_t m_dissipation_offset; //!< Offset of coupled dissipation eq

    //! Array of tk::ctr::Product used to access the mean velocity
    const std::array< tk::ctr::Product, 3 > m_U;
    //! Velocity model variant
    const ctr::VelocityVariantType m_variant;

    //! Selected inverse hydrodynamics time scale (if used)
    //! \details This is only used if the coefficients policy is
    //!   VelocityCoeffHydroTimeScale. See constructor.
    tk::Table m_hts;

    //! Coefficients
    kw::sde_c0::info::expect::type m_c0;
    std::array< tk::real, 9 > m_G;

    //! Coefficients policy
    Coefficients m_coeff;

    //! (Optionally) prescribed mean velocity gradient
    std::array< tk::real, 9 > m_dU;

    //! Optional gravity body force
    std::array< tk::real, 3 > m_gravity;
};

} // walker::

#endif // Velocity_h
