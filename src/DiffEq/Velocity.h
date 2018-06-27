// *****************************************************************************
/*!
  \file      src/DiffEq/Velocity.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
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

#include "InitPolicy.h"
#include "VelocityCoeffPolicy.h"
#include "RNG.h"
#include "Particles.h"
#include "SystemComponents.h"

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
    using ncomp_t = tk::ctr::ncomp_type;

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of beta SDEs to construct.
    //!   There can be multiple beta ... end blocks in a control file. This
    //!   index specifies which beta SDE system to instantiate. The index
    //!   corresponds to the order in which the beta ... end blocks are given
    //!   the control file.
    explicit Velocity( ncomp_t c ) :
      m_c( c ),
      m_depvar(
        g_inputdeck.get< tag::param, tag::velocity, tag::depvar >().at(c) ),
      m_ncomp(
        g_inputdeck.get< tag::component >().get< tag::velocity >().at(c) ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< tag::velocity >(c) ),
      m_dissipation_depvar(
        g_inputdeck.get< tag::param, tag::velocity, tag::dissipation >().at(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::velocity, tag::rng >().at(c) ) ) ),
      m_solve(g_inputdeck.get< tag::param, tag::velocity, tag::solve >().at(c)),
      m_U( {{ tk::ctr::mean( m_depvar, 0 ),
              tk::ctr::mean( m_depvar, 1 ),
              tk::ctr::mean( m_depvar, 2 ) }} ),
      m_variant(
       g_inputdeck.get< tag::param, tag::velocity, tag::variant >().at(c) ),
      m_c0(),
      m_G(),
      m_coeff( g_inputdeck.get< tag::param, tag::velocity, tag::c0 >().at(c),
               m_c0,
               m_dU )
    {
      Assert( m_ncomp == 3, "Velocity eq number of components must be 3" );
      // Zero prescribed mean velocity gradient if full variable is solved for
      if (m_solve == ctr::DepvarType::FULLVAR) m_dU.fill( 0.0 );
      // Populate inverse hydrodynamics time scales extracted from DNS
      if ( Coefficients::type() == ctr::CoeffPolicyType::HYDROTIMESCALE ) {
        // Configure inverse hydrodyanmics time scale from DNS
        const auto& hts = g_inputdeck.get< tag::param,
                                           tag::velocity,
                                           tag::hydrotimescales >().at(c);
        Assert( hts.size() == 1,
                "Velocity eq Hydrotimescales vector size must be 1" );
        m_hts = ctr::HydroTimeScales().table( hts[0] );
      }
    }

    //! Initalize SDE, prepare for time integration
    //! \param[in] stream Thread (or more precisely stream) ID
    //! \param[in,out] particles Array of particle properties
    void initialize( int stream, tk::Particles& particles ) {
      // Set initial conditions using initialization policy
      Init::template
        init< tag::velocity >
            ( g_inputdeck, m_rng, stream, particles, m_c, m_ncomp, m_offset );
    }

    //! \brief Advance particles according to the system of beta SDEs
    //! \param[in,out] particles Array of particle properties
    //! \param[in] stream Thread (or more precisely stream) ID
    //! \param[in] dt Time step size
    void advance( tk::Particles& particles,
                  int stream,
                  tk::real dt,
                  tk::real t,
                  const std::map< tk::ctr::Product, tk::real >& moments )
    {
      // Update coefficients
      tk::real eps = 0.0;
      m_coeff.update( m_depvar, m_dissipation_depvar, moments, m_hts, m_solve,
                      m_variant, m_c0, t, eps, m_G );

      // Access mean velocity (if needed)
      std::array< tk::real, 3 > U{{ 0.0, 0.0, 0.0 }};
      if (m_solve == ctr::DepvarType::FULLVAR) {
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
        // Acces particle velocity
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
      }
    }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const char m_depvar;                //!< Dependent variable
    const ncomp_t m_ncomp;              //!< Number of components
    const ncomp_t m_offset;             //!< Offset SDE operates from
    const char m_dissipation_depvar;    //!< Depvar of coupled dissipation eq
    const tk::RNG& m_rng;               //!< Random number generator
    const ctr::DepvarType m_solve;      //!< Depndent variable to solve for
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
};

} // walker::

#endif // Velocity_h
