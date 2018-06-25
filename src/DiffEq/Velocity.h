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
      m_position_offset(
        g_inputdeck.get< tag::param, tag::velocity, tag::position_id >().at(c)),
      m_dissipation_depvar(
        g_inputdeck.get< tag::param, tag::velocity, tag::dissipation >().at(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::velocity, tag::rng >().at(c) ) ) ),
      m_U( {{ tk::ctr::mean( m_depvar, 0 ),
              tk::ctr::mean( m_depvar, 1 ),
              tk::ctr::mean( m_depvar, 2 ) }} ),
      m_c0(),
      m_g(),
      m_coeff( g_inputdeck.get< tag::param, tag::velocity, tag::c0 >().at(c),
               m_c0 )
    {
      Assert( m_ncomp == 3, "Velocity eq number of components must be 3" );
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
      using tk::ctr::lookup;

      // Update coefficients
      tk::real eps = 0.0;
      m_coeff.update( m_depvar, m_dissipation_depvar, moments, m_hts, m_c0, t,
                      eps, m_g );

      // Access mean velocity
      std::array< tk::real, 3 > U{{ lookup(m_U[0],moments),
                                    lookup(m_U[1],moments),
                                    lookup(m_U[2],moments) }};

      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        std::vector< tk::real > dW( m_ncomp );
        m_rng.gaussian( stream, m_ncomp, dW.data() );
        // Acces particle velocity
        std::array< const tk::real*, 3 > Up{{ particles.cptr( 0, m_offset ),
                                              particles.cptr( 1, m_offset ),
                                              particles.cptr( 2, m_offset ) }};
        tk::real d = m_c0 * eps * dt;
        d = (d > 0.0 ? std::sqrt(d) : 0.0);
        tk::real u = particles.var(Up[0],p) - U[0];
        tk::real v = particles.var(Up[1],p) - U[1];
        tk::real w = particles.var(Up[2],p) - U[2];
        particles.var(Up[0],p) += (m_g[0]*u + m_g[1]*v + m_g[2]*w)*dt + d*dW[0];
        particles.var(Up[1],p) += (m_g[3]*u + m_g[4]*v + m_g[5]*w)*dt + d*dW[1];
        particles.var(Up[2],p) += (m_g[6]*u + m_g[7]*v + m_g[8]*w)*dt + d*dW[2];
      }
    }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const char m_depvar;                //!< Dependent variable
    const ncomp_t m_ncomp;              //!< Number of components
    const ncomp_t m_offset;             //!< Offset SDE operates from
    const ncomp_t m_position_offset;    //!< Offset of coupled position eq
    const char m_dissipation_depvar;    //!< Depvar of coupled dissipation eq
    const tk::RNG& m_rng;               //!< Random number generator
    //! Array of tk::ctr::Product used to access the mean velocity
    const std::array< tk::ctr::Product, 3 > m_U;

    //! Selected inverse hydrodynamics time scale (if used)
    //! \details This is only used if the coefficients policy is
    //!   VelocityCoeffHydroTimeScale. See constructor.
    tk::Table m_hts;

    //! Coefficients
    kw::sde_c0::info::expect::type m_c0;
    std::array< tk::real, 9 > m_g;

    //! Coefficients policy
    Coefficients m_coeff;
};

} // walker::

#endif // Velocity_h
