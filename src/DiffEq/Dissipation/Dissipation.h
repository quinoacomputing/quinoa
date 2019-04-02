// *****************************************************************************
/*!
  \file      src/DiffEq/Dissipation/Dissipation.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     A dissipation model for Lagrangian particles
  \details   This file implements the time integration of a system of
    stochastic differential equations to model fluctuating
    particle time scales (from which the turbulent kinetic energy dissipation
    rate can be computed).
*/
// *****************************************************************************
#ifndef Dissipation_h
#define Dissipation_h

#include <array>
#include <vector>
#include <cmath>

#include "InitPolicy.h"
#include "DissipationCoeffPolicy.h"
#include "RNG.h"
#include "Particles.h"
#include "CoupledEq.h"

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief Dissipation equation used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see DiffEq/DissipationCoeffPolicy.h
template< class Init, class Coefficients >
class Dissipation {

  private:
    using ncomp_t = tk::ctr::ncomp_type;
    using eq = tag::dissipation;

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of dissipation SDEs to
    //!   construct. There can be multiple dissipation ... end blocks in a
    //!   control file. This index specifies which dissipation SDE system to
    //!   instantiate. The index corresponds to the order in which the
    //!   dissipation ... end blocks are given the control file.
    explicit Dissipation( ncomp_t c ) :
      m_c( c ),
      m_depvar( g_inputdeck.get< tag::param, eq, tag::depvar >().at(c) ),
      m_ncomp( g_inputdeck.get< tag::component >().get< eq >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< eq >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, eq, tag::rng >().at(c) ) ) ),
      m_velocity_coupled( coupled< eq, tag::velocity >( c ) ),
      m_velocity_depvar( depvar< eq, tag::velocity >( c ) ),
      m_velocity_offset( offset< eq, tag::velocity, tag::velocity_id >( c ) ),
      m_coeff(
        g_inputdeck.get< tag::param, eq, tag::c3 >().at(c),
        g_inputdeck.get< tag::param, eq, tag::c4 >().at(c),
        g_inputdeck.get< tag::param, eq, tag::com1 >().at(c),
        g_inputdeck.get< tag::param, eq, tag::com2 >().at(c),
        m_c3, m_c4, m_com1, m_com2 ),
      m_O( tk::ctr::mean( m_depvar, 0 ) ),
      m_R( {{ tk::ctr::variance( m_velocity_depvar, 0 ),
              tk::ctr::variance( m_velocity_depvar, 1 ),
              tk::ctr::variance( m_velocity_depvar, 2 ),
              tk::ctr::covariance( m_velocity_depvar, 0, m_velocity_depvar, 1 )
           }} )
    {
      Assert( m_ncomp == 1, "Dissipation eq number of components must be 1" );
    }

    //! Initalize SDE, prepare for time integration
    //! \param[in] stream Thread (or more precisely stream) ID
    //! \param[in,out] particles Array of particle properties
    void initialize( int stream, tk::Particles& particles ) {
      //! Set initial conditions using initialization policy
      Init::template init< eq >
        ( g_inputdeck, m_rng, stream, particles, m_c, m_ncomp, m_offset );
    }

    //! \brief Advance particles according to the dissipation SDE
    //! \param[in,out] particles Array of particle properties
    //! \param[in] stream Thread (or more precisely stream) ID
    //! \param[in] dt Time step size
    //! \param[in] moments Map of statistical moments
    void advance( tk::Particles& particles,
                  int stream,
                  tk::real dt,
                  tk::real,
                  const std::map< tk::ctr::Product, tk::real >& moments )
    {
      using tk::ctr::lookup;

      // Access mean turbulence frequency
      tk::real O = lookup( m_O, moments );

      tk::ctr::Term u( static_cast<char>(std::toupper(m_velocity_depvar)), 0,
                       tk::ctr::Moment::ORDINARY );
      tk::ctr::Term v( static_cast<char>(std::toupper(m_velocity_depvar)), 1,
                       tk::ctr::Moment::ORDINARY );
      tk::ctr::Term w( static_cast<char>(std::toupper(m_velocity_depvar)), 2,
                       tk::ctr::Moment::ORDINARY );
      auto r11 = tk::ctr::Product( { u, u } );
      auto r22 = tk::ctr::Product( { v, v } );
      auto r33 = tk::ctr::Product( { w, w } );
      auto r12 = tk::ctr::Product( { u, v } );

      // Compute turbulent kinetic energy
      tk::real k = ( lookup(r11,moments) +
                     lookup(r22,moments) +
                     lookup(r33,moments) ) / 2.0;

      // Production of turbulent kinetic energy
      tk::real S = 1.0; // prescribed shear: hard-coded in a single direction
      tk::real P = -lookup(r12,moments)*S;

      // Source for turbulent frequency
      tk::real Som = m_com2 - m_com1*P/(O*k);

      // Update source based on coefficients policy
      Coefficients::src( Som );

      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate a Gaussian random number with zero mean and unit variance
        tk::real dW;
        m_rng.gaussian( stream, m_ncomp, &dW );
        // Advance particle frequency
        tk::real& Op = particles( p, 0, m_offset );
        tk::real d = 2.0*m_c3*m_c4*O*O*Op*dt;
        d = (d > 0.0 ? std::sqrt(d) : 0.0);
        Op += (-m_c3*(Op-O) - Som*Op)*O*dt + d*dW;
      }
    }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const char m_depvar;                //!< Dependent variable
    const ncomp_t m_ncomp;              //!< Number of components
    const ncomp_t m_offset;             //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator

    const bool m_velocity_coupled;      //!< True if coupled to velocity
    const char m_velocity_depvar;       //!< Coupled velocity dependent variable
    const ncomp_t m_velocity_offset;    //!< Offset of coupled velocity eq


    //! Coefficients policy
    Coefficients m_coeff;

    // Model coefficients
    tk::real m_c3;
    tk::real m_c4;
    tk::real m_com1;
    tk::real m_com2;

    //! tk::Product used to access the mean of turbulence frequency
    tk::ctr::Product m_O;
    //! Array of tk::Product used to access the the Reynolds stress
    std::array< tk::ctr::Product, 4 > m_R;
};

} // walker::

#endif // Dissipation_h
