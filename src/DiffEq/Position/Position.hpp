// *****************************************************************************
/*!
  \file      src/DiffEq/Position/Position.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     A position model for Lagrangian particles
  \details   This file implements the time integration of a system of
    deterministic or stochastic differential equations to model fluctuating
    particle postitions.
*/
// *****************************************************************************
#ifndef Position_h
#define Position_h

#include <array>

#include "InitPolicy.hpp"
#include "PositionCoeffPolicy.hpp"
#include "RNG.hpp"
#include "Particles.hpp"
#include "CoupledEq.hpp"

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief Position equation used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see DiffEq/PositionCoeffPolicy.h
template< class Init, class Coefficients >
class Position {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::position;

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of position SDEs to construct
    //!   There can be multiple position ... end blocks in a control file. This
    //!   index specifies which position SDE system to instantiate. The index
    //!   corresponds to the order in which the position ... end blocks are
    //!   given the control file.
    explicit Position( ncomp_t c ) :
      m_c( c ),
      m_depvar( g_inputdeck.get< tag::param, eq, tag::depvar >().at(c) ),
      m_ncomp( g_inputdeck.get< tag::component, eq >().at(c) ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< eq >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, eq, tag::rng >().at(c) ) ) ),
      m_velocity_coupled( coupled< eq, tag::velocity >( c ) ),
      m_velocity_depvar( depvar< eq, tag::velocity >( c ) ),
      m_velocity_offset( offset< eq, tag::velocity, tag::velocity_id >( c ) ),
      m_coeff( m_dU )
    {
      // Zero prescribed mean velocity gradient if full variable is solved for
      if (g_inputdeck.get< tag::param, eq, tag::solve >().at(c) ==
            ctr::DepvarType::FULLVAR) {
        m_dU.fill( 0.0 );
      }
      Assert( m_ncomp == 3, "Position eq number of components must be 3" );
    }

    //! Initalize SDE, prepare for time integration
    //! \param[in] stream Thread (or more precisely stream) ID
    //! \param[in,out] particles Array of particle properties
    void initialize( int stream, tk::Particles& particles ) {
      //! Set initial conditions using initialization policy
      Init::template init< eq >
        ( g_inputdeck, m_rng, stream, particles, m_c, m_ncomp, m_offset );
    }

    //! \brief Advance particles according to the system of position SDEs
    //! \param[in,out] particles Array of particle properties
    //! \param[in] dt Time step size
    void advance( tk::Particles& particles,
                  int,
                  tk::real dt,
                  tk::real,
                  const std::map< tk::ctr::Product, tk::real >& )
    {
      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Access particle velocity
        tk::real u = particles( p, 0, m_velocity_offset );
        tk::real v = particles( p, 1, m_velocity_offset );
        tk::real w = particles( p, 2, m_velocity_offset );
        // Advance all particle positions
        tk::real& Xp = particles( p, 0, m_offset );
        tk::real& Yp = particles( p, 1, m_offset );
        tk::real& Zp = particles( p, 2, m_offset );
        // Advance particle position
        Xp += (m_dU[0]*Xp + m_dU[1]*Yp + m_dU[2]*Zp + u)*dt;
        Yp += (m_dU[3]*Xp + m_dU[4]*Yp + m_dU[5]*Zp + v)*dt;
        Zp += (m_dU[6]*Xp + m_dU[7]*Yp + m_dU[8]*Zp + w)*dt;
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

    //! (Optionally) prescribed mean velocity gradient
    std::array< tk::real, 9 > m_dU;
};

} // walker::

#endif // Position_h
