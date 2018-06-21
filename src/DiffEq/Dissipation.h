// *****************************************************************************
/*!
  \file      src/DiffEq/Dissipation.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
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
#include "SystemComponents.h"

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

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of beta SDEs to construct.
    //!   There can be multiple beta ... end blocks in a control file. This
    //!   index specifies which beta SDE system to instantiate. The index
    //!   corresponds to the order in which the beta ... end blocks are given
    //!   the control file.
    explicit Dissipation( ncomp_t c ) :
      m_c( c ),
      m_depvar(
        g_inputdeck.get< tag::param, tag::dissipation, tag::depvar >().at(c) ),
      m_ncomp(
        g_inputdeck.get< tag::component >().get< tag::dissipation >().at(c) ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< tag::dissipation >(c) ),
      m_velocity_offset( g_inputdeck.get< tag::param, tag::dissipation,
                                          tag::velocity_id >().at(c)),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::dissipation, tag::rng >().at(c) ) ) )
    {
      Assert( m_ncomp == 1, "Dissipation eq number of components must be 1" );
    }

    //! Initalize SDE, prepare for time integration
    //! \param[in] stream Thread (or more precisely stream) ID
    //! \param[in,out] particles Array of particle properties
    void initialize( int stream, tk::Particles& particles ) {
      //! Set initial conditions using initialization policy
      Init::template
        init< tag::dissipation >
            ( g_inputdeck, m_rng, stream, particles, m_c, m_ncomp, m_offset );
    }

    //! \brief Advance particles according to the system of beta SDEs
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
        tk::real Up = particles( p, 0, m_velocity_offset );
        tk::real Vp = particles( p, 1, m_velocity_offset );
        tk::real Wp = particles( p, 2, m_velocity_offset );
        // Advance all particle frequency
        tk::real& Op = particles( p, 0, m_offset );
        //Op += Op*dt;
      }
    }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const char m_depvar;                //!< Dependent variable
    const ncomp_t m_ncomp;              //!< Number of components
    const ncomp_t m_offset;             //!< Offset SDE operates from
    const ncomp_t m_velocity_offset;    //!< Offset for coupled velocity eq
    const tk::RNG& m_rng;               //!< Random number generator
};

} // walker::

#endif // Dissipation_h
