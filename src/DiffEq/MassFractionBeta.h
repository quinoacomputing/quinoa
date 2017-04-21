// *****************************************************************************
/*!
  \file      src/DiffEq/MassFractionBeta.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     System of mass-fraction beta SDEs
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs) with linear drift and quadratic diagonal
    diffusion, whose invariant is the joint [beta
    distribution](http://en.wikipedia.org/wiki/Beta_distribution). The main
    difference compared to the plain beta SDE (see DiffEq/Beta.h), is that in
    the mass-fraction beta SDE the dependent variable, there are two
    additional stochastic variables computed from the beta variables.

    In a nutshell, the equation integrated governs a set of scalars,
    \f$0\!\le\!Y_\alpha\f$, \f$\alpha\!=\!1,\dots,N\f$, as
    \f[
       \mathrm{d}Y_\alpha(t) = \frac{b_\alpha}{2}\left(S_\alpha - Y_\alpha\right)
       \mathrm{d}t + \sqrt{\kappa_\alpha Y_\alpha(1-Y_\alpha)}
       \mathrm{d}W_\alpha(t), \qquad \alpha=1,\dots,N
    \f]
    with parameter vectors \f$b_\alpha > 0\f$, \f$\kappa_\alpha > 0\f$, and \f$0
    < S_\alpha < 1\f$. This is the same as in DiffEq/Beta.h. Here
    \f$\mathrm{d}W_\alpha(t)\f$ is an isotropic vector-valued [Wiener
    process](http://en.wikipedia.org/wiki/Wiener_process) with independent
    increments. The invariant distribution is the joint beta distribution. This
    system of SDEs consists of N independent equations. For
    more on the beta SDE, see http://dx.doi.org/10.1080/14685248.2010.510843.

    In addition to integrating the above SDE, there are two additional functions
    of \f$ Y_\alpha \f$ are computed as
    \f[ \begin{align}
      \rho(Y_\alpha) & = \frac{ \rho_{2\alpha} }{ 1 + r_\alpha Y_\alpha } \\
      V(Y_\alpha) & = \frac{1}{ \rho(Y_\alpha) }
    \end{align} \f]
    These equations compute the instantaneous mixture density, \f$ \rho \f$, and
    instantaneous specific volume, \f$ V_\alpha \f$, for equation \f$ \alpha \f$
    in the system. These quantities are used in binary mixing of
    variable-density turbulence between two fluids with constant densities, \f$
    \rho_1, \f$ and \f$ \rho_2 \f$. The additional parameters, \f$ \rho_2 \f$
    and \f$ r' \f$ are user input parameters and kept constant during
    integration. Since we compute the above variables, \f$\rho,\f$ and \f$V\f$,
    and call them mixture density and specific volume, respectively, \f$Y\f$,
    governed by the beta SDE is a mass fraction, hence the name mass-fraction
    beta.

    _All of this is unpublished, but will be linked in here once published_.
*/
// *****************************************************************************
#ifndef MassFractionBeta_h
#define MassFractionBeta_h

#include <vector>
#include <cmath>

#include "InitPolicy.h"
#include "MassFractionBetaCoeffPolicy.h"
#include "RNG.h"
#include "Particles.h"
#include "SystemComponents.h"

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief MassFractionBeta SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see
//!     DiffEq/MassFractionBetaCoeffPolicy.h
template< class Init, class Coefficients >
class MassFractionBeta {

  private:
    using ncomp_t = tk::ctr::ncomp_type;

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of mass-fraction beta SDEs
    //!   to construct. There can be multiple massfracbeta ... end blocks in a
    //!   control file. This index specifies which mass-fraction beta SDE
    //!   system to instantiate. The index corresponds to the order in which the
    //!   massfracbeta ... end blocks are given the control file.
    //! \author J. Bakosi
    explicit MassFractionBeta( ncomp_t c ) :
      m_c( c ),
      m_depvar(
        g_inputdeck.get< tag::param, tag::massfracbeta, tag::depvar >().at(c) ),
      m_ncomp(
        g_inputdeck.get< tag::component >().get< tag::massfracbeta >().at(c) / 3 ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< tag::massfracbeta >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::massfracbeta, tag::rng >().at(c) ) ) ),
      m_b(),
      m_S(),
      m_k(),
      m_rho2(),
      m_r(),
      coeff(
        m_ncomp,
        g_inputdeck.get< tag::param, tag::massfracbeta, tag::b >().at(c),
        g_inputdeck.get< tag::param, tag::massfracbeta, tag::S >().at(c),
        g_inputdeck.get< tag::param, tag::massfracbeta, tag::kappa >().at(c),
        g_inputdeck.get< tag::param, tag::massfracbeta, tag::rho2 >().at(c),
        g_inputdeck.get< tag::param, tag::massfracbeta, tag::r >().at(c),
        m_b, m_S, m_k, m_rho2, m_r ) {}

    //! Initalize SDE, prepare for time integration
    //! \param[in] stream Thread (or more precisely stream) ID 
    //! \param[in,out] particles Array of particle properties 
    //! \author J. Bakosi
    void initialize( int stream, tk::Particles& particles ) {
      //! Set initial conditions using initialization policy
      Init::template
        init< tag::massfracbeta >
            ( g_inputdeck, m_rng, stream, particles, m_c, m_ncomp, m_offset );
    }

    //! \brief Advance particles according to the system of mass-fraction beta
    //!    SDEs
    //! \param[in,out] particles Array of particle properties
    //! \param[in] stream Thread (or more precisely stream) ID
    //! \param[in] dt Time step size
    //! \author J. Bakosi
    void advance( tk::Particles& particles,
                  int stream,
                  tk::real dt,
                  tk::real,
                  const std::map< tk::ctr::Product, tk::real >& )
    {
      // Advance particles
      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        std::vector< tk::real > dW( m_ncomp );
        m_rng.gaussian( stream, m_ncomp, dW.data() );
        // Advance all m_ncomp scalars
        for (ncomp_t i=0; i<m_ncomp; ++i) {
          tk::real& Y = particles( p, i, m_offset );
          tk::real d = m_k[i] * Y * (1.0 - Y) * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          Y += 0.5*m_b[i]*(m_S[i] - Y)*dt + d*dW[i];
          // Compute instantaneous values derived from updated Y
          particles( p, m_ncomp+i, m_offset ) = rho( Y, i );
          particles( p, m_ncomp*2+i, m_offset ) = vol( Y, i );
        }
      }
    }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const char m_depvar;                //!< Dependent variable
    const ncomp_t m_ncomp;              //!< Number of components
    const ncomp_t m_offset;             //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator

    //! Coefficients
    std::vector< kw::sde_b::info::expect::type > m_b;
    std::vector< kw::sde_S::info::expect::type > m_S;
    std::vector< kw::sde_kappa::info::expect::type > m_k;
    std::vector< kw::sde_rho2::info::expect::type > m_rho2;
    std::vector< kw::sde_r::info::expect::type > m_r;

    //! Coefficients policy
    Coefficients coeff;

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
      return 1.0 / rho( Y, i );
    }
};

} // walker::

#endif // MassFractionBeta_h
