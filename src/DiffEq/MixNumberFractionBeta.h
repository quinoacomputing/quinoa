// *****************************************************************************
/*!
  \file      src/DiffEq/MixNumberFractionBeta.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     System of mix number-fraction beta SDEs
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs) with linear drift and quadratic diagonal
    diffusion, whose invariant is the joint [beta
    distribution](http://en.wikipedia.org/wiki/Beta_distribution). There are two
    differences compared to the plain beta SDE (see DiffEq/Beta.h):

    - First, the parameters, b, and kappa are specified via functions that
    constrain the beta SDE to be consistent with the turbulent mixing process.
    In particular, the SDE is made consistent with the no-mix and fully mixed
    limits. See, e.g., MixNumberFractionBetaCoeffConst::update().

    - Second, there two additional random variables computed, the same as also
    computed by the number-fraction beta equation, see also
    DiffEq/NumberFractionBeta.h.

    In a nutshell, the equation integrated governs a set of scalars,
    \f$0\!\le\!X_\alpha\f$, \f$\alpha\!=\!1,\dots,N\f$, as
    \f[
       \mathrm{d}X_\alpha(t) = \frac{b_\alpha}{2}\left(S_\alpha - X_\alpha\right)
       \mathrm{d}t + \sqrt{\kappa_\alpha X_\alpha(1-X_\alpha)}
       \mathrm{d}W_\alpha(t), \qquad \alpha=1,\dots,N
    \f]
    with parameter vectors \f$b_\alpha = \Theta b'_\alpha > 0\f$, \f$
    \newcommand{\irv}[1]{\langle{#1^2}\rangle} \kappa_\alpha = \kappa' \irv{x} >
    0\f$, and \f$0 < S_\alpha < 1\f$. This is similar to DiffEq/Beta.h, but the
    parameters, \f$b\f$ and \f$\kappa\f$ constrained. Here \f$
    \newcommand{\irv}[1]{\langle{#1^2}\rangle}
    \newcommand{\irmean}[1]{{\langle{#1}\rangle}} \Theta = 1 - \irv{x} /
    [ \irmean{X} (1-\irmean{X}) ]\f$. The fluctuation about the mean, \f$
    \newcommand{\irmean}[1]{{\langle{#1}\rangle}} \irmean{X} \f$, is defined as
    usual: \f$ \newcommand{\irmean}[1]{{\langle{#1}\rangle}} x = X - \irmean{X}
    \f$, and \f$b'\f$ and \f$ \kappa'\f$ are user-specified constants. Also,
    \f$\mathrm{d}W_\alpha(t)\f$ is an isotropic vector-valued
    [Wiener process](http://en.wikipedia.org/wiki/Wiener_process) with
    independent increments. The invariant distribution is the joint beta
    distribution. This system of SDEs consists of N independent equations. For
    more on the beta SDE, see http://dx.doi.org/10.1080/14685248.2010.510843.

    Similar to the number-fraction beta SDE (DiffEq/NumberFractionBeta.h), in
    addition to integrating the above SDE, there are two additional functions
    of \f$ X_\alpha \f$ are computed as
    \f[ \begin{align}
      \rho(X_\alpha) & = \rho_{2\alpha} ( 1 - r'_\alpha X_\alpha ) \\
      V(X_\alpha) & = \frac{1}{ \rho_{2\alpha} ( 1 - r'_\alpha X_\alpha ) }
    \end{align} \f]
    These equations compute the instantaneous mixture density, \f$ \rho \f$, and
    instantaneous specific volume, \f$ V_\alpha \f$, for equation \f$ \alpha \f$
    in the system. These quantities are used in binary mixing of
    variable-density turbulence between two fluids with constant densities, \f$
    \rho_1, \f$ and \f$ \rho_2 \f$. The additional parameters, \f$ \rho_2 \f$
    and \f$ r' \f$ are user input parameters and kept constant during
    integration. Since we compute the above variables, \f$\rho,\f$ and \f$V\f$,
    and call them mixture density and specific volume, respectively, \f$X\f$,
    governed by the beta SDE is a number (or mole) fraction.

    _All of this is unpublished, but will be linked in here once published_.
*/
// *****************************************************************************
#ifndef MixNumberFractionBeta_h
#define MixNumberFractionBeta_h

#include <vector>
#include <cmath>

#include "InitPolicy.h"
#include "MixNumberFractionBetaCoeffPolicy.h"
#include "RNG.h"
#include "Particles.h"
#include "SystemComponents.h"

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief MixNumberFractionBeta SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see
//!     DiffEq/MixNumberFractionBetaCoeffPolicy.h
template< class Init, class Coefficients >
class MixNumberFractionBeta {

  private:
    using ncomp_t = tk::ctr::ncomp_type;

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of mix number-fraction beta
    //!   SDEs to construct. There can be multiple mixnumfracbeta ... end blocks
    //!   in a control file. This index specifies which mix number-fraction beta
    //!   SDE system to instantiate. The index corresponds to the order in which
    //!   the mixnumfracbeta ... end blocks are given the control file.
    //! \author J. Bakosi
    explicit MixNumberFractionBeta( ncomp_t c ) :
      m_c( c ),
      m_depvar(
        g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::depvar >().at(c)
      ),
      m_ncomp(
        g_inputdeck.get< tag::component >().get< tag::mixnumfracbeta >().at(c)/3
      ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< tag::mixnumfracbeta >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::rng >().at(c) ) )
      ),
      m_bprime(),
      m_S(),
      m_kprime(),
      m_rho2(),
      m_rcomma(),
      m_b(),
      m_k(),
      coeff(
        m_ncomp,
        g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::bprime >().at(c),
        g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::S >().at(c),
        g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::kappaprime >().at(c),
        g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::rho2 >().at(c),
        g_inputdeck.get< tag::param, tag::mixnumfracbeta, tag::rcomma >().at(c),
        m_bprime, m_S, m_kprime, m_rho2, m_rcomma, m_b, m_k ) {}

    //! Initalize SDE, prepare for time integration
    //! \param[in] stream Thread (or more precisely stream) ID 
    //! \param[in,out] particles Array of particle properties 
    //! \author J. Bakosi
    void initialize( int stream, tk::Particles& particles ) {
      //! Set initial conditions using initialization policy
      Init::template
        init< tag::mixnumfracbeta >
            ( g_inputdeck, m_rng, stream, particles, m_c, m_ncomp, m_offset );
    }

    //! \brief Advance particles according to the system of mix number-fraction
    //!   beta SDEs
    //! \param[in,out] particles Array of particle properties
    //! \param[in] stream Thread (or more precisely stream) ID
    //! \param[in] dt Time step size
    //! \param[in] moments Map of statistical moments
    //! \author J. Bakosi
    void advance( tk::Particles& particles,
                  int stream,
                  tk::real dt,
                  tk::real,
                  const std::map< tk::ctr::Product, tk::real >& moments )
    {
      // Update SDE coefficients
      coeff.update( m_depvar, m_ncomp, moments, m_bprime, m_kprime, m_b, m_k );
      // Advance particles
      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        std::vector< tk::real > dW( m_ncomp );
        m_rng.gaussian( stream, m_ncomp, dW.data() );
        // Advance all m_ncomp scalars
        for (ncomp_t i=0; i<m_ncomp; ++i) {
          tk::real& X = particles( p, i, m_offset );
          tk::real d = m_k[i] * X * (1.0 - X) * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          X += 0.5*m_b[i]*(m_S[i] - X)*dt + d*dW[i];
          // Compute instantaneous values derived from updated X
          particles( p, m_ncomp+i, m_offset ) = rho( X, i );
          particles( p, m_ncomp*2+i, m_offset ) = vol( X, i );
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
    std::vector< kw::sde_bprime::info::expect::type > m_bprime;
    std::vector< kw::sde_S::info::expect::type > m_S;
    std::vector< kw::sde_kappaprime::info::expect::type > m_kprime;
    std::vector< kw::sde_rho2::info::expect::type > m_rho2;
    std::vector< kw::sde_rcomma::info::expect::type > m_rcomma;
    std::vector< kw::sde_b::info::expect::type > m_b;
    std::vector< kw::sde_kappa::info::expect::type > m_k;

    //! Coefficients policy
    Coefficients coeff;

    //! \brief Return density for mole fraction
    //! \details Functional wrapper around the dependent variable of the beta
    //!   SDE. This function returns the instantaneous density, rho,
    //!   based on the number fraction, X, and parameters rho2 and r'.
    //! \param[in] X Instantaneous value of the mole fraction, X
    //! \param[in] i Index specifying which (of multiple) parameters to use
    //! \return Instantaneous value of the density, rho
    tk::real rho( tk::real X, ncomp_t i ) const {
      return m_rho2[i] * ( 1.0 - m_rcomma[i] * X );
    }

    //! \brief Return specific volume for mole fraction
    //! \details Functional wrapper around the dependent variable of the beta
    //!   SDE. This function returns the instantaneous specific volume, V,
    //!   based on the number fraction, X, and parameters rho2 and r'.
    //! \param[in] X Instantaneous value of the mole fraction, X
    //! \param[in] i Index specifying which (of multiple) parameters to use
    //! \return Instantaneous value of the specific volume, V
    tk::real vol( tk::real X, ncomp_t i ) const {
      return 1.0 / rho( X, i );
    }
};

} // walker::

#endif // MixNumberFractionBeta_h
