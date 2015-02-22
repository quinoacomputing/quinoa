//******************************************************************************
/*!
  \file      src/DiffEq/FunctionalBeta.h
  \author    J. Bakosi
  \date      Sun 22 Feb 2015 11:45:48 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     System of functional beta SDEs
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs) with linear drift and quadratic diagonal
    diffusion, whose invariant is the joint [beta
    distribution](http://en.wikipedia.org/wiki/Beta_distribution). The main
    difference compared to the plain beta SDE (see DiffEq/Beta.h), is that in
    the functional beta SDE the dependent variable, in addition to being
    integrated in the same exact way as in the plain beta SDE, is also an
    independent variable of additional mathematical functions. Effectively, this
    is a beta SDE but there are additinal functional wrappers around the
    dependent variable governed by the beta SDE. In other words, if Y is
    governed by the beta SDE, then the functional beta SDE also governs X =
    f(Y), where both Y and X are random variables. The functions are currently
    hard-coded, but eventually may be abstracted away to allow for different
    mathematical functions.

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

    In addition to integrating the above SDE, currently there are two additional
    functions of \f$ Y_\alpha \f$ are computed as
    \f[ \begin{align}
      \rho(Y_\alpha) & = \rho_{2\alpha} ( 1 - r'_\alpha Y_\alpha ) \\
      V(Y_\alpha) & = \frac{1}{ \rho_{2\alpha} ( 1 - r'_\alpha Y_\alpha ) }
    \end{align} \f]
    These equations compute the instantaneous mixture density, \f$ \rho \f$, and
    instantaneous specific volume, \f$ V_\alpha \f$, for equation \f$ \alpha \f$
    in the system. These quantities are used in binary mixing of
    variable-density turbulence between two fluids with constant densities, \f$
    \rho_1, \f$ and \f$ \rho_2 \f$. The additional parameters, \f$ \rho_2 \f$
    and \f$ r' \f$ are user input parameters and kept constant during
    integration.

    _All of this is unpublished, but will be linked in here once published_.
*/
//******************************************************************************
#ifndef FunctionalBeta_h
#define FunctionalBeta_h

#include <cmath>

#include <InitPolicy.h>
#include <FunctionalBetaCoeffPolicy.h>
#include <RNG.h>

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief FunctionalBeta SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see
//!     DiffEq/FunctionalBetaCoeffPolicy.h
template< class Init, class Coefficients >
class FunctionalBeta {

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of functional beta SDEs to
    //!   construct. There can be multiple funcbeta ... end blocks in a control
    //!   file. This index specifies which functional beta SDE system to
    //!   instantiate. The index corresponds to the order in which the funcbeta
    //!   ... end blocks are given the control file.
    //! \author J. Bakosi
    explicit FunctionalBeta( tk::ctr::ncomp_type c ) :
      m_depvar(
        g_inputdeck.get< tag::param, tag::funcbeta, tag::depvar >().at(c) ),
      m_ncomp(
        g_inputdeck.get< tag::component >().get< tag::funcbeta >().at(c) / 3 ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< tag::funcbeta >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::funcbeta, tag::rng >().at(c) ) ) ),
      coeff( m_ncomp,
             g_inputdeck.get< tag::param, tag::funcbeta, tag::b >().at(c),
             g_inputdeck.get< tag::param, tag::funcbeta, tag::S >().at(c),
             g_inputdeck.get< tag::param, tag::funcbeta, tag::kappa >().at(c),
             g_inputdeck.get< tag::param, tag::funcbeta, tag::rho2 >().at(c),
             g_inputdeck.get< tag::param, tag::funcbeta, tag::rcomma >().at(c),
             m_b, m_S, m_k, m_rho2, m_rcomma ) {}

    //! Initalize SDE, prepare for time integration
    //! \param[inout] particles Array of particle properties 
    //! \param[in] stat Statistics object for accessing moments 
    //! \author J. Bakosi
    void initialize( tk::ParProps& particles, const tk::Statistics& stat ) {
      //! Set initial conditions using initialization policy
      Init( { particles } );
      //! Pre-lookup positions of statistical moments required by the
      //! coefficients policy
      coeff.lookup( stat, m_depvar, m_ncomp );
    }

    //! \brief Advance particles according to the system of functional beta SDEs
    //! \author J. Bakosi
    void advance( tk::ParProps& particles, int stream, tk::real dt ) {
      // Update SDE coefficients
      coeff.update( m_b, m_S, m_k, m_rho2, m_rcomma );
      // Advance particles
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        tk::real dW[m_ncomp];
        m_rng.gaussian( stream, m_ncomp, dW );
        // Advance all m_ncomp scalars
        for (tk::ctr::ncomp_type i=0; i<m_ncomp; ++i) {
          tk::real& X = particles( p, i, m_offset );
          tk::real d = m_k[i] * X * (1.0 - X) * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          X += 0.5*m_b[i]*(m_S[i] - X)*dt + d*dW[i];
          // Compute instantaneous values derivedderived  from updated X
          particles( p, m_ncomp+i, m_offset ) = rho( X, i );
          particles( p, m_ncomp*2+i, m_offset ) = vol( X, i );
        }
      }
    }

  private:
    const char m_depvar;                //!< Dependent variable
    const tk::ctr::ncomp_type m_ncomp;  //!< Number of components
    const int m_offset;                 //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator

    //! Coefficients
    std::vector< kw::sde_b::info::expect::type > m_b;
    std::vector< kw::sde_S::info::expect::type > m_S;
    std::vector< kw::sde_kappa::info::expect::type > m_k;
    std::vector< kw::sde_rho2::info::expect::type > m_rho2;
    std::vector< kw::sde_rcomma::info::expect::type > m_rcomma;

    //! Coefficients policy
    Coefficients coeff;

    //! \brief Return density for mole fraction
    //! \details Functional wrapper around the dependent variable of the beta
    //!   SDE. This function returns the instantaneous density, rho,
    //!   based on the number fraction, X, and parameters rho2 and r'.
    //! \param[in] p Instantaneous value of the mole fraction, X
    //! \param[in] i Index specifying which (of multiple) parameters to use
    //! \return Instantaneous value of the density, rho
    tk::real rho( tk::real X, tk::ctr::ncomp_type i ) const {
      return m_rho2[i] * ( 1.0 - m_rcomma[i] * X );
    }

    //! \brief Return specific volume for mole fraction
    //! \details Functional wrapper around the dependent variable of the beta
    //!   SDE. This function returns the instantaneous specific volume, V,
    //!   based on the number fraction, X, and parameters rho2 and r'.
    //! \param[in] p Instantaneous value of the mole fraction, X
    //! \param[in] i Index specifying which (of multiple) parameters to use
    //! \return Instantaneous value of the specific volume, V
    tk::real vol( tk::real X, tk::ctr::ncomp_type i ) const {
      return 1.0 / m_rho2[i] / ( 1.0 - m_rcomma[i] * X );
    }
};

} // walker::

#endif // FunctionalBeta_h
