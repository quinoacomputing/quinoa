// *****************************************************************************
/*!
  \file      src/DiffEq/Dirichlet/MixDirichlet.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mixture Dirichlet SDE
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs), whose invariant is the [Dirichlet
    distribution](http://en.wikipedia.org/wiki/Dirichlet_distribution), with
    various constraints to model multi-material mixing process in turbulent
    flows.

    In a nutshell, the equation integrated governs a set of scalars,
    \f$0\!\le\!Y_\alpha\f$, \f$\alpha\!=\!1,\dots,N-1\f$,
    \f$\sum_{\alpha=1}^{N-1}Y_\alpha\!\le\!1\f$, as

    @m_class{m-show-m}

    \f[
       \mathrm{d}Y_\alpha(t) = \frac{b_\alpha}{2}\big[S_\alpha Y_N -
       (1-S_\alpha)Y_\alpha\big]\mathrm{d}t + \sqrt{\kappa_\alpha Y_\alpha
       Y_N}\mathrm{d}W_\alpha(t), \qquad \alpha=1,\dots,N-1
    \f]

    @m_class{m-hide-m}

    \f[
       \begin{split}
       \mathrm{d}Y_\alpha(t) = \frac{b_\alpha}{2}\big[S_\alpha Y_N -
       (1-S_\alpha)Y_\alpha\big]\mathrm{d}t + \sqrt{\kappa_\alpha Y_\alpha
       Y_N}\mathrm{d}W_\alpha(t), \\ \alpha=1,\dots,N-1
       \end{split}
    \f]

    with parameter vectors \f$b_\alpha > 0\f$, \f$\kappa_\alpha > 0\f$, and \f$0
    < S_\alpha < 1\f$, and \f$Y_N = 1-\sum_{\alpha=1}^{N-1}Y_\alpha\f$. Here
    \f$\mathrm{d}W_\alpha(t)\f$ is an isotropic vector-valued [Wiener
    process](http://en.wikipedia.org/wiki/Wiener_process) with independent
    increments. The invariant distribution is the Dirichlet distribution,
    provided the parameters of the drift and diffusion terms satisfy
    \f[
      (1-S_1) b_1 / \kappa_1 = \dots = (1-S_{N-1}) b_{N-1} / \kappa_{N-1}.
    \f]
    To keep the invariant distribution Dirichlet, the above constraint on the
    coefficients must be satisfied. For more details on the Dirichlet SDE, see
    https://doi.org/10.1155/2013/842981.
*/
// *****************************************************************************
#ifndef MixDirichlet_h
#define MixDirichlet_h

#include <vector>
#include <cmath>

#include "InitPolicy.hpp"
#include "MixDirichletCoeffPolicy.hpp"
#include "RNG.hpp"
#include "Particles.hpp"

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! Number of derived variables computed by the MixDirichlet SDE
const std::size_t MIXDIR_NUMDERIVED = 2;

//! \brief MixDirichlet SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see DiffEq/DirCoeffPolicy.h
template< class Init, class Coefficients >
class MixDirichlet {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::mixdirichlet;

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which MixDirichlet SDE to construct. There
    //!   can be multiple dirichlet ... end blocks in a control file. This index
    //!   specifies which MixDirichlet SDE to instantiate. The index corresponds
    //!   to the order in which the dirichlet ... end blocks are given the
    //!   control file.
    explicit MixDirichlet( ncomp_t c ) :
      m_c( c ),
      m_depvar( g_inputdeck.get< tag::param, eq, tag::depvar >().at(c) ),
      // subtract the number of derived variables computed, see advance()
      m_ncomp( g_inputdeck.get< tag::component >().get< eq >().at(c) -
               MIXDIR_NUMDERIVED ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< eq >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, eq, tag::rng >().at(c) ) ) ),
      m_b(),
      m_S(),
      m_kprime(),
      m_k(),
      m_rho(),
      m_r(),
      coeff( m_ncomp,
             g_inputdeck.get< tag::param, eq, tag::b >().at(c),
             g_inputdeck.get< tag::param, eq, tag::S >().at(c),
             g_inputdeck.get< tag::param, eq, tag::kappaprime >().at(c),
             g_inputdeck.get< tag::param, eq, tag::rho >().at(c),
             m_b, m_S, m_kprime, m_rho, m_r, m_k ) {}

    //! Initalize SDE, prepare for time integration
    //! \param[in] stream Thread (or more precisely stream) ID 
    //! \param[in,out] particles Array of particle properties 
    void initialize( int stream, tk::Particles& particles ) {
      //! Set initial conditions using initialization policy
      Init::template init< eq >( g_inputdeck, m_rng, stream, particles, m_c,
                                 m_ncomp, m_offset );

    tk::real eps = std::numeric_limits<tk::real>::epsilon() * 100.0;

//std::cout << "rho: " << m_rho[0] << ", " << m_rho[1] << ", " << m_rho[2] << '\n';
//std::cout << "r: " << m_r[0] << ", " << m_r[1] << ", " << m_r[2] << '\n';
//std::cout << "Y: ";
      // Initialize derived instantaneous variables
      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // start computing specific volume
        tk::real v = -1.0;
        // start computing density
        tk::real rho = 0.0;
        for (ncomp_t i=0; i<m_ncomp; ++i) {
          auto Y = particles( p, i, m_offset );
//std::cout << Y << ", ";
          if (Y < 0.0 || Y > 1.0) {
            std::cout << "IC Y out of bounds: " << Y << '\n';
          }
          v += m_r[i] * Y;
          rho += Y / m_rho[i];
        }
        // compute Nth scalar
        tk::real Yn = 1.0 - particles( p, 0, m_offset );
        for (ncomp_t i=1; i<m_ncomp; ++i) Yn -= particles( p, i, m_offset );
//std::cout << Yn << '\n';
        // finish computing specific volume
        v += m_r[m_ncomp] * Yn;
        v /= m_rho[m_ncomp];
        // finish computing density
        rho += Yn / m_rho[m_ncomp];
        rho = 1.0 / rho;
        if (std::abs(rho-1.0/v) > eps)
          std::cout << "IC: rho != 1/v, rho = " << rho << ", v = " << v
                    << ", rho-1/v = " << rho - 1.0/v << '\n';
        // Compute and store instantaneous density
        particles( p, m_ncomp, m_offset ) = rho; //1.0/v;
        // Store instantaneous specific volume
        particles( p, m_ncomp+1, m_offset ) = 1.0/rho; //v;
      }
    }

    //! \brief Advance particles according to the MixDirichlet SDE
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
      tk::real eps = std::numeric_limits<tk::real>::epsilon() * 100.0;

      // Update SDE coefficients
      coeff.update( m_depvar, m_ncomp, moments, m_rho, m_r, m_kprime, m_b, m_k,
                    m_S );
      // Advance particles
      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        std::vector< tk::real > dW( m_ncomp );
        m_rng.gaussian( stream, m_ncomp, dW.data() );
        // compute Nth scalar
        tk::real Yn = 1.0 - particles(p, 0, m_offset);
        for (ncomp_t i=1; i<m_ncomp; ++i) Yn -= particles( p, i, m_offset );
        // start computing specific volume
        tk::real v = -1.0;
        // start computing density
        tk::real rho = 0.0;
        // Advance first m_ncomp (K=N-1) scalars
        for (ncomp_t i=0; i<m_ncomp; ++i) {
          tk::real& Y = particles( p, i, m_offset );
          tk::real d = m_k[i] * Y * Yn * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          Y += 0.5*m_b[i]*( m_S[i]*Yn - (1.0-m_S[i]) * Y )*dt + d*dW[i];
          if (Y < eps) Y = eps;
          if (Y > 1.0-eps) Y = 1.0-eps;
          v += m_r[i] * Y;
          rho += Y / m_rho[i];
        }
        // recompute Nth scalar after time step just taken
        Yn = 1.0 - particles(p, 0, m_offset);
        for (ncomp_t i=1; i<m_ncomp; ++i) Yn -= particles( p, i, m_offset );
        // Finish computing specific volume
        v += m_r[m_ncomp] * Yn;
        v /= m_rho[m_ncomp];
        // Finish computing density
        rho += Yn / m_rho[m_ncomp];
        rho = 1.0 / rho;
        if (std::abs(rho-1.0/v) > eps)
          std::cout << "rho != 1/v, rho = " << rho << ", v = " << v
                    << ", rho-1/v = " << rho - 1.0/v << '\n';
        // Compute and store instantaneous density
        particles( p, m_ncomp, m_offset ) = rho; //1.0/v;
        // Store instantaneous specific volume
        particles( p, m_ncomp+1, m_offset ) = 1.0/rho; //v;
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
    std::vector< kw::sde_kappa::info::expect::type > m_kprime;
    std::vector< kw::sde_kappa::info::expect::type > m_k;
    std::vector< kw::sde_rho::info::expect::type > m_rho;
    std::vector< kw::sde_r::info::expect::type > m_r;

    //! Coefficients policy
    Coefficients coeff;
};

} // walker::

#endif // MixDirichlet_h
