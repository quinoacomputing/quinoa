// *****************************************************************************
/*!
  \file      src/DiffEq/GeneralizedDirichlet.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Lochner's generalized Dirichlet SDE
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs) whose invariant is Lochner's [generalized
    Dirichlet distribution](http://en.wikipedia.org/wiki/Generalized_Dirichlet_distribution).

    In a nutshell, the equation integrated governs a set of scalars, \f$0 \le
    Y_i\f$, \f$i=1,\dots,K\f$, \f$\sum_{i=1}^KY_i\le1\f$, as
    \f[ \begin{split}
      \mathrm{d}Y_i(t) = \frac{\mathcal{U}_i}{2}\left\{ b_i\Big[S_i
      \mathcal{Y}_K - (1-S_i)Y_i\Big] + Y_i\mathcal{Y}_K
      \sum_{j=i}^{K-1}\frac{c_{ij}}{\mathcal{Y}_j}\right\}\mathrm{d}t +
      \sqrt{\kappa_i Y_i \mathcal{Y}_K \mathcal{U}_i}\mathrm{d}W_i(t), \\
      \qquad i=1,\dots,K, \end{split}
    \f]
    where \f$\mathrm{d}W_i(t)\f$ is an isotropic vector-valued [Wiener
    process](http://en.wikipedia.org/wiki/Wiener_process) with independent
    increments. The statistically stationary solution of the above coupled
    system of nonlinear stochastic differential equations is the generalized
    Dirichlet distribution,
    \f[
       \newcommand{\bv}[1]{{\mbox{$\mathbf{#1}$}}}
       \mathscr{G}(\bv{Y},\bv{\alpha},\bv{\beta}) =
       \prod_{i=1}^K\frac{\Gamma(\alpha_i+\beta_i)}{\Gamma(\alpha_i)
       \Gamma(\beta_i)}Y_i^{\alpha_i-1} \mathcal{Y}_i^{\gamma_i} \qquad
       \mathrm{with} \qquad \mathcal{Y}_i = 1-\sum_{k=1}^i Y_k,
    \f]
    provided the coefficients, \f$b_i\!>\!0\f$, \f$\kappa_i\!>\!0\f$,
    \f$0\!<\!S_i\!<\!1\f$, and \f$c_{ij}\f$, with \f$c_{ij}\!=\!0\f$ for
    \f$i\!>\!j\f$, \f$i,j\!=\!1,\dots,K\!-\!1\f$, satisfy
    \f[ \begin{split}
       \alpha_i & = \frac{b_i}{\kappa_i}S_i, \qquad i=1,\dots,K,\\
       1-\gamma_i & = \frac{c_{1i}}{\kappa_1} = \dots = \frac{c_{ii}}{\kappa_i},
       \qquad i=1,\dots,K-1,\\
       1+\gamma_K & = \frac{b_1}{\kappa_1}(1-S_1) = \dots =
       \frac{b_K}{\kappa_K}(1-S_K). \end{split}
    \f]

    Here \f$\mathcal{U}_i = \prod_{j=1}^{K-i}\mathcal{Y}_{K-j}^{-1}\f$,
    \f$\alpha_i>0\f$, and \f$\beta_i>0\f$ are parameters, while
    \f$\gamma_i=\beta_i-\alpha_{i+1}-\beta_{i+1}\f$ for \f$i=1,\dots,K-1\f$, and
    \f$\gamma_K=\beta_K-1\f$. \f$\Gamma(\cdot)\f$ denotes the [gamma
    function](http://en.wikipedia.org/wiki/Gamma_function). To keep the
    invariant distribution generalized Dirichlet, the above set of constraints
    on the coefficients must be satisfied. For more details on the generalized
    Dirichlet SDE, see http://dx.doi.org/10.1063/1.4822416.
*/
// *****************************************************************************
#ifndef GeneralizedDirichlet_h
#define GeneralizedDirichlet_h

#include <vector>

#include "InitPolicy.h"
#include "GeneralizedDirichletCoeffPolicy.h"
#include "RNG.h"
#include "Particles.h"
#include "SystemComponents.h"

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief Lochner's generalized Dirichlet SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see
//!       DiffEq/GeneralizedDirichletCoeffPolicy.h
template< class Init, class Coefficients >
class GeneralizedDirichlet {

  private:
    using ncomp_t = tk::ctr::ncomp_type;

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which generalized Dirichlet system of SDEs
    //!   to construct. There can be multiple gendir ... end blocks in a control
    //!   file. This index specifies which generalized Dirichlet SDE system to
    //!   instantiate. The index corresponds to the order in which the gendir
    //!   ... end blocks are given the control file.
    //! \author J. Bakosi
    explicit GeneralizedDirichlet( ncomp_t c ) :
      m_c( c ),
      m_depvar( g_inputdeck.get< tag::param, tag::gendir, tag::depvar >().at(c) ),
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::gendir >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::gendir >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::gendir, tag::rng >().at(c) ) ) ),
      m_b(),
      m_S(),
      m_k(),
      m_cij(),
      coeff( m_ncomp,
             g_inputdeck.get< tag::param, tag::gendir, tag::b >().at(c),
             g_inputdeck.get< tag::param, tag::gendir, tag::S >().at(c),
             g_inputdeck.get< tag::param, tag::gendir, tag::kappa >().at(c),
             g_inputdeck.get< tag::param, tag::gendir, tag::c >().at(c),
             m_b, m_S, m_k, m_cij ) {}

    //! Initalize SDE, prepare for time integration
    //! \param[in] stream Thread (or more precisely stream) ID 
    //! \param[in,out] particles Array of particle properties 
    //! \author J. Bakosi
    void initialize( int stream, tk::Particles& particles ) {
      //! Set initial conditions using initialization policy
      Init::template
        init< tag::gendir >
            ( g_inputdeck, m_rng, stream, particles, m_c, m_ncomp, m_offset );
    }

    //! \brief Advance particles according to the generalized Dirichlet SDE
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
      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Y_i = 1 - sum_{k=1}^{i} y_k
        std::vector< tk::real > Y( m_ncomp );
        Y[0] = 1.0 - particles( p, 0, m_offset );
        for (ncomp_t i=1; i<m_ncomp; ++i)
          Y[i] = Y[i-1] - particles( p, i, m_offset );

        // U_i = prod_{j=1}^{K-i} 1/Y_{K-j}
        std::vector< tk::real > U( m_ncomp );
        U[m_ncomp-1] = 1.0;
        for (long i=static_cast<long>(m_ncomp)-2; i>=0; --i) {
          auto I = static_cast< std::size_t >( i );
          U[I] = U[I+1]/Y[I];
        }

        // Generate Gaussian random numbers with zero mean and unit variance
        std::vector< tk::real > dW( m_ncomp );
        m_rng.gaussian( stream, m_ncomp, dW.data() );

        // Advance first m_ncomp (K=N-1) scalars
        ncomp_t k=0;
        for (ncomp_t i=0; i<m_ncomp; ++i) {
          tk::real& par = particles( p, i, m_offset );
          tk::real d = m_k[i] * par * Y[m_ncomp-1] * U[i] * dt;
          d = (d > 0.0 ? sqrt(d) : 0.0);
          tk::real a=0.0;
          for (ncomp_t j=i; j<m_ncomp-1; ++j) a += m_cij[k++]/Y[j];
          par += U[i]/2.0*( m_b[i]*( m_S[i]*Y[m_ncomp-1] - (1.0-m_S[i])*par ) +
                            par*Y[m_ncomp-1]*a )*dt + d*dW[i];
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
    std::vector< kw::sde_c::info::expect::type > m_cij;

    //! Coefficients policy
    Coefficients coeff;
};

} // walker::

#endif // GeneralizedDirichlet_h
