//******************************************************************************
/*!
  \file      src/DiffEq/Dirichlet.h
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 11:30:21 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Dirichlet SDE
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs), whose invariant is the Dirichlet
    distribution. For more info on the Dirichlet SDE, see
    http://dx.doi.org/10.1155/2013/842981,
    \f[ \mathrm{d}Y_\alpha(t) = \frac{b_\alpha}{2}\big[S_\alpha Y_N -
        (1-S_\alpha)Y_\alpha\big]\mathrm{d}t + \sqrt{\kappa_\alpha Y_\alpha
        Y_N}\mathrm{d}W_\alpha(t), \qquad \alpha=1,\dots,N-1 \f]
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

#include <cmath>

#include <InitPolicy.h>
#include <DirichletCoeffPolicy.h>
#include <RNG.h>

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief Dirichlet SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see DiffEq/DirCoeffPolicy.h
template< class Init, class Coefficients >
class Dirichlet {

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which Dirichlet SDE to construct. There
    //!   can be multiple dirichlet ... end blocks in a control file. This index
    //!   specifies which Dirichlet SDE to instantiate. The index corresponds to
    //!   the order in which the dirichlet ... end blocks are given the control
    //!   file.
    //! \author J. Bakosi
    explicit Dirichlet( unsigned int c ) :
      m_depvar( g_inputdeck.get< tag::param, tag::dirichlet, tag::depvar >().at(c) ),
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::dirichlet >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::dirichlet >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::dirichlet, tag::rng >().at(c) ) ) ),
      coeff( m_ncomp,
             g_inputdeck.get< tag::param, tag::dirichlet, tag::b >().at(c),
             g_inputdeck.get< tag::param, tag::dirichlet, tag::S >().at(c),
             g_inputdeck.get< tag::param, tag::dirichlet, tag::kappa >().at(c),
             m_b, m_S, m_k ) {}

    //! Initalize SDE, prepare for time integration
    //! \param[inout] particles Array of particle properties 
    //! \param[in] stat Statistics object for accessing moments 
    //! \author J. Bakosi
    void initialize( tk::ParProps& particles, const tk::Statistics& stat ) {
      //! Set initial conditions using initialization policy
      Init( { particles } );
      //! Pre-lookup required statistical moments
      coeff.lookup( stat, m_depvar );
    }

    //! \brief Advance particles according to the Dirichlet SDE
    //! \author J. Bakosi
    void advance( tk::ParProps& particles, int stream, tk::real dt ) const {
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Compute Nth scalar
        tk::real yn = 1.0 - particles(p, 0, m_offset);
        for (tk::ctr::ncomp_type i=1; i<m_ncomp; ++i)
          yn -= particles( p, i, m_offset );

        // Generate Gaussian random numbers with zero mean and unit variance
        tk::real dW[m_ncomp];
        m_rng.gaussian( stream, m_ncomp, dW );

        // Advance first m_ncomp (K=N-1) scalars
        for (tk::ctr::ncomp_type i=0; i<m_ncomp; ++i) {
          tk::real& par = particles( p, i, m_offset );
          tk::real d = m_k[i] * par * yn * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          par += 0.5*m_b[i]*( m_S[i]*yn - (1.0-m_S[i]) * par )*dt + d*dW[i];
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

    //! Coefficients policy
    Coefficients coeff;
};

} // walker::

#endif // Dirichlet_h
