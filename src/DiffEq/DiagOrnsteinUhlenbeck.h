//******************************************************************************
/*!
  \file      src/DiffEq/DiagOrnsteinUhlenbeck.h
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 11:42:26 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     System of diagonal Ornstein-Uhlenbeck SDEs
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs), with linear drift and constant diagonal
    diffusion, whose invariant is the joint normal distribution.
*/
//******************************************************************************
#ifndef DiagOrnsteinUhlenbeck_h
#define DiagOrnsteinUhlenbeck_h

#include <cmath>

#include <InitPolicy.h>
#include <DiagOrnsteinUhlenbeckCoeffPolicy.h>
#include <RNG.h>

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief Diagonal Ornstein-Uhlenbeck SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see
//!       DiffEq/DiagOrnsteinUhlenbeckCoeffPolicy.h
template< class Init, class Coefficients >
class DiagOrnsteinUhlenbeck {

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of diagonal
    //!   Ornstein-Uhlenbeck SDEs to construct. There can be multiple diag_ou
    //!   ... end blocks in a control file. This index specifies which diagonal
    //!   Ornstein-Uhlenbeck SDE system to instantiate. The index corresponds to
    //!   the order in which the diag_ou ... end blocks are given the control
    //!   file.
    //! \author J. Bakosi
    explicit DiagOrnsteinUhlenbeck( unsigned int c ) :
      m_depvar( g_inputdeck.get< tag::param, tag::diagou, tag::depvar >().at(c) ),
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::diagou >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::diagou >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::diagou, tag::rng >().at(c) ) ) ),
      coeff( m_ncomp,
             g_inputdeck.get< tag::param, tag::diagou, tag::sigma >().at(c),
             g_inputdeck.get< tag::param, tag::diagou, tag::theta >().at(c),
             g_inputdeck.get< tag::param, tag::diagou, tag::mu >().at(c),
             m_sigma, m_theta, m_mu ) {}

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

    //! \brief Advance particles according to the system of diagonal
    //!   Orsntein-Uhlenbeck SDEs
    //! \author J. Bakosi
    void advance( tk::ParProps& particles, int stream, tk::real dt ) const {
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        tk::real dW[m_ncomp];
        m_rng.gaussian( stream, m_ncomp, dW );

        // Advance all m_ncomp scalars
        for (tk::ctr::ncomp_type i=0; i<m_ncomp; ++i) {
          tk::real& par = particles( p, i, m_offset );
          tk::real d = m_sigma[i] * m_sigma[i] * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          par += m_theta[i]*(m_mu[i] - par)*dt + d*dW[i];
        }
      }
    }

  private:
    const char m_depvar;                //!< Dependent variable
    const tk::ctr::ncomp_type m_ncomp;  //!< Number of components
    const int m_offset;                 //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator

    //! Coefficients
    std::vector< kw::sde_sigma::info::expect::type > m_sigma;
    std::vector< kw::sde_theta::info::expect::type > m_theta;
    std::vector< kw::sde_mu::info::expect::type > m_mu;

    //! Coefficients policy
    Coefficients coeff;
};

} // walker::

#endif // DiagOrnsteinUhlenbeck_h
