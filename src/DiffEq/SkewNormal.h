//******************************************************************************
/*!
  \file      src/DiffEq/SkewNormal.h
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 11:54:54 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     System of skew-normal SDEs
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs), whose invariant is the joint skew-normal
    distribution. For more details on the skew-normal distribution, see
    http://www.jstor.org/stable/2337278.
*/
//******************************************************************************
#ifndef SkewNormal_h
#define SkewNormal_h

#include <cmath>

#include <InitPolicy.h>
#include <SkewNormalCoeffPolicy.h>
#include <RNG.h>

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief Skew-normal SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see DiffEq/SkewNormalCoeffPolicy.h
template< class Init, class Coefficients >
class SkewNormal {

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of skew-normal SDEs to
    //!   construct. There can be multiple skew-normal ... end blocks in a
    //!   control file. This index specifies which skew-normal SDE system to
    //!   instantiate. The index corresponds to the order in which the
    //!   skew-normal ... end blocks are given the control file.
    //! \author J. Bakosi
    explicit SkewNormal( unsigned int c ) :
      m_depvar( g_inputdeck.get< tag::param, tag::skewnormal, tag::depvar >().at(c) ),
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::skewnormal >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::skewnormal >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::skewnormal, tag::rng >().at(c) ) ) ),
      coeff( m_ncomp,
             g_inputdeck.get< tag::param, tag::skewnormal, tag::timescale >().at(c),
             g_inputdeck.get< tag::param, tag::skewnormal, tag::sigma >().at(c),
             g_inputdeck.get< tag::param, tag::skewnormal, tag::lambda >().at(c),
             m_T, m_sigma, m_lambda ) {}

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

    //! \brief Advance particles according to the system of skew-normal SDEs
    //! \author J. Bakosi
    void advance( tk::ParProps& particles, int stream, tk::real dt ) const {
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        tk::real dW[m_ncomp];
        m_rng.gaussian( stream, m_ncomp, dW );

        // Advance all m_ncomp scalars
        for (tk::ctr::ncomp_type i=0; i<m_ncomp; ++i) {
          tk::real& x = particles( p, i, m_offset );
          tk::real d = 2.0 * m_sigma[i] * m_sigma[i] / m_T[i] * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          x += - ( x - m_lambda[i] * m_sigma[i] * m_sigma[i]
                       * std::sqrt( 2.0 / pi() )
                       * std::exp( - m_lambda[i] * m_lambda[i] * x * x / 2.0 )
                       / ( 1.0 + std::erf( m_lambda[i] * x / std::sqrt(2.0) ) )
                   ) / m_T[i] * dt
                 + d*dW[i];
        }
      }
    }

  private:
    const char m_depvar;                //!< Dependent variable
    const tk::ctr::ncomp_type m_ncomp;    //!< Number of components
    const int m_offset;                   //!< Offset SDE operates from
    const tk::RNG& m_rng;                 //!< Random number generator

    //! Coefficients
    std::vector< kw::sde_T::info::expect::type > m_T;
    std::vector< kw::sde_sigma::info::expect::type > m_sigma;
    std::vector< kw::sde_lambda::info::expect::type > m_lambda;

    //! Coefficients policy
    Coefficients coeff;

    //! Return the value of Pi
    constexpr tk::real pi() const { return std::atan(1.0) * 4.0; }
};

} // walker::

#endif // SkewNormal_h
