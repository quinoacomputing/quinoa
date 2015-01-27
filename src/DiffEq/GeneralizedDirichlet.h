//******************************************************************************
/*!
  \file      src/DiffEq/GeneralizedDirichlet.h
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 11:45:38 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Lochner's generalized Dirichlet SDE
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs) whose invariant is Lochner's generalized
    Dirichlet distribution. For more info on the generalized Dirichlet SDE, see
    http://dx.doi.org/10.1063/1.4822416.
*/
//******************************************************************************
#ifndef GeneralizedDirichlet_h
#define GeneralizedDirichlet_h

#include <InitPolicy.h>
#include <GeneralizedDirichletCoeffPolicy.h>
#include <RNG.h>

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

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which generalized Dirichlet system of SDEs
    //!   to construct. There can be multiple gendir ... end blocks in a control
    //!   file. This index specifies which generalized Dirichlet SDE system to
    //!   instantiate. The index corresponds to the order in which the gendir
    //!   ... end blocks are given the control file.
    //! \author J. Bakosi
    explicit GeneralizedDirichlet( unsigned int c ) :
      m_depvar( g_inputdeck.get< tag::param, tag::gendir, tag::depvar >().at(c) ),
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::gendir >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::gendir >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::gendir, tag::rng >().at(c) ) ) ),
      coeff( m_ncomp,
             g_inputdeck.get< tag::param, tag::gendir, tag::b >().at(c),
             g_inputdeck.get< tag::param, tag::gendir, tag::S >().at(c),
             g_inputdeck.get< tag::param, tag::gendir, tag::kappa >().at(c),
             g_inputdeck.get< tag::param, tag::gendir, tag::c >().at(c),
             m_b, m_S, m_k, m_c ) {}

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

    //! \brief Advance particles according to the generalized Dirichlet SDE
    //! \author J. Bakosi
    void advance( tk::ParProps& particles, int stream, tk::real dt ) const {
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Y_i = 1 - sum_{k=1}^{i} y_k
        tk::real Y[m_ncomp];
        Y[0] = 1.0 - particles( p, 0, m_offset );
        for (tk::ctr::ncomp_type i=1; i<m_ncomp; ++i)
          Y[i] = Y[i-1] - particles( p, i, m_offset );

        // U_i = prod_{j=1}^{K-i} 1/Y_{K-j}
        tk::real U[m_ncomp];
        U[m_ncomp-1] = 1.0;
        for (int i=m_ncomp-2; i>=0; --i) U[i] = U[i+1]/Y[i];

        // Generate Gaussian random numbers with zero mean and unit variance
        tk::real dW[m_ncomp];
        m_rng.gaussian( stream, m_ncomp, dW );

        // Advance first m_ncomp (K=N-1) scalars
        int k=0;
        for (tk::ctr::ncomp_type i=0; i<m_ncomp; ++i) {
          tk::real& par = particles( p, i, m_offset );
          tk::real d = m_k[i] * par * Y[m_ncomp-1] * U[i] * dt;
          d = (d > 0.0 ? sqrt(d) : 0.0);
          tk::real a=0.0;
          for (tk::ctr::ncomp_type j=i; j<m_ncomp-1; ++j) a += m_c[k++]/Y[j];
          par += U[i]/2.0*( m_b[i]*( m_S[i]*Y[m_ncomp-1] - (1.0-m_S[i])*par ) +
                            par*Y[m_ncomp-1]*a )*dt + d*dW[i];
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
    std::vector< kw::sde_c::info::expect::type > m_c;

    //! Coefficients policy
    Coefficients coeff;
};

} // walker::

#endif // GeneralizedDirichlet_h
