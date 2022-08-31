// *****************************************************************************
/*!
  \file      src/PDE/EoS/JWL.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Stiffened-gas equation of state
  \details   This file defines functions for the stiffened gas equation of
             state for the compressible flow equations.
*/
// *****************************************************************************
#ifndef JWL_h
#define JWL_h

#include <cmath>
#include <iostream>
#include "Data.hpp"
#include "EoS_Base.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

class JWL: public EoS_Base {

  private:
    tk::real m_w, m_cv, m_rho0, m_e0, m_a, m_b, m_r1, m_r2, m_tr, m_rhor;


    tk::real intEnergy( tk::real rho, tk::real pr )
    {
      tk::real rho0 = m_rho0;
      tk::real a = m_a;
      tk::real b = m_b;
      tk::real r1 = m_r1;
      tk::real r2 = m_r2;
      tk::real w_jwl = m_w;
      tk::real e0 = m_e0;

      tk::real e = e0 + 1.0/w_jwl/rho*( pr
                    - a*(1.0 - w_jwl*rho/r1/rho0)*exp(-r1*rho0/rho)
                    - b*(1.0 - w_jwl*rho/r2/rho0)*exp(-r2*rho0/rho) );

      return e;
    }


    tk::real bisection( tk::real a, tk::real b, tk::real p_known, tk::real t_known )
    {
      tk::real tol = 1e-10;
      std::size_t maxiter = 1000;
      std::size_t i(0);
      tk::real c;
      tk::real root(0);

      // function to minimize = p_known - PfromRT
      // bounds b > a

      while (i < maxiter)
      {
        c = (a + b)/2.0;
        if ( p_known - PfromRT( c, t_known) <= 1e-16 or (b - a)/2.0 < tol )
        {
          root = c;
          break;
        }

        i++;
        if ( static_cast< int > (std::copysign( 1.0, p_known - PfromRT( c, t_known) )) ==
             static_cast< int > (std::copysign( 1.0, p_known - PfromRT( a, t_known) )) )
        {
          a = c;
        }
        else
        {
          b = c;
        }

      }
      return root;
    }


    tk::real PfromRT( tk::real rho, tk::real T)
    {
      tk::real rho0 = m_rho0;
      tk::real a = m_a;
      tk::real b = m_b;
      tk::real r1 = m_r1;
      tk::real r2 = m_r2;
      tk::real w_jwl = m_w;
      tk::real c_v = m_cv;
//      tk::real t_r = m_tr;      // reference temperature
//      tk::real rho_r = m_rhor;  // reference density
      tk::real t_r = 1200.0;      // reference temperature
      tk::real rho_r = 5.0e3;  // reference density

      tk::real pr;

      pr = a*exp(-r1*rho0/rho) + b*exp(-r2*rho0/rho) + w_jwl*(c_v*T*rho
         - c_v*t_r*std::pow(rho0/rho_r, w_jwl)*rho0/std::pow(rho0/rho, w_jwl+1));

      return pr;
    }

  public:
    // *************************************************************************
    //  Constructor
    //! \param[in] w Grueneisen coefficient
    //! \param[in] cv Specific heat at constant volume
    //! \param[in] rho0 Density of reference state
    //! \param[in] e0 Internal energy of reference state
    //! \param[in] A Parameter A
    //! \param[in] B Parameter B
    //! \param[in] R1 Parameter R1
    //! \param[in] R2 Parameter R2
    // *************************************************************************
    JWL( tk::real w, tk::real cv, tk::real rho0, tk::real e0, tk::real A,
         tk::real B, tk::real R1, tk::real R2 ) :
      m_w(w),
      m_cv(cv),
      m_rho0(rho0),
      m_e0(e0),
      m_a(A),
      m_b(B),
      m_r1(R1),
      m_r2(R2)
    { }


    tk::real eos_density( tk::real pr,
                          tk::real temp ) override
    // *************************************************************************
    //! \brief Calculate density from the material pressure and temperature 
    //!   using the stiffened-gas equation of state
    //! \param[in] pr Material pressure
    //! \param[in] temp Material temperature
    //! \return Material density calculated using the stiffened-gas EoS
    // *************************************************************************
    {
      tk::real rho_r = m_rhor;  // reference density
      tk::real r_guessL = 1e-2*rho_r;  // left density bound
      tk::real r_guessR = 1e2*rho_r;   // right density bound
      tk::real rho;

      rho = bisection( r_guessL, r_guessR, pr, temp );
    
      return rho;
    }


    tk::real eos_pressure( tk::real arho,
                           tk::real u,
                           tk::real v,
                           tk::real w,
                           tk::real arhoE,
                           tk::real alpha=1.0,
                           std::size_t imat=0 ) override
    // *************************************************************************
    //! \brief Calculate pressure from the material density, momentum and total
    //!   energy using the stiffened-gas equation of state
    //! \param[in] arho Material partial density (alpha_k * rho_k)
    //! \param[in] u X-velocity
    //! \param[in] v Y-velocity
    //! \param[in] w Z-velocity
    //! \param[in] arhoE Material total energy (alpha_k * rho_k * E_k)
    //! \param[in] alpha Material volume fraction. Default is 1.0, so that for
    //!   the single-material system, this argument can be left unspecified by
    //!   the calling code
    //! \param[in] imat Material-id who's EoS is required. Default is 0, so that
    //!   for the single-material system, this argument can be left unspecified
    //!   by the calling code
    //! \return Material partial pressure (alpha_k * p_k) calculated using the
    //!   stiffened-gas EoS
    // *************************************************************************
    {
      tk::real rho0 = m_rho0;
      tk::real a = m_a;
      tk::real b = m_b;
      tk::real r1 = m_r1;
      tk::real r2 = m_r2;
      tk::real w_jwl = m_w;
      tk::real e0 = m_e0;

      // reference energy (input quantity, might need for calculation)
//      tk::real e0 = a/r1*exp(-r1*rho0/rho) + b/r2*exp(-r2*rho0/rho);
      // internal energy
      tk::real rhoe = (arhoE - 0.5*arho*(u*u + v*v + w*w))/alpha;

      tk::real partpressure = a*(alpha*1.0 - w_jwl*arho/(rho0*r1))*exp(-r1*alpha*rho0/arho)
                            + b*(alpha*1.0 - w_jwl*arho/(rho0*r2))*exp(-r2*alpha*rho0/arho)
                            + w_jwl*(rhoe - e0)*arho/rho0;

      // check partial pressure divergence
      if (!std::isfinite(partpressure)) {
        std::cout << "Material-id:      " << imat << std::endl;
        std::cout << "Volume-fraction:  " << alpha << std::endl;
        std::cout << "Partial density:  " << arho << std::endl;
        std::cout << "Total energy:     " << arhoE << std::endl;
        std::cout << "Velocity:         " << u << ", " << v << ", " << w
          << std::endl;
        Throw("Material-" + std::to_string(imat) +
          " has nan/inf partial pressure: " + std::to_string(partpressure) +
          ", material volume fraction: " + std::to_string(alpha));
      }

      return partpressure;
    }


    tk::real eos_soundspeed( tk::real arho,
                             tk::real apr,
                             tk::real alpha=1.0,
                             std::size_t imat=0 ) override
    // *************************************************************************
    //! Calculate speed of sound from the material density and material pressure
    //! \param[in] arho Material partial density (alpha_k * rho_k)
    //! \param[in] apr Material partial pressure (alpha_k * p_k)
    //! \param[in] alpha Material volume fraction. Default is 1.0, so that for
    //!   the single-material system, this argument can be left unspecified by
    //!   the calling code
    //! \param[in] imat Material-id who's EoS is required. Default is 0, so that
    //!   for the single-material system, this argument can be left unspecified
    //!   by the calling code
    //! \return Material speed of sound using the stiffened-gas EoS
    // *************************************************************************
    {
      tk::real rho0 = m_rho0;
      tk::real a = m_a;
      tk::real b = m_b;
      tk::real r1 = m_r1;
      tk::real r2 = m_r2;
      tk::real w_jwl = m_w;

      // limiting pressure to near-zero
      auto apr_eff = std::max( 1.0e-15, apr );

      auto co1 = rho0*alpha*alpha/(arho*arho);
      auto co2 = alpha*(1.0+w_jwl)/arho;

      tk::real ss = a*(r1*co1 - co2) * exp(-r1*alpha*rho0/arho)
                  + b*(r2*co1 - co2) * exp(-r2*alpha*rho0/arho)
                  + (1.0+w_jwl)*apr_eff/arho;

      ss = std::sqrt(ss);
    
      // check sound speed divergence
      if (!std::isfinite(ss)) {
        std::cout << "Material-id:      " << imat << std::endl;
        std::cout << "Volume-fraction:  " << alpha << std::endl;
        std::cout << "Partial density:  " << arho << std::endl;
        std::cout << "Partial pressure: " << apr << std::endl;
        Throw("Material-" + std::to_string(imat) + " has nan/inf sound speed: "
          + std::to_string(ss) + ", material volume fraction: " +
          std::to_string(alpha));
      }
    
      return ss;
    }


    tk::real eos_totalenergy( tk::real rho,
                              tk::real u,
                              tk::real v,
                              tk::real w,
                              tk::real pr ) override
    // *************************************************************************
    //! \brief Calculate material specific total energy from the material
    //!   density, momentum and material pressure
    //! \param[in] rho Material density
    //! \param[in] u X-velocity
    //! \param[in] v Y-velocity
    //! \param[in] w Z-velocity
    //! \param[in] pr Material pressure
    //! \return Material specific total energy using the stiffened-gas EoS
    // *************************************************************************
    {

      // reference energy (input quantity, might need for calculation)
//      tk::real e0 = a/r1*exp(-r1*rho0/rho) + b/r2*exp(-r2*rho0/rho);
    
      tk::real rhoE = intEnergy( rho, pr )
                    + 0.5*rho*(u*u + v*v + w*w);

      return rhoE;
    }


    tk::real eos_temperature( tk::real arho,
                              tk::real u,
                              tk::real v,
                              tk::real w,
                              tk::real arhoE,
                              tk::real alpha=1.0 ) override
    // *************************************************************************
    //! \brief Calculate material temperature from the material density, and
    //!   material specific total energy
    //! \param[in] arho Material partial density (alpha_k * rho_k)
    //! \param[in] u X-velocity
    //! \param[in] v Y-velocity
    //! \param[in] w Z-velocity
    //! \param[in] arhoE Material total energy (alpha_k * rho_k * E_k)
    //! \param[in] alpha Material volume fraction. Default is 1.0, so that for
    //!   the single-material system, this argument can be left unspecified by
    //!   the calling code
    //! \return Material temperature using the stiffened-gas EoS
    // *************************************************************************
    {
      tk::real rho0 = m_rho0;
      tk::real a = m_a;
      tk::real b = m_b;
      tk::real r1 = m_r1;
      tk::real r2 = m_r2;
      tk::real w_jwl = m_w;
      tk::real c_v = m_cv;      // constant specific heat
//      tk::real t_r = m_tr;      // reference temperature
//      tk::real rho_r = m_rhor;  // reference density
      tk::real t_r = 1200.0;      // reference temperature
      tk::real rho_r = 5.0e3;  // reference density

      tk::real rho = arho/alpha;

      // reference energy (input quantity, might need for calculation)
//      tk::real e0 = a/r1*exp(-r1*rho0/rho) + b/r2*exp(-r2*rho0/rho);
    
      tk::real t = ((arhoE - 0.5*arho*(u*u + v*v + w*w))/arho - 1.0/rho0*(
                   a/r1*exp(-r1*rho0/rho)
                 + b/r2*exp(-r2*rho0/rho) ))/c_v
                 + ( t_r*std::pow(rho0/rho_r, w_jwl) )/std::pow(rho0/rho, w_jwl);

      return t;
    }


    // Destructor
    ~JWL() override {}
};

} //inciter::

#endif // JWL_h

