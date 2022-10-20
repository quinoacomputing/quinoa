// *****************************************************************************
/*!
  \file      src/PDE/EoS/SGclass.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Stiffened-gas equation of state
  \details   This file defines functions for the stiffened gas equation of
             state for the compressible flow equations.
*/
// *****************************************************************************
#ifndef SGclass_h
#define SGclass_h

#include <cmath>
#include <iostream>
#include "Data.hpp"

namespace inciter {

class SGclass {

  private:
    tk::real m_gamma, m_pstiff, m_cv;

  public:
    // *************************************************************************
    //  Default constructor
    // *************************************************************************
    SGclass(){}


    // *************************************************************************
    //  Constructor
    //! \param[in] gamma Ratio of specific heats
    //! \param[in] pstiff Stiffened pressure term
    //! \param[in] cv Specific heat at constant volume
    // *************************************************************************
    SGclass(tk::real gamma, tk::real pstiff, tk::real cv ) :
      m_gamma(gamma), m_pstiff(pstiff), m_cv(cv)
    { }


    tk::real eos_density( tk::real pr,
                          tk::real temp ) const
    // *************************************************************************
    //! \brief Calculate density from the material pressure and temperature 
    //!   using the stiffened-gas equation of state
    //! \param[in] pr Material pressure
    //! \param[in] temp Material temperature
    //! \return Material density calculated using the stiffened-gas EoS
    // *************************************************************************
    {
      tk::real g = m_gamma;
      tk::real p_c = m_pstiff;
      tk::real c_v = m_cv;
    
      tk::real rho = (pr + p_c) / ((g-1.0) * c_v * temp);
      return rho;
    }


    tk::real eos_pressure( tk::real arho,
                           tk::real u,
                           tk::real v,
                           tk::real w,
                           tk::real arhoE,
                           tk::real alpha=1.0,
                           std::size_t imat=0 ) const
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
      tk::real g = m_gamma;
      tk::real p_c = m_pstiff;

      tk::real partpressure = (arhoE - 0.5 * arho * (u*u + v*v + w*w) -
        alpha*p_c) * (g-1.0) - alpha*p_c;

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
                             std::size_t imat=0 ) const
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
      auto g = m_gamma;
      auto p_c = m_pstiff;
    
      auto p_eff = std::max( 1.0e-15, apr+(alpha*p_c) );
    
      tk::real a = std::sqrt( g * p_eff / arho );
    
      // check sound speed divergence
      if (!std::isfinite(a)) {
        std::cout << "Material-id:      " << imat << std::endl;
        std::cout << "Volume-fraction:  " << alpha << std::endl;
        std::cout << "Partial density:  " << arho << std::endl;
        std::cout << "Partial pressure: " << apr << std::endl;
        Throw("Material-" + std::to_string(imat) + " has nan/inf sound speed: "
          + std::to_string(a) + ", material volume fraction: " +
          std::to_string(alpha));
      }
    
      return a;
    }


    tk::real eos_totalenergy( tk::real rho,
                              tk::real u,
                              tk::real v,
                              tk::real w,
                              tk::real pr ) const
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
      auto g = m_gamma;
      auto p_c = m_pstiff;
    
      tk::real rhoE = (pr + p_c) / (g-1.0) + 0.5 * rho * (u*u + v*v + w*w) + p_c;
      return rhoE;
    }


    tk::real eos_temperature( tk::real arho,
                              tk::real u,
                              tk::real v,
                              tk::real w,
                              tk::real arhoE,
                              tk::real alpha=1.0 ) const
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
      auto c_v = m_cv;
      auto p_c = m_pstiff;
    
      tk::real t = (arhoE - 0.5 * arho * (u*u + v*v + w*w) - alpha*p_c) 
                   / (arho*c_v);
      return t;
    }


    tk::real min_eff_pressure( tk::real min ) const
    // *************************************************************************
    //! Compute the minimum effective pressure
    //! \param[in] min Minimum threshold in positivity preserving limiting
    //! \return Minimum effective pressure
    // *************************************************************************
    {
      auto p_c = m_pstiff;

      return (min - p_c);
    }


    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) /*override*/ {
      p | m_gamma;
      p | m_pstiff;
      p | m_cv;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i SGclass object reference
    friend void operator|( PUP::er& p, SGclass& i ) { i.pup(p); }
    //@}
};

} //inciter::

#endif // SGclass_h
