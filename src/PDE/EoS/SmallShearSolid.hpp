// *****************************************************************************
/*!
  \file      src/PDE/EoS/SmallShearSolid.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Small shear strain equation of state for solids
  \details   This file declares functions for the SmallShearSolid equation of
             state for the compressible flow equations. These functions are
             taken from Plohr, J. N., & Plohr, B. J. (2005). Linearized analysis
             of Richtmyer–Meshkov flow for elastic materials. Journal of Fluid
             Mechanics, 537, 55-89.
*/
// *****************************************************************************
#ifndef SmallShearSolid_h
#define SmallShearSolid_h

#include "Data.hpp"

namespace inciter {

class SmallShearSolid {

  private:
    tk::real m_gamma, m_pstiff, m_cv, m_mu;

  public:
    //! Default constructor
    SmallShearSolid() = default;

    //! Constructor
    SmallShearSolid(tk::real gamma, tk::real pstiff, tk::real cv, tk::real mu );


    //! Calculate density from the material pressure and temperature
    tk::real density( tk::real pr,
                      tk::real temp ) const;

    //! Calculate pressure from the material density, momentum and total energy
    tk::real pressure( tk::real arho,
                       tk::real u,
                       tk::real v,
                       tk::real w,
                       tk::real arhoE,
                       tk::real alpha=1.0,
                       std::size_t imat=0 ) const;

    //! Calculate speed of sound from the material density and material pressure
    tk::real soundspeed( tk::real arho,
                         tk::real apr,
                         tk::real alpha=1.0,
                         std::size_t imat=0 ) const;

    //! \brief Calculate material specific total energy from the material
    //!   density, momentum and material pressure
    tk::real totalenergy( tk::real rho,
                          tk::real u,
                          tk::real v,
                          tk::real w,
                          tk::real pr ) const;

    //! \brief Calculate material temperature from the material density, and
    //!   material specific total energy
    tk::real temperature( tk::real arho,
                          tk::real u,
                          tk::real v,
                          tk::real w,
                          tk::real arhoE,
                          tk::real alpha=1.0 ) const;

    //! Compute the minimum allowed pressure
    tk::real min_eff_pressure(
      tk::real min,
      tk::real,
      tk::real ) const;

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) /*override*/ {
      p | m_gamma;
      p | m_pstiff;
      p | m_cv;
      p | m_mu;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i SmallShearSolid object reference
    friend void operator|( PUP::er& p, SmallShearSolid& i ) { i.pup(p); }
    //@}
};

} //inciter::

#endif // SmallShearSolid_h
