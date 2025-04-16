// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Mixture.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Multispecies mixture function
  \details   This file declares functions for computing mixture flow quantities
*/
// *****************************************************************************
#ifndef Mixture_h
#define Mixture_h

namespace inciter {

class Mixture {

  private:
    std::size_t m_nspec;
    std::vector< EOS > m_mat_blk;
    tk::real m_mix_density;
    tk::real m_mix_R;
    std::vector< tk::real > m_Ys;

  public:
    //! Constructor
    Mixture(const std::size_t nspec,
            const std::vector< EOS > mat_blk);

    //! Set mixture properties based off given state vector
    void set_state(std::vector< tk::real > ugp) const;

    //! Set mixture properties based off given mass frac, pres, temp
    void set_massfrac(std::vector< tk::real > Ys,
                      tk::real mix_pressure,
                      tk::real temperature) const;

    //! Return mixture density
    tk::real get_mix_density() { return m_mix_density; }

    //! Compute mixture frozen speed of sound.
    tk::real frozen_soundspeed(tk::real mix_density,
                               tk::real mix_pressure) const;

    //! Compute mixture total energy
    tk::real totalenergy(tk::real mix_density,
                         tk::real u,
                         tk::real v,
                         tk::real w,
                         tk::real mix_pressure) const;

    //! Compute mixture pressure
    tk::real pressure(tk::real mix_density,
                      tk::real u,
                      tk::real v,
                      tk::real w,
                      tk::real rhoE) const;

    //! Compute mixture temperature
    tk::real temperature(tk::real mix_density,
                         tk::real u,
                         tk::real v,
                         tk::real w,
                         tk::real rhoE) const;

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) /*override*/ {
      p | m_nspec;
      p | m_mat_blk;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i Mixture object reference
    friend void operator|( PUP::er& p, Mixture& i ) { i.pup(p); }
    //@}

};

} //inciter::

#endif // Mixture_h
