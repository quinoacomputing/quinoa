//******************************************************************************
/*!
  \file      src/MonteCarlo/Physics.h
  \author    J. Bakosi
  \date      Tue 31 Dec 2013 12:51:25 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************
#ifndef Physics_h
#define Physics_h

#include <Base.h>
#include <Statistics.h>
#include <GlobWriter.h>
#include <TxtStatWriter.h>
#include <MonteCarlo.h>
#include <Mass/Mass.h>
#include <Hydro/Hydro.h>
#include <Mix/Mix.h>

namespace quinoa {

//! Physics : MonteCarlo
class Physics : public MonteCarlo {

  public:
    //! Destructor
    ~Physics() override = default;

    //! Accessor to mass model
    //! \return Pointer to mass model
    Mass* mass() const noexcept { return m_mass.get(); }

    //! Accessor to hydro model
    //! \return Pointer to hydro model
    Hydro* hydro() const noexcept { return m_hydro.get(); }

    //! Accessor to mix model
    //! \return Pointer to mix model
    Mix* mix() const noexcept { return m_mix.get(); }

  protected:
    //! Constructor: protected, designed to be base-only
    explicit Physics( const Base& base );

    const int m_nposition;                //!< Number of position components
    const int m_ndensity;                 //!< Number of density components
    const int m_nvelocity;                //!< Number of velocity components
    const int m_nscalar;                  //!< Number of scalar components

  private:
    //! Don't permit copy constructor
    Physics(const Physics&) = delete;
    //! Don't permit copy assigment
    Physics& operator=(const Physics&) = delete;
    //! Don't permit move constructor
    Physics(Physics&&) = delete;
    //! Don't permit move assigment
    Physics& operator=(Physics&&) = delete;

    //! Initialize factories
    void initFactories(const tk::Print& print);

    //! Echo information on physics
    void echo();

    //! Factories
    ctr::MassFactory m_massFactory;             //!< Mass model factory
    ctr::HydroFactory m_hydroFactory;           //!< Hydrodynamics model factory
    ctr::MixFactory m_mixFactory;               //!< Material mix model factory

    //! Pointers to selected options
    std::unique_ptr< tk::RNG > m_rng;           //!< Random number generator
    std::unique_ptr< Mass > m_mass;             //!< Mass model
    std::unique_ptr< Hydro > m_hydro;           //!< Hydro model
    std::unique_ptr< Mix > m_mix;               //!< Mix model
};

} // quinoa::

#endif // Physics_h
