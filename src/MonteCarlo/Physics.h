//******************************************************************************
/*!
  \file      src/MonteCarlo/Physics.h
  \author    J. Bakosi
  \date      Wed Mar 19 16:10:46 2014
  \copyright 2005-2014, Jozsef Bakosi.
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

namespace quinoa {

//! Physics : MonteCarlo
class Physics : public MonteCarlo {

  public:
    //! Destructor
    ~Physics() override = default;

    //! Accessor to mass model
    //! \return Pointer to mass model
    Model* mass() const noexcept { return m_mass.get(); }

    //! Accessor to hydro model
    //! \return Pointer to hydro model
    Model* hydro() const noexcept { return m_hydro.get(); }

    //! Accessor to mix model
    //! \return Pointer to mix model
    Model* mix() const noexcept { return m_mix.get(); }

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
    void initFactories();

    //! Echo information on physics
    void echo();

    //! Factories
    ctr::MassFactory m_massFactory;             //!< Mass model factory
    ctr::HydroFactory m_hydroFactory;           //!< Hydrodynamics model factory
    ctr::MixFactory m_mixFactory;               //!< Material mix model factory

    //! Pointers to selected options
    std::unique_ptr< Model > m_mass;            //!< Mass model
    std::unique_ptr< Model > m_hydro;           //!< Hydro model
    std::unique_ptr< Model > m_mix;             //!< Mix model
};

} // quinoa::

#endif // Physics_h
