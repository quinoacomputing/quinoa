//******************************************************************************
/*!
  \file      src/Physics/Physics.h
  \author    J. Bakosi
  \date      Thu Nov 14 08:18:30 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************
#ifndef Physics_h
#define Physics_h

#include <Base.h>
#include <RNG.h>
#include <Mass/Mass.h>
#include <Hydro/Hydro.h>
#include <Mix/Mix.h>
#include <Statistics.h>
#include <GlobWriter.h>
#include <TxtStatWriter.h>

namespace quinoa {

//! Random number generator factory type
using RNGFactory = std::map< quinoa::ctr::RNGType, std::function<tk::RNG*()> >;

//! Physics base
class Physics {

  public:
    //! Destructor
    virtual ~Physics() = default;

    //! Initialize physics
    virtual void init() = 0;

    //! Solve physics
    virtual void solve() = 0;

    //! Constant accessor to control object
    //! \return Control object
    const ctr::InputDeck& control() const noexcept { return m_base.control; }

    //! Constant accessor to timer object pointer
    //! \return Pointer to timer object
    const tk::Timer& timer() const noexcept { return m_base.timer; }

    //! Constant accessor to mass model
    //! \return Pointer to mass model
    Mass* mass() const noexcept { return m_mass.get(); }

    //! Constant accessor to hydro model
    //! \return Pointer to hydro model
    Hydro* hydro() const noexcept { return m_hydro.get(); }

    //! Constant accessor to mix model
    //! \return Pointer to mix model
    Mix* mix() const noexcept { return m_mix.get(); }

    //! Constant accessor to statistics estimator
    //! \return Pointer to statistics estimator
    Statistics& statistics() noexcept { return m_statistics; }

    //! Constant accessor to glob file writer
    //! \return Pointer to glob file writer
    GlobWriter& globWriter() noexcept { return m_glob; }

    //! Constant accessor to statistics file writer
    //! \return Pointer to statistics file writer
    TxtStatWriter& statWriter() noexcept { return m_stat; }

    //! Constant accessor to particle properties pointer
    //! \return Raw pointer to particle properties array
    tk::real* particles() const noexcept { return m_particles.get(); }

  protected:
    //! Constructor: protected, designed to be base-only
    explicit Physics(const Base& base);

    const int m_nposition;                //!< Number of position components
    const int m_ndensity;                 //!< Number of density components
    const int m_nvelocity;                //!< Number of velocity components
    const int m_nscalar;                  //!< Number of scalar components
    const uint64_t m_npar;                //!< Number of particles
    const tk::real m_term;                //!< Maximum time to simulate

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

    //! Register random number generators into factory
    void initRNGFactory( const quinoa::ctr::RNG& opt,
                         std::list< quinoa::ctr::RNGType >& reg,
                         int nthreads,
                         const quinoa::ctr::MKLRNGParameters& mklparam );

    //! Echo information on physics
    void echo();

    const Base& m_base;                             //!< Essentials
    const std::unique_ptr<tk::real[]> m_particles;  //!< Particle properties

    Statistics m_statistics;                        //!< Statistics estimator

    //! Factories
    RNGFactory m_RNGFactory;                    //!< RNG factory
    ctr::MassFactory m_massFactory;             //!< Mass model factory
    ctr::HydroFactory m_hydroFactory;           //!< Hydrodynamics model factory
    ctr::MixFactory m_mixFactory;               //!< Material mix model factory

    //! Pointers to selected options
    std::unique_ptr< tk::RNG > m_rng;           //!< Random number generator
    std::unique_ptr< Mass > m_mass;             //!< Mass model
    std::unique_ptr< Hydro > m_hydro;           //!< Hydro model
    std::unique_ptr< Mix > m_mix;               //!< Mix model

    GlobWriter m_glob;                          //!< Glob file writer
    TxtStatWriter m_stat;                       //!< Statistics file writer
};

} // quinoa::

#endif // Physics_h
