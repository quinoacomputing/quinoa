//******************************************************************************
/*!
  \file      src/Physics/Physics.h
  \author    J. Bakosi
  \date      Mon Sep  9 08:21:28 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************
#ifndef Physics_h
#define Physics_h

#include <QuinoaConfig.h>
#include <Mass.h>
#include <Mix.h>
#include <Hydro.h>

namespace quinoa {

class Memory;
class Paradigm;
class QuinoaControl;
class Timer;
class Statistics;
class GlobWriter;
class TxtStatWriter;
class Dirichlet;
class GeneralizedDirichlet;
class SimplifiedLangevin;
class GeneralizedLangevin;
class Beta;

//! Physics base
class Physics {

  private:
  // Select models based on <build>/Base/QuinoaConfig.h filled by CMake based
  // on src/MainQuinoaConfig.h.in

  // Mass models
  #ifdef QUINOA_BETA
    using MassType = Beta;
  #else
    #error "No mass model defined in Base/QuinaConfig.h"
  #endif

  // Hydrodynamics models
  #ifdef QUINOA_SLM
    using HydroType = SimplifiedLangevin;
  #elif QUINOA_GLM
    using HydroType = GeneralizedLangevin;
  #else
    #error "No hydrodynamics model defined in Base/QuinaConfig.h"
  #endif

  // Mix models
  #ifdef QUINOA_DIRICHLET
    using MixType = Dirichlet;
  #elif QUINOA_GENERALIZED_DIRICHLET
    using MixType = GeneralizedDirichlet;
  #else
    #error "No mix model defined in Base/QuinaConfig.h"
  #endif

  public:
    //! Destructor
    virtual ~Physics() noexcept;

    //! Initialize physics
    virtual void init() = 0;

    //! Solve physics
    virtual void solve() = 0;

    //! Constant accessor to control object
    //! \return Control object
    const QuinoaControl& control() const noexcept { return m_control; }

    //! Constant accessor to timer object pointer
    //! \return Pointer to timer object
    Timer* timer() const noexcept { return m_timer; }

    //! Constant accessor to mass model
    //! \return Pointer to mass model
    Mass<MassType>* mass() const noexcept { return m_mass; }

    //! Constant accessor to hydro model
    //! \return Pointer to hydro model
    Hydro<HydroType>* hydro() const noexcept { return m_hydro; }

    //! Constant accessor to mix model
    //! \return Pointer to mix model
    Mix<MixType>* mix() const noexcept { return m_mix; }

    //! Constant accessor to statistics estimator
    //! \return Pointer to statistics estimator
    Statistics* statistics() const noexcept { return m_statistics; }

    //! Constant accessor to glob file writer
    //! \return Pointer to glob file writer
    GlobWriter* globWriter() const noexcept { return m_glob; }

    //! Constant accessor to statistics file writer
    //! \return Pointer to statistics file writer
    TxtStatWriter* statWriter() const noexcept { return m_stat; }

    //! Constant accessor to particle properties pointer
    //! \return Raw pointer to particle properties array
    real* particles() const noexcept { return m_particles.ptr; }

  protected:
    //! Constructor: protected, designed to be base-only
    explicit Physics(Memory* const memory,
                     Paradigm* const paradigm,
                     const QuinoaControl& control,
                     Timer* const timer);

    //! Echo information on physics
    void echo();

    const int m_nposition;                //!< Number of position components
    const int m_ndensity;                 //!< Number of density components
    const int m_nvelocity;                //!< Number of velocity components
    const int m_nscalar;                  //!< Number of scalar components
    const uint64_t m_npar;                //!< Number of particles
    const real m_term;                    //!< Maximum time to simulate

  private:
    //! Don't permit copy constructor
    Physics(const Physics&) = delete;
    //! Don't permit copy assigment
    Physics& operator=(const Physics&) = delete;
    //! Don't permit move constructor
    Physics(Physics&&) = delete;
    //! Don't permit move assigment
    Physics& operator=(Physics&&) = delete;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

    Memory* const m_memory;               //!< Memory object
    Paradigm* const m_paradigm;           //!< Parallel programming object
    const QuinoaControl& m_control;       //!< Control object
    Timer* const m_timer;                 //!< Timer object

    Mass<MassType>* m_mass;               //!< Mass model object    
    Hydro<HydroType>* m_hydro;            //!< Hydro model object    
    Mix<MixType>* m_mix;                  //!< Mix model object
    Statistics* m_statistics;             //!< Statistics estimator object
    GlobWriter* m_glob;                   //!< Glob file writer
    TxtStatWriter* m_stat;                //!< Statistics file writer
    Data<real> m_particles;               //!< Particle properties
};

} // namespace quinoa

#endif // Physics_h
