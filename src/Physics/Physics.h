//******************************************************************************
/*!
  \file      src/Physics/Physics.h
  \author    J. Bakosi
  \date      Fri May 10 18:00:28 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************
#ifndef Physics_h
#define Physics_h

#include <QuinoaConfig.h>
#include <Mix.h>
#include <Hydro.h>

using namespace std;

namespace Quinoa {

class Memory;
class Paradigm;
class Control;
class Timer;
class Statistics;
class GlobWriter;
class TxtPlotWriter;
class Dirichlet;
class GeneralizedDirichlet;
class SimplifiedLangevin;
class GeneralizedLangevin;

//! Physics base
class Physics {

  public:
    //! Constructor
    explicit Physics(Memory* const memory,
                     Paradigm* const paradigm,
                     Control* const control,
                     Timer* const timer);

    //! Destructor
    virtual ~Physics() noexcept;

    //! Echo informaion on physics
    virtual void echo() const = 0;

    //! Initialize physics
    virtual void init() = 0;

    //! Solve physics
    virtual void solve() = 0;

    //! Constant accessor to control object pointer
    //! \return Pointer to control object
    const Control* control() const noexcept { return m_control; }

    //! Constant accessor to timer object pointer
    //! \return Pointer to timer object
    Timer* timer() const noexcept { return m_timer; }

    //! Constant accessor to hydro model
    //! \return Pointer to hydro model
    //Hydro* hydro() const noexcept { return m_hydro; }

    //! Constant accessor to mix model
    //! \return Pointer to mix model
    //Mix* mix() const noexcept { return m_mix; }

    //! Constant accessor to statistics estimator
    //! \return Pointer to statistics estimator
    Statistics* statistics() const noexcept { return m_statistics; }

    //! Constant accessor to glob file writer
    //! \return Pointer to glob file writer
    GlobWriter* globWriter() const noexcept { return m_glob; }

    //! Constant accessor to plot file writer
    //! \return Pointer to plot file writer
    TxtPlotWriter* plotWriter() const noexcept { return m_plot; }

    //! Constant accessor to particle properties pointer
    //! \return Raw pointer to particle properties array
    const real* particles() const noexcept { return m_particles.ptr; }

    //! Accessor to number of particle (velocity+scalar) properties
    //! \return Number of particle (velocity+scalar) components
    int nprop() const noexcept { return m_nvelocity + m_nscalar; }

  protected:
    const int m_nthread;                  //!< Number of threads
    const int m_nvelocity;                //!< Number of velocity components
    const int m_nscalar;                  //!< Number of scalar components
    const int m_npar;                     //!< Numer of particles
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
    Control* const m_control;             //!< Control object
    Timer* const m_timer;                 //!< Timer object


  // Select models based on <build>/Base/QuinoaConfig.h filled by CMake based
  // on src/MainQuinoaConfig.h.in
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

    Mix<MixType>* m_mix;                  //!< Mix model object
    Hydro<HydroType>* m_hydro;            //!< Hydro model object    
    Statistics* m_statistics;             //!< Statistics estimator object
    GlobWriter* m_glob;                   //!< Glob file writer
    TxtPlotWriter* m_plot;                //!< Plot file writer
    Data<real> m_particles;               //!< Particle properties
};

} // namespace Quinoa

#endif // Physics_h
