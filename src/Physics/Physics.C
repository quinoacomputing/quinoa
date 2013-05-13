//******************************************************************************
/*!
  \file      src/Physics/Physics.C
  \author    J. Bakosi
  \date      Sun 12 May 2013 09:39:11 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <QuinoaConfig.h>
#include <Memory.h>
#include <Paradigm.h>
#include <Control.h>
#include <Physics.h>
#include <Beta.h>
#include <Dirichlet.h>
#include <GeneralizedDirichlet.h>
#include <SimplifiedLangevin.h>
#include <GeneralizedLangevin.h>
#include <Statistics.h>
#include <GlobWriter.h>
#include <TxtPlotWriter.h>

using namespace std;
using namespace Quinoa;

Physics::Physics(Memory* const memory,
                 Paradigm* const paradigm,
                 Control* const control,
                 Timer* const timer)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object
//! \param[in]  paradigm Parallel programming object
//! \param[in]  control  Control object
//! \param[in]  timer    Timer object
//! \author  J. Bakosi
//******************************************************************************
try :
  m_nthread(paradigm->nthread()),
  m_ndensity(control->get<control::NDENSITY>()),
  m_nvelocity(control->get<control::NVELOCITY>()),
  m_nscalar(control->get<control::NSCALAR>()),
  m_npar(control->get<control::NPAR>()),
  m_term(control->get<control::TERM>()),
  m_memory(memory),
  m_paradigm(paradigm),
  m_control(control),
  m_timer(timer),
  m_hydro(nullptr),
  m_mix(nullptr),
  m_statistics(nullptr),
  m_glob(nullptr),
  m_plot(nullptr),
  m_particles()
{

  // Allocate memory to store all particle properties
  m_particles =
    m_memory->newEntry<real>((m_ndensity + m_nvelocity + m_nscalar) * m_npar,
                             REAL,
                             SCALAR,
                             "Particles");

  // Instantiate mass model
  if (m_ndensity) {
    m_mass = new (nothrow)
      MassType(memory, paradigm, control, m_particles + 0);
    ErrChk(m_mix != nullptr, FATAL, "Cannot allocate memory");
  }

  // Instantiate hydrodynamics model
  if (m_nvelocity) {
    m_hydro = new (nothrow)
      HydroType(memory,
                paradigm,
                control,
                m_particles + m_ndensity * m_npar);
    ErrChk(m_hydro != nullptr, FATAL, "Cannot allocate memory");
  }

  // Instantiate mix model
  if (m_nscalar) {
    m_mix = new (nothrow)
      MixType(memory,
              paradigm,
              control,
              m_particles + (m_ndensity + m_nvelocity) * m_npar);
    ErrChk(m_mix != nullptr, FATAL, "Cannot allocate memory");
  }

  // Instantiate statistics estimator
  m_statistics = new (nothrow) Statistics(memory, paradigm, control, this);
  ErrChk(m_statistics != nullptr, FATAL,"Cannot allocate memory");

  // Instantiate glob file writer
  m_glob = new (nothrow) GlobWriter(m_control->get<control::GLOBNAME>());
  ErrChk(m_glob != nullptr, FATAL,"Cannot allocate memory");

  // Instantiate plot file writer
  m_plot = new (nothrow) TxtPlotWriter(m_control->get<control::PLOTNAME>(),
                                       m_statistics);
  ErrChk(m_plot != nullptr, FATAL,"Cannot allocate memory");

} // Roll back changes and rethrow on error
  catch (exception&) {
    finalize();
    throw;
  }
  // Catch uncaught exceptions
  catch (...) {
    finalize();
    Throw(UNCAUGHT, "Non-standard exception");
  }

Physics::~Physics() noexcept
//******************************************************************************
//  Destructor
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  finalize();
}

void
Physics::finalize() noexcept
//******************************************************************************
//! \details Single exit point, called implicitly from destructor or explicitly
//!          from anywhere else. Exception safety: no-throw guarantee: never
//!          throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  m_memory->freeEntry(m_particles);  
  if (m_mass) { delete m_mass; m_mass = nullptr; }
  if (m_hydro) { delete m_hydro; m_hydro = nullptr; }
  if (m_mix) { delete m_mix; m_mix = nullptr; }
  if (m_statistics) { delete m_statistics; m_statistics = nullptr; }
  if (m_plot) { delete m_plot; m_plot = nullptr; }
  if (m_glob) { delete m_glob; m_glob = nullptr; }
}
