//******************************************************************************
/*!
  \file      src/Physics/Physics.C
  \author    J. Bakosi
  \date      Fri May 10 17:28:44 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <type_traits>

#include <QuinoaConfig.h>
#include <Memory.h>
#include <Paradigm.h>
#include <Control.h>
#include <Physics.h>
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
  m_nvelocity(control->get<control::NVELOCITY>()),
  m_nscalar(control->get<control::NSCALAR>()),
  m_npar(control->get<control::NPAR>()),
  m_term(control->get<control::TERM>()),
  m_memory(memory),
  m_paradigm(paradigm),
  m_control(control),
  m_timer(timer),
  //m_mix(nullptr),
  m_hydro(nullptr),
  m_statistics(nullptr),
  m_glob(nullptr),
  m_plot(nullptr),
  m_particles()
{

  // Compute number of particle properties
  int nprop = m_nvelocity + m_nscalar;
  // Compute offset where scalars start
  int offset = m_npar * m_nvelocity;

  // Allocate memory to store all particle properties
  m_particles =
    m_memory->newEntry<real>(m_npar*nprop, REAL, SCALAR, "Particles");  

//   // Instantiate selected hydrodynamics model
//   switch (control->get<control::HYDRO>()) {
// 
//     case control::HydroType::NO_HYDRO :
//       Throw(FATAL, "No hydro model selected");
//       break;
// 
//     case control::HydroType::SLM :
//       m_hydro = new (nothrow) SimplifiedLangevin(memory,
//                                                  paradigm,
//                                                  control,
//                                                  m_particles + 0);
//       ErrChk(m_hydro != nullptr, FATAL, "Cannot allocate memory");
//       break;
// 
//     case control::HydroType::GLM :
//       m_hydro = new (nothrow) GeneralizedLangevin(memory,
//                                                   paradigm,
//                                                   control,
//                                                   m_particles + 0);
//       ErrChk(m_hydro != nullptr, FATAL, "Cannot allocate memory");
//       break;
// 
//     default :
//       Throw(FATAL, "Hydro model not implemented");
//   }

  // Instantiate selected mix model
  m_mix = new (nothrow)
    Mix<MixType>(memory, paradigm, control, m_nscalar, m_particles + offset);
  ErrChk(m_mix != nullptr, FATAL, "Cannot allocate memory");

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
  //if (m_mix) { delete m_mix; m_mix = nullptr; }
  if (m_hydro) { delete m_hydro; m_hydro = nullptr; }
  if (m_statistics) { delete m_statistics; m_statistics = nullptr; }
  if (m_plot) { delete m_plot; m_plot = nullptr; }
  if (m_glob) { delete m_glob; m_glob = nullptr; }
}
