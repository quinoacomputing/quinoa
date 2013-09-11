//******************************************************************************
/*!
  \file      src/Physics/Physics.C
  \author    J. Bakosi
  \date      Wed Sep 11 17:21:48 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <QuinoaConfig.h>
#include <Memory.h>
#include <Paradigm.h>
#include <QuinoaControl.h>
#include <Physics.h>
#include <Beta.h>
#include <Dirichlet.h>
#include <GenDirichlet.h>
#include <SLM.h>
#include <GLM.h>
#include <Statistics.h>
#include <GlobWriter.h>
#include <TxtStatWriter.h>

using namespace std;
using namespace quinoa;

Physics::Physics(Memory* const memory,
                 Paradigm* const paradigm,
                 const QuinoaControl& control,
                 Timer* const timer,
                 const QuinoaPrinter& print)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object
//! \param[in]  paradigm Parallel programming object
//! \param[in]  control  Control object
//! \param[in]  timer    Timer object
//! \param[in]  print    Quinoa's pretty printer
//! \author  J. Bakosi
//******************************************************************************
try :
  m_nposition(control.get<control::component, control::nposition>()),
  m_ndensity(control.get<control::component, control::ndensity>()),
  m_nvelocity(control.get<control::component, control::nvelocity>()),
  m_nscalar(control.get<control::component, control::nscalar>()),
  m_npar(control.get<control::component, control::npar>()),
  m_term(control.get<control::incpar, control::term>()),
  m_memory(memory),
  m_paradigm(paradigm),
  m_control(control),
  m_print(print),
  m_timer(timer),
  m_mass(nullptr),
  m_hydro(nullptr),
  m_mix(nullptr),
  m_statistics(nullptr),
  m_glob(nullptr),
  m_stat(nullptr),
  m_particles()
{
IGNORE(m_paradigm);

  using namespace control;

  //! Echo information on physics
  echo();

  ErrChk(control.nprop() != 0, ExceptType::FATAL, "No need for physics?");

  // Allocate memory to store all particle properties
  m_particles =
    m_memory->newEntry<real>(m_npar * control.nprop(),
                             REAL,
                             SCALAR,
                             "Particles");

  // Instantiate mass model
  if (m_ndensity) {
    m_mass = new (nothrow) MassType(memory, paradigm, control, m_particles.ptr);
    ErrChk(m_mass != nullptr, ExceptType::FATAL, "Cannot allocate memory");
  }

  // Instantiate hydrodynamics model
  if (m_nvelocity) {
    m_hydro = new (nothrow)
              HydroType(memory, paradigm, control, m_particles.ptr);
    ErrChk(m_hydro != nullptr, ExceptType::FATAL, "Cannot allocate memory");
  }

  // Instantiate mix model
  if (m_nscalar) {
    m_mix = new (nothrow) MixType(memory, paradigm, control, m_particles.ptr);
    ErrChk(m_mix != nullptr, ExceptType::FATAL, "Cannot allocate memory");
  }

  // Instantiate statistics estimator
  m_statistics = new (nothrow) Statistics(memory, paradigm, control, this);
  ErrChk(m_statistics != nullptr, ExceptType::FATAL,"Cannot allocate memory");

  // Instantiate glob file writer
  m_glob = new (nothrow) GlobWriter(control.get<io,glob>());
  ErrChk(m_glob != nullptr, ExceptType::FATAL,"Cannot allocate memory");

  // Instantiate statistics plot file writer
  m_stat = new (nothrow) TxtStatWriter(control.get<io,stats>(), m_statistics);
  ErrChk(m_stat != nullptr, ExceptType::FATAL,"Cannot allocate memory");

} // Roll back changes and rethrow on error
  catch (exception&) {
    finalize();
    throw;
  }
  // Catch uncaught exceptions
  catch (...) {
    finalize();
    Throw(ExceptType::UNCAUGHT, "Non-standard exception");
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
  if (m_glob) { delete m_glob; m_glob = nullptr; }
  if (m_stat) { delete m_stat; m_stat = nullptr; }
}

void
Physics::echo()
//******************************************************************************
//  Echo information on physics
//! \author J. Bakosi
//******************************************************************************
{
  control::Option<select::Physics> ph;
  control::Option<select::Position> po;
  control::Option<select::Mass> ms;
  control::Option<select::Hydro> hy;
  control::Option<select::Energy> en;
  control::Option<select::Mix> mx;
  control::Option<select::Frequency> fr;
  control::Option<select::MixRate> mr;

  m_print.section("Physics",
                  ph.name(m_control.get<control::selected,control::physics>()));

  m_print.subsection("Models");
  m_print.item("Position",
               po.name(m_control.get<control::selected,control::position>()));
  m_print.item("Mass",
               ms.name(m_control.get<control::selected,control::mass>()));
  m_print.item("Hydrodynamics",
               hy.name(m_control.get<control::selected,control::hydro>()));
  m_print.item("Internal energy",
               en.name(m_control.get<control::selected,control::energy>()));
  m_print.item("Material mixing",
               mx.name(m_control.get<control::selected,control::mix>()));
  m_print.item("Turbulence frequency",
               fr.name(m_control.get<control::selected,control::frequency>()));
  m_print.item("Material mix rate",
               mr.name(m_control.get<control::selected,control::mixrate>()));

  m_print.subsection("Incrementation parameters");
  m_print.item("Number of time steps",
               m_control.get<control::incpar,control::nstep>());
  m_print.item("Terminate time",
               m_control.get<control::incpar,control::term>());
  m_print.item("Initial time step size",
               m_control.get<control::incpar,control::dt>());

  m_print.subsection("Number of components");
  m_print.item("Positions",
               m_control.get<control::component,control::nposition>());
  m_print.item("Densities",
               m_control.get<control::component,control::ndensity>());
  m_print.item("Velocities",
               m_control.get<control::component,control::nvelocity>());
  m_print.item("Scalars",
               m_control.get<control::component,control::nscalar>());
  m_print.item("Turbulent frequencies",
               m_control.get<control::component,control::nfrequency>());
  m_print.item("Particles",
               m_control.get<control::component,control::npar>());

  m_print.subsection("Output intervals");
  m_print.item("TTY", m_control.get<control::interval,control::tty>());
  m_print.item("Dump", m_control.get<control::interval,control::dump>());
  m_print.item("Glob", m_control.get<control::interval,control::glob>());
  m_print.item("Statistics", m_control.get<control::interval,control::plot>());
  m_print.item("PDF", m_control.get<control::interval,control::pdf>());

  m_print.subsection("I/O filenames");
  m_print.item("Input", m_control.get<control::io,control::input>());
  m_print.item("Output", m_control.get<control::io,control::output>());
  m_print.item("Glob", m_control.get<control::io,control::glob>());
  m_print.item("Statistics", m_control.get<control::io,control::stats>());
  m_print.item("PDF", m_control.get<control::io,control::pdf>());

  m_print.subsection("Statistics");
  //m_control.echoVecVecNames<control::stats>("Requested statistics",true);
  //m_control.echoVecVecNames<control::stats>("Estimated statistics");
}
