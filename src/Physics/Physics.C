//******************************************************************************
/*!
  \file      src/Physics/Physics.C
  \author    J. Bakosi
  \date      Sun 15 Sep 2013 10:26:08 PM MDT
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

Physics::Physics(const Base& base)
//******************************************************************************
//  Constructor
//! \param[in]  base     Essentials
//! \author  J. Bakosi
//******************************************************************************
try :
  m_nposition(base.control.get<control::component, control::nposition>()),
  m_ndensity(base.control.get<control::component, control::ndensity>()),
  m_nvelocity(base.control.get<control::component, control::nvelocity>()),
  m_nscalar(base.control.get<control::component, control::nscalar>()),
  m_npar(base.control.get<control::component, control::npar>()),
  m_term(base.control.get<control::incpar, control::term>()),
  m_base(base),
  m_mass(nullptr),
  m_hydro(nullptr),
  m_mix(nullptr),
  m_statistics(nullptr),
  m_glob(nullptr),
  m_stat(nullptr),
  m_particles()
{
  using namespace control;

  //! Echo information on physics
  echo();

  ErrChk(m_base.control.nprop() != 0, ExceptType::FATAL, "No need for physics?");

  // Allocate memory to store all particle properties
  m_particles =
    m_base.memory.newEntry<real>(m_npar * m_base.control.nprop(),
                                 REAL,
                                 SCALAR,
                                 "Particles");

  // Instantiate mass model
  if (m_ndensity) {
    m_mass = new (nothrow) MassType(m_base, m_particles.ptr);
    ErrChk(m_mass != nullptr, ExceptType::FATAL, "Cannot allocate memory");
  }

  // Instantiate hydrodynamics model
  if (m_nvelocity) {
    m_hydro = new (nothrow) HydroType(m_base, m_particles.ptr);
    ErrChk(m_hydro != nullptr, ExceptType::FATAL, "Cannot allocate memory");
  }

  // Instantiate mix model
  if (m_nscalar) {
    m_mix = new (nothrow) MixType(m_base, m_particles.ptr);
    ErrChk(m_mix != nullptr, ExceptType::FATAL, "Cannot allocate memory");
  }

  // Instantiate statistics estimator
  m_statistics = new (nothrow) Statistics(m_base, this);
  ErrChk(m_statistics != nullptr, ExceptType::FATAL,"Cannot allocate memory");

  // Instantiate glob file writer
  m_glob = new (nothrow) GlobWriter(m_base.control.get<io,glob>());
  ErrChk(m_glob != nullptr, ExceptType::FATAL,"Cannot allocate memory");

  // Instantiate statistics plot file writer
  m_stat = new (nothrow) TxtStatWriter(m_base.control.get<io,stats>(), m_statistics);
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
  m_base.memory.freeEntry(m_particles);  
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

  m_base.print.section("Physics",
                  ph.name(m_base.control.get<control::selected,control::physics>()));

  m_base.print.subsection("I/O filenames");
  m_base.print.item("Input", m_base.control.get<control::io,control::input>());
  m_base.print.item("Output", m_base.control.get<control::io,control::output>());
  m_base.print.item("Glob", m_base.control.get<control::io,control::glob>());
  m_base.print.item("Statistics", m_base.control.get<control::io,control::stats>());
  m_base.print.item("PDF", m_base.control.get<control::io,control::pdf>());
  m_base.print.endsubsection();

  m_base.print.subsection("Models");
  m_base.print.item("Position",
               po.name(m_base.control.get<control::selected,control::position>()));
  m_base.print.item("Mass",
               ms.name(m_base.control.get<control::selected,control::mass>()));
  m_base.print.item("Hydrodynamics",
               hy.name(m_base.control.get<control::selected,control::hydro>()));
  m_base.print.item("Internal energy",
               en.name(m_base.control.get<control::selected,control::energy>()));
  m_base.print.item("Material mixing",
               mx.name(m_base.control.get<control::selected,control::mix>()));
  m_base.print.item("Turbulence frequency",
               fr.name(m_base.control.get<control::selected,control::frequency>()));
  m_base.print.item("Material mix rate",
               mr.name(m_base.control.get<control::selected,control::mixrate>()));
  m_base.print.endsubsection();

  m_base.print.subsection("Number of components");
  m_base.print.item("Positions",
               m_base.control.get<control::component,control::nposition>());
  m_base.print.item("Densities",
               m_base.control.get<control::component,control::ndensity>());
  m_base.print.item("Velocities",
               m_base.control.get<control::component,control::nvelocity>());
  m_base.print.item("Scalars",
               m_base.control.get<control::component,control::nscalar>());
  m_base.print.item("Turbulent frequencies",
               m_base.control.get<control::component,control::nfrequency>());
  m_base.print.item("Particles",
               m_base.control.get<control::component,control::npar>());
  m_base.print.endsubsection();

  m_base.print.subsection("Incrementation parameters");
  m_base.print.item("Number of time steps",
               m_base.control.get<control::incpar,control::nstep>());
  m_base.print.item("Terminate time",
               m_base.control.get<control::incpar,control::term>());
  m_base.print.item("Initial time step size",
               m_base.control.get<control::incpar,control::dt>());
  m_base.print.endsubsection();

  m_base.print.subsection("Output intervals");
  m_base.print.item("TTY", m_base.control.get<control::interval,control::tty>());
  m_base.print.item("Dump", m_base.control.get<control::interval,control::dump>());
  m_base.print.item("Glob", m_base.control.get<control::interval,control::glob>());
  m_base.print.item("Statistics", m_base.control.get<control::interval,control::plot>());
  m_base.print.item("PDF", m_base.control.get<control::interval,control::pdf>());
  m_base.print.endsubsection();

  m_base.print.subsection("Statistics");
  m_base.print.vecvecNames<control::stats>(m_base.control,"Requested statistics",true);
  m_base.print.vecvecNames<control::stats>(m_base.control,"Estimated statistics");
  m_base.print.endpart();
}
