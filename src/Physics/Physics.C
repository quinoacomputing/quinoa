//******************************************************************************
/*!
  \file      src/Physics/Physics.C
  \author    J. Bakosi
  \date      Fri Sep 27 09:06:17 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

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
#include <GlobWriter.h>
#include <TxtStatWriter.h>

using namespace std;
using namespace quinoa;

Physics::Physics(const Base& base) :
  m_nposition(base.control.get<ctr::component, ctr::nposition>()),
  m_ndensity(base.control.get<ctr::component, ctr::ndensity>()),
  m_nvelocity(base.control.get<ctr::component, ctr::nvelocity>()),
  m_nscalar(base.control.get<ctr::component, ctr::nscalar>()),
  m_npar(base.control.get<ctr::component, ctr::npar>()),
  m_term(base.control.get<ctr::incpar, ctr::term>()),
  m_base(base),
  m_particles(new real [m_npar * base.control.nprop()]),
  m_statistics(base, m_particles.get()),
  m_glob(base.control.get<ctr::io, ctr::glob>()),
  m_stat(base.control.get<ctr::io, ctr::stats>(), m_statistics)
//******************************************************************************
//  Constructor
//! \param[in]  base     Essentials
//! \author  J. Bakosi
//******************************************************************************
{
  //! Echo information on physics to be created
  echo();

  //! Initialize model factories
  initFactory();

  // Instantiate mass model
  if (m_ndensity) {
    sel::MassType m = m_base.control.get<ctr::selected, ctr::mass>();
    m_mass = std::unique_ptr<Mass>(m_massFactory[m]());
  }

  // Instantiate hydrodynamics model
  if (m_nvelocity) {
    sel::HydroType m = m_base.control.get<ctr::selected, ctr::hydro>();
    m_hydro = std::unique_ptr<Hydro>(m_hydroFactory[m]());
  }

  // Instantiate mix model
  if (m_nscalar) {
    sel::MixType m = m_base.control.get<ctr::selected, ctr::mix>();
    m_mix = std::unique_ptr<Mix>(m_mixFactory[m]());
  }
}

void
Physics::initFactory()
//******************************************************************************
//  Initialize model factory
//! \author  J. Bakosi
//******************************************************************************
{
  // Get raw pointer to particle properties
  real* const p = m_particles.get();

  // Register mass models
  m_massFactory[sel::MassType::BETA] =
    std::bind(boost::factory<Beta*>(), m_base, p);

  // Register hydro models
  m_hydroFactory[sel::HydroType::SLM] =
    std::bind(boost::factory<SimplifiedLangevin*>(), m_base, p);
  m_hydroFactory[sel::HydroType::GLM] =
    std::bind(boost::factory<GeneralizedLangevin*>(), m_base, p);

  // Register mix models
  m_mixFactory[sel::MixType::DIRICHLET] =
    std::bind(boost::factory<Dirichlet*>(), m_base, p);
  m_mixFactory[sel::MixType::GENERALIZED_DIRICHLET] =
    std::bind(boost::factory<GeneralizedDirichlet*>(), m_base, p);
}

void
Physics::echo()
//******************************************************************************
//  Echo information on physics
//! \author J. Bakosi
//******************************************************************************
{
  ctr::Option<sel::Physics> ph;
  ctr::Option<sel::Position> po;
  ctr::Option<sel::Mass> ms;
  ctr::Option<sel::Hydro> hy;
  ctr::Option<sel::Energy> en;
  ctr::Option<sel::Mix> mx;
  ctr::Option<sel::Frequency> fr;
  ctr::Option<sel::MixRate> mr;

  m_base.print.section("Physics",
                  ph.name(m_base.control.get<ctr::selected,ctr::physics>()));

  m_base.print.subsection("I/O filenames");
  m_base.print.item("Input", m_base.control.get<ctr::io,ctr::input>());
  m_base.print.item("Output", m_base.control.get<ctr::io,ctr::output>());
  m_base.print.item("Glob", m_base.control.get<ctr::io,ctr::glob>());
  m_base.print.item("Statistics", m_base.control.get<ctr::io,ctr::stats>());
  m_base.print.item("PDF", m_base.control.get<ctr::io,ctr::pdf>());
  m_base.print.endsubsection();

  m_base.print.subsection("Models");
  m_base.print.option<sel::Position, ctr::selected, ctr::position>();
  m_base.print.option<sel::Mass, ctr::selected, ctr::mass>();
  m_base.print.option<sel::Hydro, ctr::selected, ctr::hydro>();
  m_base.print.option<sel::Energy, ctr::selected, ctr::energy>();
  m_base.print.option<sel::Mix, ctr::selected, ctr::mix>();
  m_base.print.option<sel::Frequency, ctr::selected, ctr::frequency>();
  m_base.print.option<sel::MixRate, ctr::selected, ctr::mixrate>();
  m_base.print.endsubsection();

  m_base.print.subsection("Number of components");
  m_base.print.item("Positions",
               m_base.control.get<ctr::component,ctr::nposition>());
  m_base.print.item("Densities",
               m_base.control.get<ctr::component,ctr::ndensity>());
  m_base.print.item("Velocities",
               m_base.control.get<ctr::component,ctr::nvelocity>());
  m_base.print.item("Scalars",
               m_base.control.get<ctr::component,ctr::nscalar>());
  m_base.print.item("Turbulent frequencies",
               m_base.control.get<ctr::component,ctr::nfrequency>());
  m_base.print.item("Particles",
               m_base.control.get<ctr::component,ctr::npar>());
  m_base.print.endsubsection();

  m_base.print.subsection("Incrementation parameters");
  m_base.print.item("Number of time steps",
               m_base.control.get<ctr::incpar,ctr::nstep>());
  m_base.print.item("Terminate time",
               m_base.control.get<ctr::incpar,ctr::term>());
  m_base.print.item("Initial time step size",
               m_base.control.get<ctr::incpar,ctr::dt>());
  m_base.print.endsubsection();

  m_base.print.subsection("Output intervals");
  m_base.print.item("TTY", m_base.control.get<ctr::interval,ctr::tty>());
  m_base.print.item("Dump", m_base.control.get<ctr::interval,ctr::dump>());
  m_base.print.item("Glob", m_base.control.get<ctr::interval,ctr::glob>());
  m_base.print.item("Statistics", m_base.control.get<ctr::interval,ctr::plot>());
  m_base.print.item("PDF", m_base.control.get<ctr::interval,ctr::pdf>());
  m_base.print.endsubsection();

  m_base.print.subsection("Statistics");
  m_base.print.vecvecNames<ctr::stats>(m_base.control,"Requested statistics",true);
  m_base.print.vecvecNames<ctr::stats>(m_base.control,"Estimated statistics");
  m_base.print.endpart();
}
