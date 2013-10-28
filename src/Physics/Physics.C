//******************************************************************************
/*!
  \file      src/Physics/Physics.C
  \author    J. Bakosi
  \date      Mon Oct 28 13:54:30 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Physics.h>
#include <Quinoa/InputDeck/InputDeck.h>
#include <Mass/Beta/Beta.h>
#include <Mix/Dirichlet/Dirichlet.h>
#include <Mix/GenDirichlet/GenDirichlet.h>
#include <Hydro/SLM/SLM.h>
#include <Hydro/GLM/GLM.h>
#include <GlobWriter.h>
#include <TxtStatWriter.h>

using namespace quinoa;

Physics::Physics(const Base& base) :
  m_nposition(base.control.get<ctr::component, ctr::nposition>()),
  m_ndensity(base.control.get<ctr::component, ctr::ndensity>()),
  m_nvelocity(base.control.get<ctr::component, ctr::nvelocity>()),
  m_nscalar(base.control.get<ctr::component, ctr::nscalar>()),
  m_npar(base.control.get<ctr::component, ctr::npar>()),
  m_term(base.control.get<ctr::incpar, ctr::term>()),
  m_base(base),
  m_particles(new tk::real [m_npar * base.control.nprop()]),
  m_statistics(base, m_particles.get()),
  m_glob(base.control.get<ctr::cmd, ctr::io, ctr::glob>()),
  m_stat(base.control.get<ctr::cmd, ctr::io, ctr::stat>(), m_statistics)
//******************************************************************************
//  Constructor
//! \param[in]  base     Essentials
//! \author  J. Bakosi
//******************************************************************************
{
  //! Echo information on physics to be created
  echo();

  //! Initialize factories
  initFactories();

  // Instantiate mass model
  if (m_ndensity) {
    ctr::MassType m = m_base.control.get<ctr::selected, ctr::mass>();
    m_mass = std::unique_ptr<Mass>( m_massFactory[m]() );
  }

  // Instantiate hydrodynamics model
  if (m_nvelocity) {
    ctr::HydroType m = m_base.control.get<ctr::selected, ctr::hydro>();
    m_hydro = std::unique_ptr<Hydro>( m_hydroFactory[m]() );
  }

  // Instantiate mix model
  if (m_nscalar) {
    ctr::MixType m = m_base.control.get<ctr::selected, ctr::mix>();
    m_mix = std::unique_ptr<Mix>( m_mixFactory[m]() );
  }
}

void
Physics::initFactories()
//******************************************************************************
//  Initialize factories
//! \author  J. Bakosi
//******************************************************************************
{
  // TODO: Parse in and from control
  // TODO: Echo number of registered stuff
  unsigned int seed = 0;

  // Register random number generators
  ctr::RNG rng;
  rng.initFactory(m_RNGFactory, m_base.paradigm.nthreads(), seed);

  // Register mass models
  ctr::Mass mass;
  mass.initFactory(m_massFactory);

  // Register hydro models
  ctr::Hydro hydro;
  hydro.initFactory(m_hydroFactory);

  // Register mix models
  ctr::Mix mix;
  mix.initFactory(m_mixFactory);
}

void
Physics::echo()
//******************************************************************************
//  Echo information on physics
//! \author J. Bakosi
//******************************************************************************
{
  m_base.print.section<ctr::Physics, ctr::selected, ctr::physics>();

  m_base.print.subsection("Output filenames");
  m_base.print.item( "Input",
                     m_base.control.get< ctr::cmd, ctr::io, ctr::input >() );
  m_base.print.item( "Output",
                     m_base.control.get< ctr::cmd, ctr::io, ctr::output >() );
  m_base.print.item( "Glob",
                     m_base.control.get< ctr::cmd, ctr::io, ctr::glob >() );
  m_base.print.item( "Statistics",
                     m_base.control.get< ctr::cmd, ctr::io, ctr::stat >() );
  m_base.print.item( "PDF",
                     m_base.control.get< ctr::cmd, ctr::io, ctr::pdf >() );
  m_base.print.endsubsection();

  m_base.print.subsection("Models");
  m_base.print.item<ctr::Position, ctr::selected, ctr::position>();
  m_base.print.item<ctr::Mass, ctr::selected, ctr::mass>();
  m_base.print.item<ctr::Hydro, ctr::selected, ctr::hydro>();
  m_base.print.item<ctr::Energy, ctr::selected, ctr::energy>();
  m_base.print.item<ctr::Mix, ctr::selected, ctr::mix>();
  m_base.print.item<ctr::Frequency, ctr::selected, ctr::frequency>();
  m_base.print.item<ctr::MixRate, ctr::selected, ctr::mixrate>();
  m_base.print.endsubsection();

  m_base.print.subsection("Number of components");
  m_base.print.item<ctr::component,ctr::nposition>("Positions");
  m_base.print.item<ctr::component,ctr::ndensity>("Densities");
  m_base.print.item<ctr::component,ctr::nvelocity>("Velocities");
  m_base.print.item<ctr::component,ctr::nscalar>("Scalars");
  m_base.print.item<ctr::component,ctr::nfrequency>("Turbulent frequencies");
  m_base.print.item<ctr::component,ctr::npar>("Particles");
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
  m_base.print.requestedStats("Requested");
  m_base.print.estimatedStats("Estimated");
  m_base.print.endpart();
}
