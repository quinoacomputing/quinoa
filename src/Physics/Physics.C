//******************************************************************************
/*!
  \file      src/Physics/Physics.C
  \author    J. Bakosi
  \date      Tue Oct 29 15:37:09 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Physics.h>
#include <Quinoa/InputDeck/InputDeck.h>
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

  // Instantiate random number generator
  ctr::RNGType r = m_base.control.get<ctr::selected, ctr::rng>();
  if (r != ctr::RNGType::NO_RNG) {
    m_rng = std::unique_ptr<tk::RNG>( m_RNGFactory[r]() );
  }

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
  // Register random number generators
  ctr::RNG rng;
  unsigned int seed = m_base.control.get<ctr::param, ctr::rng, ctr::seed>();
  std::list< std::string > registeredRNG;
  rng.initFactory(m_RNGFactory, registeredRNG, m_base.paradigm.nthreads(), seed);

  // Register mass models
  ctr::Mass mass;
  std::list< std::string > registeredMass;
  mass.initFactory(m_massFactory, registeredMass);

  // Register hydro models
  ctr::Hydro hydro;
  std::list< std::string > registeredHydro;
  hydro.initFactory(m_hydroFactory, registeredHydro);

  // Register mix models
  ctr::Mix mix;
  std::list< std::string > registeredMix;
  mix.initFactory(m_mixFactory, registeredMix);
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

  m_base.print.subsection("Selected");
  m_base.print.item<ctr::RNG, ctr::selected, ctr::rng>();
  if (m_base.control.get<ctr::selected, ctr::rng>() != ctr::RNGType::NO_RNG) {
    m_base.print.item("Seed",
                       m_base.control.get<ctr::param, ctr::rng, ctr::seed>());
  }
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
