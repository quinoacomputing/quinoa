//******************************************************************************
/*!
  \file      src/Physics/Physics.C
  \author    J. Bakosi
  \date      Mon 04 Nov 2013 10:14:59 PM MST
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
  //! Initialize factories
  initFactories(m_base.print);

  //! Echo information on physics to be created
  echo();

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
Physics::initFactories(const tk::Print& print)
//******************************************************************************
//  Initialize factories
//! \author  J. Bakosi
//******************************************************************************
{
  // Register random number generators
  ctr::RNG rng;
  std::list< std::string > regRNG;
  //auto& seed = m_base.control.get<ctr::param, ctr::rng, ctr::seed>();
  //rng.initFactory(m_RNGFactory, regRNG, m_base.paradigm.nthreads(), seed);
  print.list("Registered random number generators", regRNG);

  // Register mass models
  ctr::Mass mass;
  std::list< std::string > regMass;
  mass.initFactory(m_massFactory, regMass);
  print.list("Registered mass models", regMass);

  // Register hydro models
  ctr::Hydro hydro;
  std::list< std::string > regHydro;
  hydro.initFactory(m_hydroFactory, regHydro);
  print.list("Registered hydrodyanmics models", regHydro);

  // Register mix models
  ctr::Mix mix;
  std::list< std::string > regMix;
  mix.initFactory(m_mixFactory, regMix);
  print.list("Registered material mix models", regMix);
}

void
Physics::echo()
//******************************************************************************
//  Echo information on physics
//! \author J. Bakosi
//******************************************************************************
{
  const QuinoaPrint& print = m_base.print;
  const ctr::InputDeck& control = m_base.control;

  print.endpart();
  print.part("Problem");
  print.section("Title", control.get<ctr::title>());
  print.section<ctr::Physics, ctr::selected, ctr::physics>();

  print.subsection("Output filenames");
  print.item( "Input", control.get< ctr::cmd, ctr::io, ctr::input >() );
  print.item( "Output", control.get< ctr::cmd, ctr::io, ctr::output >() );
  print.item( "Glob", control.get< ctr::cmd, ctr::io, ctr::glob >() );
  print.item( "Statistics", control.get< ctr::cmd, ctr::io, ctr::stat >() );
  print.item( "PDF", control.get< ctr::cmd, ctr::io, ctr::pdf >() );
  print.endsubsection();

  print.subsection("Selected");
  print.item<ctr::RNG, ctr::selected, ctr::rng>();
  if (control.get<ctr::selected, ctr::rng>() != ctr::RNGType::NO_RNG) {
    print.item("Seed", control.get<ctr::param, ctr::rng, ctr::seed>());
  }
  print.item<ctr::Position, ctr::selected, ctr::position>();
  print.item<ctr::Mass, ctr::selected, ctr::mass>();
  print.item<ctr::Hydro, ctr::selected, ctr::hydro>();
  print.item<ctr::Energy, ctr::selected, ctr::energy>();
  print.item<ctr::Mix, ctr::selected, ctr::mix>();
  print.item<ctr::Frequency, ctr::selected, ctr::frequency>();
  print.item<ctr::MixRate, ctr::selected, ctr::mixrate>();
  print.endsubsection();

  print.subsection("Number of components");
  print.item<ctr::component,ctr::nposition>("Positions");
  print.item<ctr::component,ctr::ndensity>("Densities");
  print.item<ctr::component,ctr::nvelocity>("Velocities");
  print.item<ctr::component,ctr::nscalar>("Scalars");
  print.item<ctr::component,ctr::nfrequency>("Turbulent frequencies");
  print.item<ctr::component,ctr::npar>("Particles");
  print.endsubsection();

  print.subsection("Incrementation parameters");
  print.item("Number of time steps", control.get<ctr::incpar,ctr::nstep>());
  print.item("Terminate time", control.get<ctr::incpar,ctr::term>());
  print.item("Initial time step size", control.get<ctr::incpar,ctr::dt>());
  print.endsubsection();

  print.subsection("Output intervals");
  print.item("TTY", control.get<ctr::interval,ctr::tty>());
  print.item("Dump", control.get<ctr::interval,ctr::dump>());
  print.item("Glob", control.get<ctr::interval,ctr::glob>());
  print.item("Statistics", control.get<ctr::interval,ctr::plot>());
  print.item("PDF", control.get<ctr::interval,ctr::pdf>());
  print.endsubsection();

  print.subsection("Statistics");
  print.requestedStats("Requested");
  print.estimatedStats("Estimated");
  print.endpart();
}
