//******************************************************************************
/*!
  \file      src/MonteCarlo/Physics.C
  \author    J. Bakosi
  \date      Tue Jan 14 09:07:40 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <Factory.h>
#include <Physics.h>
#include <Quinoa/InputDeck/InputDeck.h>
#include <InitPolicy.h>
#include <DirCoeffPolicy.h>
#include <Dirichlet.h>

using quinoa::Physics;

Physics::Physics( const Base& base ) : MonteCarlo( base ),
  m_nposition( base.control.get<ctr::component, ctr::nposition>() ),
  m_ndensity( base.control.get<ctr::component, ctr::ndensity>() ),
  m_nvelocity( base.control.get<ctr::component, ctr::nvelocity>() ),
  m_nscalar( base.control.get<ctr::component, ctr::nscalar>() )
//******************************************************************************
//  Constructor
//! \param[in]  base     Essentials
//! \author  J. Bakosi
//******************************************************************************
{
  //! Initialize factories
  initFactories( print() );

  //! Echo information on physics to be created
  echo();

//   // Instantiate random number generator
//   ctr::RNGType r = m_base.control.get<ctr::selected, ctr::rng>();
//   if (r != ctr::RNGType::NO_RNG) {
//     m_rng = std::unique_ptr<tk::RNG>( m_RNGFactory[r]() );
//   }

//   // Instantiate mass model
//   if (m_ndensity) {
//     ctr::MassType m = control().get<ctr::selected, ctr::mass>();
//     m_mass = std::unique_ptr<Mass>( m_massFactory[m]() );
//   }
// 
//   // Instantiate hydrodynamics model
//   if (m_nvelocity) {
//     ctr::HydroType m = control().get<ctr::selected, ctr::hydro>();
//     m_hydro = std::unique_ptr<Hydro>( m_hydroFactory[m]() );
//   }
// 
//   // Instantiate mix model
//   if (m_nscalar) {
//     ctr::MixType m = control().get<ctr::selected, ctr::mix>();
//     m_mix = std::unique_ptr<Mix>( m_mixFactory[m]() );
//   }
}

void
Physics::initFactories(const tk::Print& print)
//******************************************************************************
//  Initialize factories
//! \author  J. Bakosi
//******************************************************************************
{
//   // Register mass models
//   ctr::Mass mass;
//   std::list< ctr::MassType > regMass;
//   mass.initFactory(m_massFactory, regMass);
//   print.list("Registered mass models", mass, regMass);
// 
//   // Register hydro models
//   ctr::Hydro hydro;
//   std::list< ctr::HydroType > regHydro;
//   hydro.initFactory(m_hydroFactory, regHydro);
//   print.list("Registered hydrodyanmics models", hydro, regHydro);

  // Register mix models
  ctr::Mix mix;
  std::list< ctr::MixType > regMix;
  tk::regist< Dirichlet< InitRaw, DirCoeffConst > >
            ( m_mixFactory, regMix, mix, ctr::MixType::DIRICHLET,
              base(), particles() );
  mix.initFactory(m_mixFactory, regMix);
  print.list("Registered material mix models", mix, regMix);
}

void
Physics::echo()
//******************************************************************************
//  Echo information on physics
//! \author J. Bakosi
//******************************************************************************
{
  const QuinoaPrint& print = this->print();
  const ctr::InputDeck& control = this->control();

  print.endpart();
  print.part( "Problem" );
  print.section( "Title", control.get< ctr::title >() );
  print.Section< ctr::Physics, ctr::selected, ctr::physics >();

  print.subsection( "Output filenames" );
  print.item( "Input", control.get< ctr::cmd, ctr::io, ctr::input >() );
  print.item( "Output", control.get< ctr::cmd, ctr::io, ctr::output >() );
  print.item( "Glob", control.get< ctr::cmd, ctr::io, ctr::glob >() );
  print.item( "Statistics", control.get< ctr::cmd, ctr::io, ctr::stat >() );
  print.item( "PDF", control.get< ctr::cmd, ctr::io, ctr::pdf >() );
  print.endsubsection();

  print.subsection("Selected");
//  print.Item< ctr::RNG, ctr::selected, ctr::rng >();
//   if (control.get< ctr::selected, ctr::rng>() != ctr::RNGType::NO_RNG ) {
//     print.item( "Seed", control.get< ctr::param, ctr::rng, ctr::seed >() );
//   }
  print.Item< ctr::Position, ctr::selected, ctr::position >();
  print.Item< ctr::Mass, ctr::selected, ctr::mass >();
  print.Item< ctr::Hydro, ctr::selected, ctr::hydro >();
  print.Item< ctr::Energy, ctr::selected, ctr::energy >();
  print.Item< ctr::Mix, ctr::selected, ctr::mix >();
  print.Item< ctr::Frequency, ctr::selected, ctr::frequency >();
  print.Item< ctr::MixRate, ctr::selected, ctr::mixrate >();
  print.endsubsection();

  print.subsection( "Number of components" );
  print.Item< ctr::component, ctr::nposition >( "Positions" );
  print.Item< ctr::component, ctr::ndensity >( "Densities" );
  print.Item< ctr::component, ctr::nvelocity >( "Velocities" );
  print.Item< ctr::component, ctr::nscalar >( "Scalars" );
  print.Item< ctr::component, ctr::nfrequency >( "Turbulent frequencies" );
  print.Item< ctr::component, ctr::npar >( "Particles" );
  print.endsubsection();

  print.subsection( "Incrementation parameters" );
  print.item( "Number of time steps", control.get< ctr::incpar,ctr::nstep >() );
  print.item( "Terminate time", control.get< ctr::incpar,ctr::term >() );
  print.item( "Initial time step size", control.get< ctr::incpar, ctr::dt >() );
  print.endsubsection();

  print.subsection( "Output intervals" );
  print.item( "TTY", control.get< ctr::interval, ctr::tty>() );
  print.item( "Dump", control.get< ctr::interval, ctr::dump>() );
  print.item( "Glob", control.get< ctr::interval, ctr::glob >() );
  print.item( "Statistics", control.get< ctr::interval, ctr::plot >() );
  print.item( "PDF", control.get< ctr::interval, ctr::pdf >() );
  print.endsubsection();

  print.subsection( "Statistics" );
  print.RequestedStats( "Requested" );
  print.EstimatedStats( "Estimated" );
  print.endpart();
}
