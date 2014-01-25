//******************************************************************************
/*!
  \file      src/MonteCarlo/Physics.C
  \author    J. Bakosi
  \date      Sat 25 Jan 2014 03:53:26 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <Factory.h>
#include <Physics.h>
#include <Quinoa/InputDeck/InputDeck.h>
#include <Dirichlet.h>

using quinoa::Physics;

Physics::Physics( const Base& base ) : MonteCarlo( base ),
  m_nposition( base.control.get<tag::component, tag::nposition>() ),
  m_ndensity( base.control.get<tag::component, tag::ndensity>() ),
  m_nvelocity( base.control.get<tag::component, tag::nvelocity>() ),
  m_nscalar( base.control.get<tag::component, tag::nscalar>() )
//******************************************************************************
//  Constructor
//! \param[in]  base     Essentials
//! \author  J. Bakosi
//******************************************************************************
{
  //! Initialize factories
  initFactories( print() );

  // Instantiate mass model
  if (m_ndensity) {
    ctr::MassType m = control().get<tag::selected, tag::mass>();
    m_mass = std::unique_ptr< Model >( m_massFactory[m]() );
  }

  // Instantiate hydrodynamics model
  if (m_nvelocity) {
    ctr::HydroType m = control().get<tag::selected, tag::hydro>();
    m_hydro = std::unique_ptr< Model >( m_hydroFactory[m]() );
  }

  // Instantiate mix model
  if (m_nscalar) {
    ctr::MixType m = control().get<tag::selected, tag::mix>();
    m_mix = std::unique_ptr< Model >( m_mixFactory[m]() );
  }

  //! Echo information on physics to be created
  echo();
}

void
Physics::initFactories(const tk::Print& print)
//******************************************************************************
//  Initialize factories
//! \author  J. Bakosi
//******************************************************************************
{
  // Register mass models
  ctr::Mass mass;
  std::list< ctr::MassType > regMass;
  // tk::regist()
  print.list("Registered mass models", mass, regMass);

  // Register hydro models
  ctr::Hydro hydro;
  std::list< ctr::HydroType > regHydro;
  // tk::regist()
  print.list("Registered hydrodyanmics models", hydro, regHydro);

  // Register mix models
  ctr::Mix mix;
  std::list< ctr::MixType > regMix;
  tk::regist< Dirichlet< InitRaw, DirCoeffConst > >
            ( m_mixFactory, regMix, mix, ctr::MixType::DIRICHLET,
              base(), std::cref(particles()) );
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

  print.section( "Title", control.get< tag::title >() );

  print.section("Random number generators");

  print.MKLParams( control.get< tag::selected, tk::tag::rng >(),
                   control.get< tag::param, tk::tag::mklrng >() );
  print.RNGSSEParams( control.get< tag::selected, tk::tag::rng >(),
                      control.get< tag::param, tk::tag::rngsse >() );

  print.Section< ctr::Physics, tag::selected, tag::physics >();

  print.subsection( "Output filenames" );
  print.item( "Input", control.get< tag::cmd, tag::io, tag::input >() );
  print.item( "Output", control.get< tag::cmd, tag::io, tag::output >() );
  print.item( "Glob", control.get< tag::cmd, tag::io, tag::glob >() );
  print.item( "Statistics", control.get< tag::cmd, tag::io, tag::stat >() );
  print.item( "PDF", control.get< tag::cmd, tag::io, tag::pdf >() );
  print.endsubsection();

//   print.Item< ctr::Position, tag::selected, tag::position >();
//   print.Item< ctr::Mass, tag::selected, tag::mass >();
//   print.Item< ctr::Hydro, tag::selected, tag::hydro >();
//   print.Item< ctr::Energy, tag::selected, tag::energy >();
  print.Model< ctr::Mix, tag::selected, tag::mix >( m_mix.get() );
//   print.Item< ctr::Frequency, tag::selected, tag::frequency >();
//   print.Item< ctr::MixRate, tag::selected, tag::mixrate >();
  print.endsubsection();

  print.subsection( "Number of components" );
//   print.Item< tag::component, tag::nposition >( "Positions" );
//   print.Item< tag::component, tag::ndensity >( "Densities" );
//   print.Item< tag::component, tag::nvelocity >( "Velocities" );
//   print.Item< tag::component, tag::nfrequency >( "Turbulent frequencies" );
  print.Item< tag::component, tag::npar >( "Particles" );
  print.endsubsection();

  print.subsection( "Increment parameters" );
  print.item( "Number of time steps", control.get< tag::incpar, tag::nstep >() );
  print.item( "Terminate time", control.get< tag::incpar, tag::term >() );
  print.item( "Initial time step size", control.get< tag::incpar, tag::dt >() );
  print.endsubsection();

  print.subsection( "Output intervals" );
  print.item( "TTY", control.get< tag::interval, tag::tty>() );
  print.item( "Dump", control.get< tag::interval, tag::dump>() );
  print.item( "Glob", control.get< tag::interval, tag::glob >() );
  print.item( "Statistics", control.get< tag::interval, tag::plot >() );
  print.item( "PDF", control.get< tag::interval, tag::pdf >() );
  print.endsubsection();

  print.subsection( "Statistics" );
  print.RequestedStats( "Requested" );
  print.EstimatedStats( "Estimated" );
  print.endpart();
}
