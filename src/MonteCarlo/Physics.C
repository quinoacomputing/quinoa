//******************************************************************************
/*!
  \file      src/MonteCarlo/Physics.C
  \author    J. Bakosi
  \date      Fri 21 Feb 2014 06:30:52 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <Factory.h>
#include <Physics.h>
#include <Quinoa/InputDeck/InputDeck.h>
#include <Dirichlet.h>
#include <GenDirichlet.h>

using quinoa::Physics;

Physics::Physics( const Base& base ) : MonteCarlo( base ),
  m_nposition( base.control.get<tag::component, tag::position>() ),
  m_ndensity( base.control.get<tag::component, tag::mass>() ),
  m_nvelocity( base.control.get<tag::component, tag::hydro>() ),
  m_nscalar( base.control.get<tag::component, tag::mix>() )
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
  // tk::record<>()
  // print.list< ctr::Mass >( "Registered mass models", m_MassFactory );

  // Register hydro models
  // tk::record<>()
  // print.list< ctr::Hydro >( "Registered hydrod models", m_hydroFactory );

  // Register mix models
  const auto& comp = control().get< tag::component >();
  tk::record< Dirichlet< InitZero, DirCoeffConst > >
            ( m_mixFactory, ctr::MixType::DIRICHLET,
              base(),
              std::cref( particles() ),
              comp.offset< tag::dirichlet >(),
              comp.get< tag::dirichlet >() );
  tk::record< GenDirichlet< InitZero, GenDirCoeffConst > >
            ( m_mixFactory, ctr::MixType::GENDIR,
              base(),
              std::cref( particles() ),
              comp.offset< tag::gendir >(),
              comp.get< tag::gendir >() );
  print.list< ctr::Mix >( "Registered material mix models", m_mixFactory );
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

  print.Section< ctr::MonteCarlo, tag::selected, tag::montecarlo >();

  print.subsection( "Output filenames" );
  print.item( "Input", control.get< tag::cmd, tag::io, tag::input >() );
  print.item( "Output", control.get< tag::cmd, tag::io, tag::output >() );
  print.item( "Glob", control.get< tag::cmd, tag::io, tag::glob >() );
  print.item( "Statistics", control.get< tag::cmd, tag::io, tag::stat >() );
  print.item( "PDF", control.get< tag::cmd, tag::io, tag::pdf >() );
  print.endsubsection();

  print.Model< ctr::Mix, tag::selected, tag::mix >( *m_mix );

  print.subsection( "Increment parameters" );
  print.item( "Number of particles", control.get< tag::incpar, tag::npar >() );
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
