//******************************************************************************
/*!
  \file      src/MonteCarlo/Physics.C
  \author    J. Bakosi
  \date      Mon 24 Feb 2014 09:01:02 PM MST
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
    m_mass = tk::instantiate( m_massFactory,
                              control().get< tag::selected, tag::mass >() );
  }

  // Instantiate hydrodynamics model
  if (m_nvelocity) {
    m_hydro = tk::instantiate( m_hydroFactory,
                               control().get< tag::selected, tag::hydro >() );
  }

  // Instantiate mix model
  if (m_nscalar) {
    m_mix = tk::instantiate( m_mixFactory,
                             control().get< tag::selected, tag::mix >() );
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
  print().endpart();
  print().part( "Problem" );

  print().section( "Title", control().get< tag::title >() );

  echoRNGs();

  print().Section< ctr::MonteCarlo, tag::selected, tag::montecarlo >();

  echoIO();

  print().Model< ctr::Mix, tag::selected, tag::mix >( *m_mix );

  echoIncpar();

  echoIntervals();

  echoStatistics();

  print().endpart();
}
