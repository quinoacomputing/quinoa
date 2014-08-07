//******************************************************************************
/*!
  \file      src/MonteCarlo/Physics.C
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 04:34:50 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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

Physics::Physics( const Base& b ) : MonteCarlo( b ),
  m_nposition( b.control.get<tag::component, tag::position>() ),
  m_ndensity( b.control.get<tag::component, tag::mass>() ),
  m_nvelocity( b.control.get<tag::component, tag::hydro>() ),
  m_nscalar( b.control.get<tag::component, tag::mix>() )
//******************************************************************************
//  Constructor
//! \param[in]  b     Essentials
//! \author  J. Bakosi
//******************************************************************************
{
  //! Initialize factories
  initFactories();

//   // Instantiate mass model
//   if (m_ndensity) {
//     m_mass = tk::instantiate( m_massFactory,
//                               control().get< tag::selected, tag::mass >() );
//   }
// 
//   // Instantiate hydrodynamics model
//   if (m_nvelocity) {
//     m_hydro = tk::instantiate( m_hydroFactory,
//                                control().get< tag::selected, tag::hydro >() );
//   }
// 
//   // Instantiate mix model
//   if (m_nscalar) {
//     m_mix = tk::instantiate( m_mixFactory,
//                              control().get< tag::selected, tag::mix >() );
//   }
// 
  //! Echo information on physics to be created
  echo();
}

void
Physics::initFactories()
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

//   // Register mix models
//   const auto& comp = control().get< tag::component >();
//   tk::record< Dirichlet< InitZero, DirCoeffConst > >
//             ( m_mixFactory, ctr::MixType::DIRICHLET,
//               base(),
//               std::cref( particles() ),
//               comp.offset< tag::dirichlet >(),
//               comp.get< tag::dirichlet >() );
//   tk::record< GenDirichlet< InitZero, GenDirCoeffConst > >
//             ( m_mixFactory, ctr::MixType::GENDIR,
//               base(),
//               std::cref( particles() ),
//               comp.offset< tag::gendir >(),
//               comp.get< tag::gendir >() );
//   print().list< ctr::Mix >( "Registered material mix models", m_mixFactory );
}

void
Physics::echo()
//******************************************************************************
//  Echo information on physics
//! \author J. Bakosi
//******************************************************************************
{
//   print().endpart();
//   print().part( "Problem" );
// 
//   print().section( "Title", control().get< tag::title >() );
// 
//   echoRNGs();
// 
//   print().Section< ctr::MonteCarlo, tag::selected, tag::montecarlo >();
// 
//   echoIO();
// 
//   print().Model< ctr::Mix, tag::selected, tag::mix >( *m_mix );
// 
//   echoIncpar();
// 
//   echoIntervals();
// 
//   echoStatistics();
// 
//   print().endpart();
}
