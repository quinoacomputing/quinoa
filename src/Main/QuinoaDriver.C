//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.C
  \author    J. Bakosi
  \date      Fri Oct 18 11:31:17 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <QuinoaDriver.h>
#include <Quinoa/InputDeck/Parser.h>
#include <Quinoa/CmdLine/Parser.h>
#include <AnalyticGeometry.h>
#include <DiscreteGeometry.h>
#include <HomMix/HomMix.h>
#include <HomHydro/HomHydro.h>
#include <HomRT/HomRT.h>
#include <SPINSFlow/SPINSFlow.h>

using namespace quinoa;

QuinoaDriver::QuinoaDriver(int argc, char** argv, const tk::Print& print)
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \param[in] print     Simple pretty printer
//! \details   Instantiate geometry, physics, set initial conditions, etc.
//! \author J. Bakosi
//******************************************************************************
{
  // Parse command line into cmdline
  std::unique_ptr< ctr::CmdLine > cmdline;
  CmdLineParser cmdParser(argc, argv, print, cmdline);

  // Parse input deck into m_control
  InputDeckParser inputdeckParser(print, cmdline, m_control);

  // Create pretty printer for Quinoa
  m_print = std::unique_ptr< QuinoaPrint >( new QuinoaPrint(m_control) );
  m_paradigm = std::unique_ptr< tk::Paradigm >( new tk::Paradigm(print) );
  m_timer = std::unique_ptr< tk::Timer >( new tk::Timer );

  // Bundle up essentials
  m_base = std::unique_ptr< Base >(
             new Base(*m_print, *m_paradigm, *m_control, *m_timer) );

  print.endpart();
  print.part("Problem setup");
  print.section("Title", m_control->get<ctr::title>());

  //! Initialize geometry and physics factories
  initFactory();

  // Instantiate geometry object
  ctr::GeometryType g = m_control->get<ctr::selected, ctr::geometry>();
  if (g != ctr::GeometryType::NO_GEOMETRY) {
    m_geometry = std::unique_ptr<Geometry>(m_geometryFactory[g]());
  }

  // Instantiate physics object
  ctr::PhysicsType p = m_control->get<ctr::selected, ctr::physics>();
  if (p != ctr::PhysicsType::NO_PHYSICS) {
    m_physics = std::unique_ptr<Physics>(m_physicsFactory[p]());
  }

  // Echo 'unspecified' if both geometry and physics are unspecified
  if (!m_geometry && !m_physics) {
    print.section("Physics", std::string("unspecified"));
    print.section("Geometry", std::string("unspecified"));
    print.endpart();
  }
}

void
QuinoaDriver::initFactory()
//******************************************************************************
//  Initialize geometry and physics factories
//! \author  J. Bakosi
//******************************************************************************
{
  // Register geometry classes
  m_geometryFactory[ctr::GeometryType::ANALYTIC] =
    std::bind(boost::factory<AnalyticGeometry*>(), *m_base);
  m_geometryFactory[ctr::GeometryType::DISCRETE] =
    std::bind(boost::factory<DiscreteGeometry*>(), *m_base);

  // Register physics classes
  m_physicsFactory[ctr::PhysicsType::HOMOGENEOUS_MIX] =
    std::bind(boost::factory<HomMix*>(), *m_base);
  m_physicsFactory[ctr::PhysicsType::HOMOGENEOUS_HYDRO] =
    std::bind(boost::factory<HomHydro*>(), *m_base);
  m_physicsFactory[ctr::PhysicsType::HOMOGENEOUS_RAYLEIGH_TAYLOR] =
    std::bind(boost::factory<HomRT*>(), *m_base);
  m_physicsFactory[ctr::PhysicsType::SPINSFLOW] =
    std::bind(boost::factory<SPINSFlow*>(), *m_base);
}

void
QuinoaDriver::execute() const
//******************************************************************************
//  Execute
//! \author J. Bakosi
//******************************************************************************
{
  //! Initialize and execute geometry (if any)
  if (m_geometry) {
    m_geometry->init();
    m_geometry->fill();
  }

  //! Initialize and execute physics (if any)
  if (m_physics) {
    m_physics->init();
    m_physics->solve();
  }
}
