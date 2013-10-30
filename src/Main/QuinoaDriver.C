//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.C
  \author    J. Bakosi
  \date      Wed Oct 30 06:59:41 2013
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

  // Parse input deck into m_control, transfer cmdline (no longer needed)
  InputDeckParser inputdeckParser(print, std::move(cmdline), m_control);

  // Create pretty printer for Quinoa
  m_print = std::unique_ptr< QuinoaPrint >( new QuinoaPrint(m_control) );
  m_paradigm = std::unique_ptr< tk::Paradigm >( new tk::Paradigm(print) );
  m_timer = std::unique_ptr< tk::Timer >( new tk::Timer );

  print.endpart();

  // Bundle up essentials
  m_base = std::unique_ptr< Base >(
             new Base(*m_print, *m_paradigm, *m_control, *m_timer) );

  print.part("Factory");

  //! Initialize factories
  initFactories(print);

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
  }
}

void
QuinoaDriver::initFactories(const tk::Print& print)
//******************************************************************************
//  Initialize geometry and physics factories
//! \author  J. Bakosi
//******************************************************************************
{
  // Register geometry types
  ctr::Geometry geometry;
  std::list< std::string > regGeometry;
  regGeometry.push_back(
    geometry.add<AnalyticGeometry>(m_geometryFactory,
                                   ctr::GeometryType::ANALYTIC, *m_base) );
  regGeometry.push_back(
    geometry.add<DiscreteGeometry>(m_geometryFactory,
                                   ctr::GeometryType::DISCRETE, *m_base) );
  print.list("Registered geometry options", regGeometry);

  // Register physics types
  ctr::Physics physics;
  std::list< std::string > regPhysics;
  regPhysics.push_back(
    physics.add<HomMix>(m_physicsFactory,
                        ctr::PhysicsType::HOMOGENEOUS_MIX, *m_base) );
  regPhysics.push_back(
    physics.add<HomHydro>(m_physicsFactory,
                          ctr::PhysicsType::HOMOGENEOUS_HYDRO, *m_base) );
  regPhysics.push_back(
    physics.add<HomRT>(m_physicsFactory,
                       ctr::PhysicsType::HOMOGENEOUS_RAYLEIGH_TAYLOR, *m_base) );
  regPhysics.push_back(
    physics.add<SPINSFlow>(m_physicsFactory,
                           ctr::PhysicsType::SPINSFLOW, *m_base) );
  print.list("Registered physics options", regPhysics);
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
