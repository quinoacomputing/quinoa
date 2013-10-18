//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.C
  \author    J. Bakosi
  \date      Fri Oct 11 16:32:16 2013
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

QuinoaDriver::QuinoaDriver(int argc, char** argv, Base& base) :
  m_base(base)
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \param[in] base      Essentials
//! \details   Instantiate geometry, physics, set initial conditions, etc.
//! \author J. Bakosi
//******************************************************************************
{
  // Parse command line
  std::unique_ptr< ctr::CmdLine > cmdline;
  CmdLineParser cmdParser(argc, argv, base, cmdline);

  // Parse input deck
  std::unique_ptr< ctr::InputDeck > inputdeck;
  InputDeckParser inputdeckParser(base, cmdline, inputdeck);

  m_base.print.endpart();
  base.print.part("Problem setup");
  base.print.section("Title", base.control.get<ctr::title>());

  //! Initialize geometry and physics factories
  initFactory();

  // Instantiate geometry object
  ctr::GeometryType g = m_base.control.get<ctr::selected, ctr::geometry>();
  if (g != ctr::GeometryType::NO_GEOMETRY) {
    m_geometry = std::unique_ptr<Geometry>(m_geometryFactory[g]());
  }

  // Instantiate physics object
  ctr::PhysicsType p = m_base.control.get<ctr::selected, ctr::physics>();
  if (p != ctr::PhysicsType::NO_PHYSICS) {
    m_physics = std::unique_ptr<Physics>(m_physicsFactory[p]());
  }

  // Echo 'unspecified' if both geometry and physics are unspecified
  if (!m_geometry && !m_physics) {
    base.print.section("Physics", std::string("unspecified"));
    base.print.section("Geometry", std::string("unspecified"));
    base.print.endpart();
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
    std::bind(boost::factory<AnalyticGeometry*>(), m_base);
  m_geometryFactory[ctr::GeometryType::DISCRETE] =
    std::bind(boost::factory<DiscreteGeometry*>(), m_base);

  // Register physics classes
  m_physicsFactory[ctr::PhysicsType::HOMOGENEOUS_MIX] =
    std::bind(boost::factory<HomMix*>(), m_base);
  m_physicsFactory[ctr::PhysicsType::HOMOGENEOUS_HYDRO] =
    std::bind(boost::factory<HomHydro*>(), m_base);
  m_physicsFactory[ctr::PhysicsType::HOMOGENEOUS_RAYLEIGH_TAYLOR] =
    std::bind(boost::factory<HomRT*>(), m_base);
  m_physicsFactory[ctr::PhysicsType::SPINSFLOW] =
    std::bind(boost::factory<SPINSFlow*>(), m_base);
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
