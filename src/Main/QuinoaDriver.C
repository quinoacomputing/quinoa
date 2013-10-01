//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.C
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 10:26:00 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************

#include <iostream>

#include <boost/functional/factory.hpp>

#include <QuinoaDriver.h>
#include <Timer.h>
#include <AnalyticGeometry.h>
#include <DiscreteGeometry.h>
#include <HomMix/HomMix.h>
#include <HomHydro/HomHydro.h>
#include <HomRT/HomRT.h>
#include <SPINSFlow/SPINSFlow.h>
#include <Quinoa/InputDeck/Parser.h>

using namespace quinoa;

QuinoaDriver::QuinoaDriver(int argc, char** argv, Base& base)
  : m_base(base)
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \param[in] base      Essentials
//! \details   Instantiate geometry, physics, set initial conditions, etc.
//! \author J. Bakosi
//******************************************************************************
{
  // Take exactly one filename argument for now
  // Will need to be extended with a more elaborate command line parser
  ErrChk(argc == 2, ExceptType::FATAL,
         "Exactly one command line argument required: filename.q");

  // Instantiate control file parser
  InputDeckParser parser(argv[1], base);

  // Parse control file
  parser.parse();

  base.print.part("Problem setup");
  base.print.section("Title", base.control.get<ctr::title>());

  //! Initialize factory
  initFactory();

  // Instantiate geometry object
  sel::GeometryType g = m_base.control.get<ctr::selected, ctr::geometry>();
  if (g != sel::GeometryType::NO_GEOMETRY) {
    m_geometry = std::unique_ptr<Geometry>(m_geometryFactory[g]());
  }

  // Instantiate physics object
  sel::PhysicsType p = m_base.control.get<ctr::selected, ctr::physics>();
  if (p != sel::PhysicsType::NO_PHYSICS) {
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
//  Initialize physics factory
//! \author  J. Bakosi
//******************************************************************************
{
  // Register geometry classes
  m_geometryFactory[sel::GeometryType::ANALYTIC] =
    std::bind(boost::factory<AnalyticGeometry*>(), m_base);
  m_geometryFactory[sel::GeometryType::DISCRETE] =
    std::bind(boost::factory<DiscreteGeometry*>(), m_base);

  // Register physics classes
  m_physicsFactory[sel::PhysicsType::HOMOGENEOUS_MIX] =
    std::bind(boost::factory<HomMix*>(), m_base);
  m_physicsFactory[sel::PhysicsType::HOMOGENEOUS_HYDRO] =
    std::bind(boost::factory<HomHydro*>(), m_base);
  m_physicsFactory[sel::PhysicsType::HOMOGENEOUS_RAYLEIGH_TAYLOR] =
    std::bind(boost::factory<HomRT*>(), m_base);
  m_physicsFactory[sel::PhysicsType::SPINSFLOW] =
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
