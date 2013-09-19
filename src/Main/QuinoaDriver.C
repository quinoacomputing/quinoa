//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.C
  \author    J. Bakosi
  \date      Thu Sep 19 09:11:30 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************

#include <iostream>

#include <QuinoaDriver.h>
#include <Timer.h>
#include <AnalyticGeometry.h>
#include <DiscreteGeometry.h>
#include <HomMix.h>
#include <HomHydro.h>
#include <HomRT.h>
#include <SPINSFlow.h>
#include <QuinoaParser.h>

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
  // Note: base is const, but allow parser to populate by removing const
  QuinoaParser parser(argv[1], base.print, base.control);

  // Parse control file
  parser.parse();

  base.print.part("Problem setup");

  // Echo information of stuff parsed
  parser.echo();

  // Instantiate geometry object
  initGeometry();

  // Instantiate physics object
  initPhysics();
}

void
QuinoaDriver::initGeometry()
//******************************************************************************
//  Instantiate geometry object (if any)
//! \author J. Bakosi
//******************************************************************************
{
  using namespace ctr;
  using namespace select;

  //  Instantiate geometry object (if any)
  if (m_base.control.get<selected,geometry>() == GeometryType::ANALYTIC) {
    m_geometry = std::unique_ptr<AnalyticGeometry>(new AnalyticGeometry(m_base));
  } else if (m_base.control.get<selected,geometry>() == GeometryType::DISCRETE) {
    m_geometry = std::unique_ptr<DiscreteGeometry>(new DiscreteGeometry(m_base));
  }

  // Initialize geometry (if any)
  if (m_geometry) m_geometry->init();
}

void
QuinoaDriver::initPhysics()
//******************************************************************************
//  Instantiate physics object (if any)
//! \author J. Bakosi
//******************************************************************************
{
  using namespace ctr;
  using namespace select;

  //  Instantiate physics object (if any)
  if (m_base.control.get<selected,physics>() == PhysicsType::HOMOGENEOUS_MIX) {
    m_physics = std::unique_ptr<HomMix>(new HomMix(m_base));
  } else if (m_base.control.get<selected,physics>() ==
             PhysicsType::HOMOGENEOUS_HYDRO) {
    m_physics = std::unique_ptr<HomHydro>(new HomHydro(m_base));
  } else if (m_base.control.get<selected,physics>() ==
             PhysicsType::HOMOGENEOUS_RAYLEIGH_TAYLOR) {
    m_physics = std::unique_ptr<HomRT>(new HomRT(m_base));
  } if (m_base.control.get<selected,physics>() == PhysicsType::SPINSFLOW) {
    m_physics = std::unique_ptr<SPINSFlow>(new SPINSFlow(m_base));
  }

  // Initialize physics (if any)
  if (m_physics) m_physics->init();
}

void
QuinoaDriver::execute() const
//******************************************************************************
//  Execute
//! \author J. Bakosi
//******************************************************************************
{
  //if (m_geometry) m_geometry->fill();
  if (m_physics) m_physics->solve();
}
