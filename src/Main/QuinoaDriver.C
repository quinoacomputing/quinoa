//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.C
  \author    J. Bakosi
  \date      Sun 15 Sep 2013 10:18:45 PM MDT
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

QuinoaDriver::~QuinoaDriver() noexcept
//******************************************************************************
//  Destructor
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  if (m_geometry) { delete m_geometry; m_geometry = nullptr; }
  if (m_physics)  { delete m_physics;  m_physics  = nullptr; }
}

void
QuinoaDriver::initGeometry()
//******************************************************************************
//  Instantiate geometry object (if any)
//! \author J. Bakosi
//******************************************************************************
{
  using namespace control;
  using namespace select;

  //  Instantiate geometry object (if any)
  if (m_base.control.get<selected,geometry>() == GeometryType::ANALYTIC) {

    m_geometry = new(std::nothrow) AnalyticGeometry(m_base);
    ErrChk(m_geometry != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for geometry object");

  } else if (m_base.control.get<selected,geometry>() == GeometryType::DISCRETE) {

    m_geometry = new(std::nothrow) DiscreteGeometry(m_base);
    ErrChk(m_geometry != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for geometry object");

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
  using namespace control;
  using namespace select;

  //  Instantiate physics object (if any)
  if (m_base.control.get<selected,physics>() == PhysicsType::HOMOGENEOUS_MIX) {

    m_physics = new(std::nothrow) HomMix(m_base);
    ErrChk(m_physics != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for physics object");

  } else if (m_base.control.get<selected,physics>() ==
             PhysicsType::HOMOGENEOUS_HYDRO) {

    m_physics = new(std::nothrow) HomHydro(m_base);
    ErrChk(m_physics != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for physics object");

  } else if (m_base.control.get<selected,physics>() ==
             PhysicsType::HOMOGENEOUS_RAYLEIGH_TAYLOR) {

    m_physics = new(std::nothrow) HomRT(m_base);
    ErrChk(m_physics != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for physics object");

  } if (m_base.control.get<selected,physics>() == PhysicsType::SPINSFLOW) {

    m_physics = new(std::nothrow) SPINSFlow(m_base);
    ErrChk(m_physics != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for physics object");

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
  if (m_geometry) m_geometry->fill();
  if (m_physics) m_physics->solve();
}
