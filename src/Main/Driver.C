//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Fri Jul 26 13:51:19 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class definition
  \details   Driver base class definition
*/
//******************************************************************************

#include <iostream>

#include <Driver.h>
#include <Control.h>
#include <Timer.h>
#include <Parser.h>
#include <AnalyticGeometry.h>
#include <DiscreteGeometry.h>
#include <HomMix.h>
#include <HomHydro.h>
#include <HomRT.h>
#include <SPINSFlow.h>

using namespace Quinoa;

Driver::Driver(int argc,
               char** argv,
               Memory* const memory,
               Paradigm* const paradigm)
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \param[in] memory    Memory oject pointer
//! \param[in] paradigm  Parallel programming paradigm object pointer
//! \details   Instantiate geometry, physics, set initial conditions, etc.
//! \author J. Bakosi
//******************************************************************************
try :
  m_memory(memory),
  m_paradigm(paradigm),
  m_geometry(nullptr),
  m_physics(nullptr),
  m_control(nullptr),
  m_timer(nullptr)
{

  // Instantiate main control category
  m_control = new(std::nothrow) Control;
  ErrChk(m_control != nullptr, ExceptType::FATAL,
         "Cannot allocate memory for control object");

  // Take exactly one filename argument for now
  // Will need to be extended with a more elaborate command line parser
  ErrChk(argc == 2, ExceptType::FATAL,
         "Exactly one command line argument required: filename.q");

  // Instantiate control file parser
  Parser parser(argv[1], m_control);

  // Parse control file
  parser.parse();

  // Echo information of stuff parsed
  parser.echo();

  // Instantiate timer object
  m_timer = new(std::nothrow) Timer;
  ErrChk(m_timer != nullptr, ExceptType::FATAL,
         "Cannot allocate memory for timer object");

  // Instantiate physics object
  initPhysics();

  // Instantiate geometry object
  initGeometry();

} // Roll back changes and rethrow on error
  catch (std::exception&) {
    finalize();
    throw;
  }
  // Catch uncaught exceptions
  catch (...) {
    finalize();
    Throw(ExceptType::UNCAUGHT, "Non-standard exception");
  }

Driver::~Driver() noexcept
//******************************************************************************
//  Destructor
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  finalize();
}

void
Driver::finalize() noexcept
//******************************************************************************
//  Finalize
//! \details Cleanup either at the end of business as usual or due to an
//!          exception. No-throw guarantee: this member function never throws
//!          exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  if (m_geometry) { delete m_geometry; m_geometry = nullptr; }
  if (m_physics)  { delete m_physics;  m_physics  = nullptr; }
  if (m_timer)    { delete m_timer;    m_timer    = nullptr; }
  if (m_control)  { delete m_control;  m_control  = nullptr; }
  m_memory->freeAllEntries();
}

void
Driver::initGeometry()
//******************************************************************************
//  Instantiate geometry object (if any)
//! \author J. Bakosi
//******************************************************************************
{
  //  Instantiate geometry object (if any)
  if (m_control->get<control::GEOMETRY>() == select::GeometryTypes::ANALYTIC) {

    m_geometry = new(std::nothrow)
                   AnalyticGeometry(m_memory, m_paradigm, m_control, m_timer);

  } else if (m_control->get<control::GEOMETRY>() ==
               select::GeometryTypes::DISCRETE) {

    m_geometry = new(std::nothrow)
                  DiscreteGeometry(m_memory, m_paradigm, m_control, m_timer);

  }

  // Initialize geometry (if any)
  if (m_geometry) m_geometry->init();
}

void
Driver::initPhysics()
//******************************************************************************
//  Instantiate physics object (if any)
//! \author J. Bakosi
//******************************************************************************
{
  //  Instantiate physics object (if any)
  if (m_control->get<control::PHYSICS>() ==
        select::PhysicsTypes::HOMOGENEOUS_MIX) {

    m_physics =
      new(std::nothrow) HomMix(m_memory, m_paradigm, m_control, m_timer);

  } else if (m_control->get<control::PHYSICS>() ==
               select::PhysicsTypes::HOMOGENEOUS_HYDRO) {

    m_physics =
      new(std::nothrow) HomHydro(m_memory, m_paradigm, m_control, m_timer);

  } else if (m_control->get<control::PHYSICS>() ==
               select::PhysicsTypes::HOMOGENEOUS_RAYLEIGH_TAYLOR) {

    m_physics =
      new(std::nothrow) HomRT(m_memory, m_paradigm, m_control, m_timer);

  } if (m_control->get<control::PHYSICS>() == select::PhysicsTypes::SPINSFLOW) {

    m_physics = new(std::nothrow)
                  SPINSFlow(m_memory, m_paradigm, m_control, m_timer,
                  "cylinder.msh");

  }

  // Initialize physics (if any)
  if (m_physics) m_physics->init();
}

void
Driver::execute() const
//******************************************************************************
//  Execute
//! \author J. Bakosi
//******************************************************************************
{
  if (m_geometry) m_geometry->fill();
  if (m_physics) m_physics->solve();
}
