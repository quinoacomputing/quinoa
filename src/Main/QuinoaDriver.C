//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.C
  \author    J. Bakosi
  \date      Fri Aug  2 16:31:39 2013
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

using namespace Quinoa;

QuinoaDriver::QuinoaDriver(int argc,
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
  Driver(argc, argv),
  m_memory(memory),
  m_paradigm(paradigm),
  m_geometry(nullptr),
  m_physics(nullptr),
  m_timer(nullptr)
{

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

QuinoaDriver::~QuinoaDriver() noexcept
//******************************************************************************
//  Destructor
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  finalize();
}

void
QuinoaDriver::finalize() noexcept
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
}

void
QuinoaDriver::initGeometry()
//******************************************************************************
//  Instantiate geometry object (if any)
//! \author J. Bakosi
//******************************************************************************
{
  //  Instantiate geometry object (if any)
  if (control()->get<control::GEOMETRY>() == select::GeometryType::ANALYTIC) {

    m_geometry = new(std::nothrow)
                   AnalyticGeometry(m_memory, m_paradigm, control(), m_timer);
    ErrChk(m_geometry != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for geometry object");

  } else if (control()->get<control::GEOMETRY>() ==
               select::GeometryType::DISCRETE) {

    m_geometry = new(std::nothrow)
                  DiscreteGeometry(m_memory, m_paradigm, control(), m_timer);
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
  //  Instantiate physics object (if any)
  if (control()->get<control::PHYSICS>() ==
      select::PhysicsType::HOMOGENEOUS_MIX) {

    m_physics =
      new(std::nothrow) HomMix(m_memory, m_paradigm, control(), m_timer);
    ErrChk(m_physics != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for physics object");

  } else if (control()->get<control::PHYSICS>() ==
             select::PhysicsType::HOMOGENEOUS_HYDRO) {

    m_physics =
      new(std::nothrow) HomHydro(m_memory, m_paradigm, control(), m_timer);
    ErrChk(m_physics != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for physics object");

  } else if (control()->get<control::PHYSICS>() ==
             select::PhysicsType::HOMOGENEOUS_RAYLEIGH_TAYLOR) {

    m_physics =
      new(std::nothrow) HomRT(m_memory, m_paradigm, control(), m_timer);
    ErrChk(m_physics != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for physics object");

  } if (control()->get<control::PHYSICS>() == select::PhysicsType::SPINSFLOW) {

    m_physics = new(std::nothrow)
                  SPINSFlow(m_memory, m_paradigm, control(), m_timer,
                  "cylinder.msh");
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
