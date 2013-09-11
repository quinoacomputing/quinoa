//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.C
  \author    J. Bakosi
  \date      Wed Sep 11 16:41:33 2013
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

QuinoaDriver::QuinoaDriver(int argc,
                           char** argv,
                           Memory* const memory,
                           Paradigm* const paradigm,
                           const QuinoaPrinter& print)
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \param[in] memory    Memory oject pointer
//! \param[in] paradigm  Parallel programming paradigm object pointer
//! \param[in] print     Quinoa's pretty printer
//! \details   Instantiate geometry, physics, set initial conditions, etc.
//! \author J. Bakosi
//******************************************************************************
try :
  Driver(),
  m_memory(memory),
  m_paradigm(paradigm),
  m_print(print),
  m_geometry(nullptr),
  m_physics(nullptr)
{

  // Take exactly one filename argument for now
  // Will need to be extended with a more elaborate command line parser
  ErrChk(argc == 2, ExceptType::FATAL,
         "Exactly one command line argument required: filename.q");

  // Instantiate control file parser
  QuinoaParser parser(argv[1], m_print, m_control);

  // Parse control file
  parser.parse();

  // Echo information of stuff parsed
  parser.echo();

  // Instantiate geometry object
  initGeometry();

  // Instantiate physics object
  initPhysics();

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
  if (m_control.get<selected,geometry>() == GeometryType::ANALYTIC) {

    m_geometry = new(std::nothrow)
                   AnalyticGeometry(m_memory, m_paradigm, m_control, m_timer);
    ErrChk(m_geometry != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for geometry object");

  } else if (m_control.get<selected,geometry>() == GeometryType::DISCRETE) {

    m_geometry = new(std::nothrow)
                  DiscreteGeometry(m_memory, m_paradigm, m_control, m_timer);
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
  if (m_control.get<selected,physics>() == PhysicsType::HOMOGENEOUS_MIX) {

    m_physics = new(std::nothrow)
                HomMix(m_memory, m_paradigm, m_control, m_timer, m_print);
    ErrChk(m_physics != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for physics object");

  } else if (m_control.get<selected,physics>() ==
             PhysicsType::HOMOGENEOUS_HYDRO) {

    m_physics = new(std::nothrow)
                HomHydro(m_memory, m_paradigm, m_control, m_timer, m_print);
    ErrChk(m_physics != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for physics object");

  } else if (m_control.get<selected,physics>() ==
             PhysicsType::HOMOGENEOUS_RAYLEIGH_TAYLOR) {

    m_physics = new(std::nothrow)
                HomRT(m_memory, m_paradigm, m_control, m_timer, m_print);
    ErrChk(m_physics != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for physics object");

  } if (m_control.get<selected,physics>() == PhysicsType::SPINSFLOW) {

    m_physics = new(std::nothrow)
                  SPINSFlow(m_memory, m_paradigm, m_control, m_timer, m_print,
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
