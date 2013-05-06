//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Mon May  6 13:28:08 2013
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
#include <HomMix.h>
#include <HomHydro.h>
#include <SPINSFlow.h>

using namespace Quinoa;

Driver::Driver(int argc,
               char** argv,
               Memory* const memory,
               Paradigm* const paradigm) noexcept :
  m_memory(memory),
  m_paradigm(paradigm),
  m_argc(argc),
  m_argv(argv),
  m_physics(nullptr),
  m_control(nullptr),
  m_timer(nullptr)
//******************************************************************************
//  Constructor
//! \param[in]  argc      Argument count from command line
//! \param[in]  argv      Argument vector from command line
//! \param[in]  memory    Memory oject pointer
//! \param[in]  paradigm  Parallel programming paradigm object pointer
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
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
Driver::setup()
//******************************************************************************
//  Setup: instantiate physics, set initial conditions
//! \author J. Bakosi
//******************************************************************************
{
  // Instantiate main control category
  m_control = new (nothrow) Control;
  if (m_control == nullptr)
    throw Exception(FATAL,
            "Cannot allocate memory for control object in Driver::setup()");

  // Take exactly one filename argument for now
  // Will need to be extended with a more elaborate command line parser
  if (m_argc != 2)
    throw Exception(FATAL,
            "Exactly one command line argument required: filename.q");

  // Instantiate control file parser
  Parser parser(m_argv[1], m_control);

  // Parse control file
  parser.parse();

  // Echo information of stuff parsed
  parser.echo();

  // Instantiate timer object
  m_timer = new (nothrow) Timer;
  if (m_timer == nullptr)
    throw Exception(FATAL,
            "Cannot allocate memory for timer object in Driver::setup()");

  // Instantiate selected physics
  switch (m_control->get<control::PHYSICS>()) {

    case control::PhysicsType::NO_PHYSICS :
      throw Exception(FATAL, "No physics selected");
      break;

    case control::PhysicsType::HOMOGENEOUS_MIX :
      m_physics = new (nothrow)
                  HomMix(m_memory, m_paradigm, m_control, m_timer);
      if (m_physics == nullptr)
        throw Exception(FATAL,
                "Cannot allocate memory for physics object in Driver::setup()");
      break;

    case control::PhysicsType::HOMOGENEOUS_HYDRO :
      m_physics = new (nothrow)
                  HomHydro(m_memory, m_paradigm, m_control, m_timer);
      if (m_physics == nullptr)
        throw Exception(FATAL,
                "Cannot allocate memory for physics object in Driver::setup()");
      break;

    case control::PhysicsType::SPINSFLOW :
      m_physics = new (nothrow)
                  SPINSFlow(m_memory, m_paradigm, m_control, m_timer,
                            "cylinder.msh");
      if (m_physics == nullptr)
        throw Exception(FATAL,
                "Cannot allocate memory for physics object in Driver::setup()");
      break;

    default :
      throw Exception(FATAL, "Selected physics not implemented");
  }

  // Echo information on physics selected
  m_physics->echo();

  // Set initial conditions
  m_physics->init();
}

void
Driver::solve() const
//******************************************************************************
//  Solve
//! \author J. Bakosi
//******************************************************************************
{
  m_physics->solve();
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
  if (m_physics) { delete m_physics; m_physics = nullptr; }
  if (m_timer)   { delete m_timer;   m_timer   = nullptr; }
  if (m_control) { delete m_control; m_control = nullptr; }
  m_memory->freeAllEntries();
}
