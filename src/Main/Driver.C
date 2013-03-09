//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Fri Mar  8 15:27:39 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class definition
  \details   Driver base class definition
*/
//******************************************************************************

#include <Driver.h>
#include <Control.h>
#include <Timer.h>
#include <Parser.h>
#include <ParserException.h>
#include <PhysicsException.h>
#include <HomMix.h>
#include <SPINSFlow.h>

using namespace Quinoa;

Driver::Driver(int argc,
               char** argv,
               Memory* const memory, Paradigm* const paradigm) :
  m_memory(memory), m_paradigm(paradigm), m_argc(argc), m_argv(argv)
//******************************************************************************
//  Constructor
//! \param[in]  argc      Argument count from command line
//! \param[in]  argv      Argument vector from command line
//! \param[in]  memory    Memory oject pointer
//! \param[in]  paradigm  Parallel programming paradigm object pointer
//! \author J. Bakosi
//******************************************************************************
{
  m_control = nullptr;
  m_timer   = nullptr;
  m_physics = nullptr;
}

Driver::~Driver()
//******************************************************************************
//  Destructor
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
  Assert(m_control != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Take exactly one filename argument for now
  // Will need to be extended with a more elaborate command line parser
  if (m_argc != 2) Throw(ParserException,FATAL,CMDLINE_EXCEPT);

  // Instantiate control file parser
  Parser parser(m_argv[1], m_control);

  // Parse control file
  parser.parse();

  // Echo information of stuff parsed
  parser.echo();

  // Instantiate timer object
  m_timer = new (nothrow) Timer;
  Assert(m_timer != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Instantiate selected physics
  switch (m_control->get<control::PHYSICS>()) {

    case control::PhysicsType::NO_PHYSICS :
      Throw(PhysicsException,FATAL,PhysicsExceptType::NO_PHYSICS);
      break;

    case control::PhysicsType::HOMOGENEOUS_MIX :
      m_physics = new (nothrow) HomMix(m_memory, m_paradigm, m_control, m_timer);
      Assert(m_physics != nullptr, MemoryException,FATAL,BAD_ALLOC);
      break;

    case control::PhysicsType::SPINSFLOW :
      m_physics = new (nothrow)
                  SPINSFlow(m_memory, m_paradigm, m_control, m_timer,
                            "cylinder.msh");
      Assert(m_physics != nullptr, MemoryException,FATAL,BAD_ALLOC);
      break;

    default :
      Throw(PhysicsException,FATAL,PHYSICS_UNIMPLEMENTED);
  }

  // Echo information on physics selected
  m_physics->echo();

  // Set initial conditions
  m_physics->init();
}

void
Driver::solve()
//******************************************************************************
//  Solve
//! \author J. Bakosi
//******************************************************************************
{
  m_physics->solve();
}

void
Driver::finalize()
//******************************************************************************
//  Finalize
//! \details Cleanup either at the end of business as usual or due to an
//!          exception
//! \author J. Bakosi
//******************************************************************************
{
  if (m_physics) { delete m_physics; m_physics = nullptr; }
  if (m_timer)   { delete m_timer;   m_timer   = nullptr; }
  if (m_control) { delete m_control; m_control = nullptr; }
  m_memory->freeAllEntries();
}
