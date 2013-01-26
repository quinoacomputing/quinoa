//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Sat 26 Jan 2013 09:46:59 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class definition
  \details   Driver base class definition
*/
//******************************************************************************

#include <Driver.h>
#include <Setup.h>
#include <Control.h>
#include <Parser.h>
#include <ParserException.h>
#include <PhysicsException.h>
#include <HomogeneousDirichlet.h>
#include <HomogeneousGeneralizedDirichlet.h>
#include <SPINSFlow.h>

using namespace Quinoa;

Driver::Driver(int argc, char** argv, Memory* memory, Paradigm* paradigm) :
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
  parser.Parse();

  // Instantiate selected physics
  // ICC: this could be a switch once ICC supports scoped enums in switch
  if (PHYSICS_TYPE == PhysicsType::HOMOGENEOUS_DIRICHLET) {
    m_physics = new (nothrow) HomogeneousDirichlet(m_memory,
                                                   m_paradigm,
                                                   NSCALAR,
                                                   NPAR,
                                                   TIME,
                                                   ECHO,
                                                   NSTEP);
    Assert(m_physics != nullptr, MemoryException,FATAL,BAD_ALLOC);
  } else if (PHYSICS_TYPE == PhysicsType::HOMOGENEOUS_GENDIRICHLET) {
    m_physics = new (nothrow) HomogeneousGeneralizedDirichlet(m_memory,
                                                              m_paradigm,
                                                              NSCALAR,
                                                              TIME,
                                                              NSTEP);
    Assert(m_physics != nullptr, MemoryException,FATAL,BAD_ALLOC);
  } else if (PHYSICS_TYPE == PhysicsType::SPINSFLOW) {
    m_physics = new (nothrow) SPINSFlow(m_memory,
                                        m_paradigm,
                                        HYDRO_TYPE,
                                        NPAR,
                                        MESH_FILENAME,
                                        TIME,
                                        ECHO,
                                        NSTEP);
    Assert(m_physics != nullptr, MemoryException,FATAL,BAD_ALLOC); 
  } else {
    Throw(PhysicsException,FATAL,NO_SUCH_PHYSICS);
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
  if (m_control) { delete m_control; m_control = nullptr; }
  if (m_physics) { delete m_physics; m_physics = nullptr; }
  m_memory->freeAllEntries();
}
