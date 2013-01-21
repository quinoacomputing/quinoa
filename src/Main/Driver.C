//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Sun 20 Jan 2013 05:53:27 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class definition
  \details   Driver base class definition
*/
//******************************************************************************

#include <Driver.h>
#include <Setup.h>
#include <PhysicsException.h>
#include <HomogeneousDirichlet.h>
#include <HomogeneousGeneralizedDirichlet.h>
#include <SPINSFlow.h>

using namespace Quinoa;

Driver::Driver(Memory* memory, Paradigm* paradigm) :
  m_memory(memory), m_paradigm(paradigm)
//******************************************************************************
//  Constructor
//! \param[in]  memory    Memory oject pointer
//! \param[in]  paradigm  Parallel programming paradigm object pointer
//! \author J. Bakosi
//******************************************************************************
{
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
                                        TIME,
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
  if (m_physics) { delete m_physics; m_physics = nullptr; }
  m_memory->freeAllEntries();
}
