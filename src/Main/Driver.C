//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Thu Nov 15 16:19:20 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class definition
  \details   Driver base class definition
*/
//******************************************************************************

#include <Driver.h>
#include <Setup.h>
#include <ModelException.h>
#include <HomDirichlet.h>
#include <HomGenDirichlet.h>

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
  m_model = nullptr;
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
//  Setup: instantiate model, set initial conditions
//! \author J. Bakosi
//******************************************************************************
{
  // Instantiate selected model
  // ICC: this could be a switch
  if (MODEL_TYPE == ModelType::HOMOGENEOUS_DIRICHLET) {
    m_model = new (nothrow) HomDirichlet(m_memory, m_paradigm, NSCALAR, NPAR);
    Assert(m_model != nullptr, MemoryException,FATAL,BAD_ALLOC);
  }
  else if (MODEL_TYPE == ModelType::HOMOGENEOUS_GENDIRICHLET) {
    m_model = new (nothrow) HomGenDirichlet(m_memory, m_paradigm, NSCALAR);
    Assert(m_model != nullptr, MemoryException,FATAL,BAD_ALLOC);
  } else {
    Throw(ModelException,FATAL,NO_SUCH_MODEL);
  }

  // Echo information on model selected
  m_model->echo();

  // Set initial conditions
  m_model->init();
}

void
Driver::solve()
//******************************************************************************
//  Solve
//! \author J. Bakosi
//******************************************************************************
{
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
  if (m_model) { delete m_model; m_model = nullptr; }
  m_memory->freeAllEntries();
}
