//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 06:36:29 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class definition
  \details   Driver base class definition
*/
//******************************************************************************

#include <Driver.h>
#include <Model.h>
#include <Setup.h>

using namespace Quinoa;

Driver::Driver(Memory* memory) : m_memory(memory)
//******************************************************************************
//  Constructor
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
  // Instantiate model
  m_model = new (nothrow) Model(MODEL_TYPE, NPEL, m_memory);
  Assert(m_model != nullptr, MemoryException,FATAL,BAD_ALLOC);

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
