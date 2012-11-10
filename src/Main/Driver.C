//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Sat 10 Nov 2012 08:23:22 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class definition
  \details   Driver base class definition
*/
//******************************************************************************

#include <Driver.h>
#include <Setup.h>
#include <Dirichlet.h>
#include <GeneralizedDirichlet.h>

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
//  Setup
//! \author J. Bakosi
//******************************************************************************
{
  // Select model
  // ICC: this could be a switch
  if (MODEL == ModelType::DIRICHLET) {
    m_model = new Dirichlet(NUM_SCALARS);
  } else if (MODEL == ModelType::GENERALIZED_DIRICHLET) {
    m_model = new GeneralizedDirichlet(NUM_SCALARS);
  } else {
    throw Exception(FATAL, "No such model");
  }
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
