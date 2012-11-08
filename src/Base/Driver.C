//******************************************************************************
/*!
  \file      src/Base/Driver.C
  \author    J. Bakosi
  \date      Wed 07 Nov 2012 09:00:40 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class definition
  \details   Driver base class definition
*/
//******************************************************************************

#include <Driver.h>

using namespace Quinoa;

Driver::~Driver()
//******************************************************************************
//  Destructor
//! \author J. Bakosi
//******************************************************************************
{
  finalize();
}

void
Driver::setup(int argc, char* argv[])
//******************************************************************************
//  Setup
//! \author J. Bakosi
//******************************************************************************
{
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
  m_memory->freeAllEntries();
}
