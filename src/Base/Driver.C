//******************************************************************************
/*!
  \file      src/Base/Driver.C
  \author    J. Bakosi
  \date      Mon 10 Sep 2012 03:57:07 AM KST
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
