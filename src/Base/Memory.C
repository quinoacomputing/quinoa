//******************************************************************************
/*!
  \file      src/Base/Memory.C
  \author    J. Bakosi
  \date      Thu 30 Aug 2012 11:29:05 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory entry definition
  \details   Memory entry definition
*/
//******************************************************************************

#include <string>
#include <cassert>

using namespace std;

#include "MemoryEntry.h"
#include "Memory.h"

using namespace Quinoa;

Memory::Memory()
//******************************************************************************
//  Constructor
//! \author    J. Bakosi
//******************************************************************************
{
}

Memory::~Memory()
//******************************************************************************
//  Destructor
//! \author    J. Bakosi
//******************************************************************************
{
}

void
Memory::allocateEntry(size_t number,
                      VariableType type,
                      int dim,
                      string name,
                      StorageLayout layout,
                      bool plot,
                      bool restart
)
//******************************************************************************
//  Allocate memory entry
//! \details
//! Allocate memory entry detailed description ...
//!
//! \param[in]  number  Number of items per dimension
//! \param[in]  type    Variable type
//! \param[in]  dim     Number of dimensions
//! \param[in]  name    Variable name
//! \param[in]  layout  Memory storage layout
//! \param[in]  plot    Variable can be plotted
//! \param[in]  restart Write to restart file
//!
//! \author    J. Bakosi
//******************************************************************************
{
  assert(number > 0);
  assert(dim > 0);
  assert(string.size() > 0);

  // Compute size of memory to allocate
  size_t nbytes = dim*number;//*VariableSize[type];

  // Allocate memory
  //ptr = (void*) new char [nbytes];
}
