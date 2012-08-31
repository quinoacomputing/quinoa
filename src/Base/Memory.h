//******************************************************************************
/*!
  \file      src/Base/Memory.h
  \author    J. Bakosi
  \date      Thu 30 Aug 2012 11:08:36 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory store declaration
  \details   Memory store base class declaration
*/
//******************************************************************************
#ifndef Memory_h
#define Memory_h

#include <set>

#include <MemoryEntry.h>

namespace Quinoa {

//! Memory store base class
class Memory {

  public:
    //! Constructor
             Memory();

    //! Destructor
    virtual ~Memory();

  private:
    //! Don't permit copy operator
    Memory(const Memory&);
    //! Don't permit assigment operator
    Memory& operator=(const Memory&);

    //! Allocate memory entry
    void allocateEntry(size_t number,
                       VariableType type,
                       int dim,
                       string name,
                       StorageLayout layout = LINEAR_ARRAY,
                       bool plot = false,
                       bool restart = false);

    //! Current number of entries in memory store
    int m_entries;

    //! Memory entries
    set<MemoryEntry> m_entry;
};

} // namespace Quinoa

#endif // Memory.h
