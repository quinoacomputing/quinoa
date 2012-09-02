//******************************************************************************
/*!
  \file      src/Base/Memory.h
  \author    J. Bakosi
  \date      Sun 02 Sep 2012 02:49:00 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory (a store for MemoryEntry objects) base class declaration
  \details   Memory (a store for MemoryEntry objects) base class declaration
*/
//******************************************************************************
#ifndef Memory_h
#define Memory_h

#include <vector>
#include <string>
#include <cassert>
#include <iostream>

#include <MemoryEntry.h>

namespace Quinoa {

//! Templated memory store for variable type V
class Memory {

  public:
    //! Constructor
    Memory();

    //! Destructor: Free all allocated memory when leaving scope
    ~Memory();

    //! Allocate memory entry
    Int newEntry(size_t number,
                 ValueType value,
                 VariableType variable,
                 string name,
                 Bool plot = false,
                 Bool restart = false);

    //! Deallocate all memory entries
    void freeAll();

    //! Return raw pointer for memory entry,
    //! template V specifies return pointer type
    template<class V> V* getPtr(Int id) {
      return static_cast<V*>(m_entry[id]->m_ptr);
    }

    //! Return number of allocated bytes
    size_t getBytes();

  private:
    //! Don't permit copy operator
    Memory(const Memory&);
    //! Don't permit assigment operator
    Memory& operator=(const Memory&);

    //! Current number of entries
    Int m_entries;

    //! Memory entries
    vector<MemoryEntry*> m_entry;
};

} // namespace Quinoa

#endif // Memory_h
