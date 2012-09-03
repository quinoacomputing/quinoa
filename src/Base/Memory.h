//******************************************************************************
/*!
  \file      src/Base/Memory.h
  \author    J. Bakosi
  \date      Sun 02 Sep 2012 06:47:13 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory (a store for MemoryEntry objects) base class declaration
  \details   Memory (a store for MemoryEntry objects) base class declaration
*/
//******************************************************************************
#ifndef Memory_h
#define Memory_h

#include <string>
#include <unordered_set>
#include <cassert>
#include <iostream>

#include <MemoryEntry.h>

namespace Quinoa {

//! Memory store
class Memory {

  //! Memory entries are stored in an STL unordered_set
  //! Compared to O(log n) in standard sets, the cost of searches, insertions,
  //! and deletions in unordered sets are amortized to O(1).
  typedef unordered_set<MemoryEntry*> MemorySet;

  public:
    //! Constructor
    Memory();

    //! Destructor: Free all allocated memory when leaving scope
    ~Memory();

    //! Allocate memory entry
    MemoryEntry* newEntry(size_t number,
                          ValueType value,
                          VariableType variable,
                          string name,
                          Bool plot = false,
                          Bool restart = false);

    //! Deallocate a memory entry
    void freeEntry(MemoryEntry* id);

    //! Deallocate all memory entries
    void freeAllEntries();

    //! Return raw pointer for memory entry,
    //! template V specifies return pointer type
    template<class V> V* getPtr(MemoryEntry* id) {
      auto it = m_entry.find(id);
      if (it!=m_entry.end())
        return static_cast<V*>((*it)->m_ptr);
      else
        return 0;
    }

    //! Return number of allocated bytes
    size_t getBytes();

  private:
    //! Don't permit copy operator
    Memory(const Memory&);
    //! Don't permit assigment operator
    Memory& operator=(const Memory&);

    //! Memory entries
    MemorySet m_entry;
};

} // namespace Quinoa

#endif // Memory_h
