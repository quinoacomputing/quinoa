//******************************************************************************
/*!
  \file      src/Base/Memory.h
  \author    J. Bakosi
  \date      Sun 02 Sep 2012 08:26:49 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory (a store for MemoryEntry objects) base class declaration
  \details   Memory (a store for MemoryEntry objects) base class declaration
*/
//******************************************************************************
#ifndef Memory_h
#define Memory_h

#include <string>
#include <unordered_set>
#include <unordered_map>

#include <MemoryEntry.h>

namespace Quinoa {

//! Memory store
class Memory {

  //! Memory entries are stored in an STL unordered_set.
  //! The keys of the set are pointers to (dynamically allocated) MemoryEntry
  //! class instances that store basic metadata on the memory entry allocated.
  //! Compared to O(log n) in standard sets, the cost of searches, insertions,
  //! and deletions in unordered sets (i.e. retrieving raw pointers, allocating,
  //! and deallocating memory entries) is amortized to O(1).
  typedef unordered_set<MemoryEntry*> MemorySet;

  //! Map of memory entry names to MemorySet keys
  typedef unordered_map<string,MemoryEntry*> MemoryNames;

  public:
    //! Constructor
    Memory();

    //! Destructor
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
      if (id == 0) return 0;
      auto it = m_entry.find(id);
      if (it!=m_entry.end())
        return static_cast<V*>((*it)->m_ptr);
      else
        return 0;
    }

    //! Return the MemorySet key based on the variable name
    MemoryEntry* getID(string name);

    //! Return the number of allocated bytes
    size_t getBytes();

  private:
    //! Don't permit copy operator
    Memory(const Memory&);
    //! Don't permit assigment operator
    Memory& operator=(const Memory&);

    //! Memory entries
    MemorySet m_entry;

    //! Memory entry names mapped to MemorySet keys
    MemoryNames m_name;
};

} // namespace Quinoa

#endif // Memory_h
