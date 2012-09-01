//******************************************************************************
/*!
  \file      src/Base/Memory.h
  \author    J. Bakosi
  \date      Sat 01 Sep 2012 03:20:20 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory store declaration
  \details   Memory store base class declaration
*/
//******************************************************************************
#ifndef Memory_h
#define Memory_h

#include <map>
#include <string>
#include <cassert>
#include <iostream>

#include <MemoryEntry.h>
#include <Macros.h>

namespace Quinoa {

//! Memory store base class
template<class V>
class Memory {

  typedef map<Int, MemoryEntry<V>*> EntryMap;

  public:
    //! Constructor
    Memory() {
      m_entries = 0;
    }

    //! Destructor: Free all allocated memory when leaving scope
    ~Memory() {
      freeAll();
    }

    //! Allocate memory entry
    Int newEntry(size_t number,
                 VarType varType,
                 string name,
                 bool plot = false,
                 bool restart = false) {
      assert(number > 0);
      assert(0 <= varType && varType < NUM_VAR_TYPES);
      assert(name.size() > 0);
      V* ptr;
      m_entry[m_entries] = new MemoryEntry<V>(number,
                                              varType,
                                              name,
                                              plot,
                                              restart,
                                              ptr = new V [number*VarComp[varType]]);
      cout << name << ": " << ptr << endl;
      return m_entries++;
    }

    //! Deallocate all memory entries
    void freeAll() {
      if (m_entries) {
        do delete m_entry[--m_entries]; while (m_entries);
        m_entry.clear();
      }
    }

    //! Return raw pointer for memory entry
    V* getPtr(Int id) {
      typename EntryMap::const_iterator it = m_entry.find(id);
      if (it!=m_entry.end()) {
        return it->second->m_ptr;
      } else {
        ERR("Variable not found!");
      }
    }

    //! Return number of allocated bytes
    size_t getBytes() {
      if (m_entries) {
        size_t bytes = 0;
        typename EntryMap::const_iterator it;
        for (it=m_entry.begin(); it!=m_entry.end(); it++) {
          MemoryEntry<V>* e = it->second;
          bytes += 3*sizeof(MemoryEntry<V>) +
                   e->m_number * VarComp[e->m_type] * sizeof(V);
        }
        return bytes;
      } else {
        return 0;
      }
    }

  private:
    //! Don't permit copy operator
    Memory(const Memory&);
    //! Don't permit assigment operator
    Memory& operator=(const Memory&);

    //! Current number of entries in the memory store
    Int m_entries;

    //! Memory entries
    EntryMap m_entry;
};

} // namespace Quinoa

#endif // Memory_h
