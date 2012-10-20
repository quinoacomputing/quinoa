//******************************************************************************
/*!
  \file      src/Base/Memory.h
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 04:41:21 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory store, container of memory entries
  \details   Memory store, container of memory entries
*/
//******************************************************************************
#ifndef Memory_h
#define Memory_h

#include <string>
#include <unordered_set>
#include <unordered_map>

#include <MemoryEntry.h>
#include <MemoryException.h>

namespace Quinoa {

//! Memory store, container of memory entries
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
    Memory(int nthreads) : m_nthreads(nthreads) {}

    //! Destructor
    ~Memory();

    //! Allocate memory entry
    MemoryEntry* newEntry(size_t number,
                          ValType value,
                          VarType variable,
                          string name,
                          bool plot = false,
                          bool restart = false);

    //! Allocate and zero memory entry
    MemoryEntry* newZeroEntry(size_t number,
                              ValType value,
                              VarType variable,
                              string name,
                              bool plot = false,
                              bool restart = false);

    //! Deallocate a memory entry
    void freeEntry(MemoryEntry* id);

    //! Deallocate all memory entries
    void freeAllEntries() noexcept;

    //! Echo all (optionally sorted) memory entries
    void echoAllEntries(MemoryEntryField crit = UNSPECIFIED);

    //! Return the number of items based on the ID
    size_t getNumber(MemoryEntry* id);

    //! Return the number of items based on the variable name
    size_t getNumber(string name);

    //! Return the value type based on the ID
    ValType getValue(MemoryEntry* id);

    //! Return the variable type based on the ID
    VarType getVariable(MemoryEntry* id);

    //! Return the variable name based on the ID
    string getName(MemoryEntry* id);

    //! Return true if the variable can be plotted based on the ID
    bool getPlot(MemoryEntry* id);

    //! Return true if the variable is written to restart file based on the ID
    bool getRestart(MemoryEntry* id);

    //! Return data pointer for memory entry based on ID,
    //! template V specifies return pointer type
    template<class V> V* getPtr(MemoryEntry* id) {
      if (id == nullptr)
        throw MemoryException(WARNING, UNDEFINED);
      auto it = m_entry.find(id);
      if (it==m_entry.end())
        throw MemoryException(WARNING, NOT_FOUND);
      return static_cast<V*>((*it)->m_ptr);
    }

    //! Return the data pointer for memory entry based on the variable name,
    //! template V specifies return pointer type
    template<class V> V* getPtr(string name) {
      return static_cast<V*>(getID(name)->m_ptr);
    }

    //! Return the pair of data pointer and number of variables for memory entry
    //! based on the variable name, template V specifies return pointer type
    template<class V> pair<size_t,V*> getNumPtr(string name) {
      MemoryEntry* entry = getID(name);
      return pair<size_t,V*>(entry->m_number, static_cast<V*>(entry->m_ptr));
    }

    //! Return the MemorySet key based on the variable name
    MemoryEntry* getID(string name);

    //! Return the number of allocated bytes
    size_t getBytes();

    //! Zero entry using multiple threads
    void zero(MemoryEntry* id);

  private:
    //! Don't permit copy constructor
    Memory(const Memory&) = delete;
    //! Don't permit copy assigment
    Memory& operator=(const Memory&) = delete;
    //! Don't permit move constructor
    Memory(Memory&&) = delete;
    //! Don't permit move assigment
    Memory& operator=(Memory&&) = delete;

    void echo();           //!< Echo unsorted entries
    void echoByBytes();    //!< Echo entries sorted by Bytes
    void echoByNumber();   //!< Echo entries sorted by Number
    void echoByValue();    //!< Echo entries sorted by Value
    void echoByVariable(); //!< Echo entries sorted by Variable
    void echoByName();     //!< Echo entries sorted by Name
    void echoByPlot();     //!< Echo entries sorted by Plot
    void echoByRestart();  //!< Echo entries sorted by Restart

    int m_nthreads;        //!< Local copy of the number of threads
    MemorySet m_entry;     //!< Memory entries
    MemoryNames m_name;    //!< Memory entry names mapped to MemorySet keys
};

} // namespace Quinoa

#endif // Memory_h
