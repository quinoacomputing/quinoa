//******************************************************************************
/*!
  \file      src/Base/Memory.h
  \author    J. Bakosi
  \date      Thu Aug 29 14:52:48 2013
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
#include <Exception.h>

namespace quinoa {

class Paradigm;

//! Data provides an interface to both Memory Entry ID and its raw pointer
template< class V >
struct Data {
  MemoryEntry* id;
  V* ptr;

  //! Null constructor
  explicit constexpr Data() noexcept : id(nullptr), ptr(nullptr) {}
  //! Data constructor
  explicit constexpr Data(MemoryEntry* const Id, V* const Ptr) noexcept :
    id(Id), ptr(Ptr) {}

  //! Operator [] is overloaded to access the elements behind ptr
  constexpr V&
  operator[] (const size_t index) const noexcept { return ptr[index]; }
  //! Operator + is overloaded to ease pointer arithmetic behind ptr
  constexpr V*
  operator+ (const size_t offset) const noexcept { return ptr + offset; }
};

//! Memory store, container of memory entries
class Memory {

  public:
    //! Constructor
    explicit Memory(Paradigm* const paradigm) noexcept;

    //! Destructor
    ~Memory() noexcept;

    //! Allocate memory entry and return both MemoryEntry* ID and raw pointer,
    //! template V specifies return pointer type
    template<class V>
    Data<V> newEntry(const size_t number,
                     const ValType value,
                     const VarType variable,
                     const std::string name,
                     const bool plot = false,
                     const bool restart = false) {
      MemoryEntry* const me =
        newEntry(number, value, variable, name, plot, restart);
      return Data<V>(me, getPtr<V>(me));
    }

    //! Allocate and zero memory entry and return both MemoryEntry* ID and raw
    //! pointer, template V specifies return pointer type
    template<class V>
    Data<V> newZeroEntry(const size_t number,
                         const ValType value,
                         const VarType variable,
                         const std::string name,
                         const bool plot = false,
                         const bool restart = false) {
      MemoryEntry* const me =
        newZeroEntry(number, value, variable, name, plot, restart);
      return Data<V>(me, getPtr<V>(me));
    }

    //! Deallocate a memory entry
    template< class V >
    void freeEntry(Data<V>& data) noexcept {
      freeEntry(data.id);
      data.ptr = nullptr;
    }

    //! Deallocate all memory entries
    void freeAllEntries() noexcept;

    //! Echo all (optionally sorted) memory entries
    void echoAllEntries(MemoryEntryField crit = UNSPECIFIED) const;

    //! Return the number of items based on the ID
    size_t getNumber(MemoryEntry* const id) const;

    //! Return the number of items based on the variable name
    size_t getNumber(const std::string& name) const;

    //! Return the value type based on the ID
    ValType getValue(MemoryEntry* const id) const;

    //! Return the variable type based on the ID
    VarType getVariable(MemoryEntry* const id) const;

    //! Return the variable name based on the ID
    std::string getName(MemoryEntry* const id) const;

    //! Return true if the variable can be plotted based on the ID
    bool getPlot(MemoryEntry* const id) const;

    //! Return true if the variable is written to restart file based on the ID
    bool getRestart(MemoryEntry* const id) const;

    //! Return the number of allocated bytes
    size_t getBytes() const;

    //! Zero entry using multiple threads
    void zero(MemoryEntry* const id) const;

  private:
    //! Memory entries are stored in an STL unordered_set.
    //! The keys of the set are pointers to (dynamically allocated) MemoryEntry
    //! class instances that store basic metadata on the memory entry allocated.
    //! Compared to O(log n) in standard sets, the cost of searches, insertions,
    //! and deletions in unordered sets (i.e. retrieving raw pointers, allocating,
    //! and deallocating memory entries) is amortized to O(1).
    using MemorySet = std::unordered_set<MemoryEntry*>;

    //! Don't permit copy constructor
    Memory(const Memory&) = delete;
    //! Don't permit copy assigment
    Memory& operator=(const Memory&) = delete;
    //! Don't permit move constructor
    Memory(Memory&&) = delete;
    //! Don't permit move assigment
    Memory& operator=(Memory&&) = delete;

    //! Allocate memory entry
    MemoryEntry* newEntry(const size_t number,
                          const ValType value,
                          const VarType variable,
                          const std::string name,
                          const bool plot = false,
                          const bool restart = false);

    //! Allocate and zero memory entry
    MemoryEntry* newZeroEntry(const size_t number,
                              const ValType value,
                              const VarType variable,
                              const std::string name,
                              const bool plot = false,
                              const bool restart = false);

    //! Return data pointer for memory entry based on ID,
    //! template V specifies return pointer type
    template<class V> V* getPtr(MemoryEntry* id) const {
      Assert(id != nullptr, ExceptType::FATAL,
             "Cannot return a raw pointer for a nullptr MemoryEntry");
      auto it = m_entry.find(id);
      Assert(it!=m_entry.end(), ExceptType::FATAL, "Cannot find memory entry");
      return static_cast<V*>((*it)->m_ptr);
    }

    //! Deallocate a memory entry
    void freeEntry(MemoryEntry* id) noexcept;

    void echo() const;          //!< Echo unsorted entries
    void echoByBytes() const;   //!< Echo entries sorted by Bytes
    void echoByNumber() const;  //!< Echo entries sorted by Number
    void echoByValue() const;   //!< Echo entries sorted by Value
    void echoByVariable() const;//!< Echo entries sorted by Variable
    void echoByName() const;    //!< Echo entries sorted by Name
    void echoByPlot() const;    //!< Echo entries sorted by Plot
    void echoByRestart() const; //!< Echo entries sorted by Restart

    const int m_nOMPthreads;    //!< Number of OpenMP threads

    MemorySet m_entry;          //!< Memory entries
};

} // namespace quinoa

#endif // Memory_h
