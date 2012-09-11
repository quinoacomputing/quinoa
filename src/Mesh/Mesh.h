//******************************************************************************
/*!
  \file      src/Mesh/Mesh.h
  \author    J. Bakosi
  \date      Tue 11 Sep 2012 04:58:34 PM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh base class declaration
  \details   Mesh base class declaration
*/
//******************************************************************************
#ifndef Mesh_h
#define Mesh_h

#include <unordered_set>

#include <MemoryEntry.h>
#include <Memory.h>

namespace Quinoa {

//! Mesh dimension
enum MeshDim { TWOD=0,
               THREED
};

//! Mesh base class
class Mesh {

  //! MeshSet stores all memory entry keys that belong to a given instance of
  //! the Mesh object.
  typedef unordered_set<MemoryEntry*> MeshSet;

  public:
    //! Constructor
    Mesh(Memory* memory) : m_memory(memory) {}

    //! Destructor
    virtual ~Mesh();

    //! Add new MeshSet entry (e.g. list of nodes, elements, node ids, etc.)
    template<class V> V* newEntry(size_t number,
                                  ValueType value,
                                  VariableType variable,
                                  string name,
                                  Bool plot = false,
                                  Bool restart = false) {
      // Allocate new mmeory entry
      MemoryEntry* entry = m_memory->newEntry(number,
                                              value,
                                              variable,
                                              name,
                                              plot,
                                              restart);
      // Store new entry
      pair<MeshSet::iterator,Bool> n = m_entry.insert(entry);
      if (!n.second) throw MemoryException(FATAL, BAD_INSERT);
      // Get pointers to element ids and 
      return m_memory->getPtr<V>(entry);
    }

  protected:
    //! Set mesh dimension
    void setDim(MeshDim dim);

  private:
    //! Don't permit copy operator
    Mesh(const Mesh&);
    //! Don't permit assigment operator
    Mesh& operator=(const Mesh&);

    //! Memory object pointer
    Memory* m_memory;

    //! Mesh dimension
    MeshDim m_dim;

    //! Node sets
    MeshSet m_entry;
};

} // namespace Quinoa

#endif // Mesh_h
