//******************************************************************************
/*!
  \file      src/Base/UnsMesh.h
  \author    J. Bakosi
  \date      Sun 16 Sep 2012 06:15:35 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Unstructured mesh class declaration
  \details   Unstructured mesh class declaration
*/
//******************************************************************************
#ifndef UnsMesh_h
#define UnsMesh_h

#include <vector>
#include <unordered_set>

#include <MemoryEntry.h>
#include <Memory.h>
#include <Mesh.h>

namespace Quinoa {

//! Base names for various mesh memory entries
const string   NODEID_NAME = "nodeID_";
const string    COORD_NAME = "coord_";
const string   ELEMID_NAME = "elmID_";
const string ELEMTYPE_NAME = "elmType_";

//! UnsMesh : Mesh
class UnsMesh : Mesh {

  //! Memory entry keys holding Mesh data
  typedef unordered_set<MemoryEntry*> MeshSet;

  public:
    //! Constructor
    UnsMesh(Memory* memory) : m_memory(memory), m_nodesets(0), m_elemsets(0) {}

    //! Destructor
    ~UnsMesh();

    //! Add new MeshSet entry (e.g. list of nodes, elements, node ids, etc.)
    template<class V> V* newEntry(size_t number,
                                  ValType value,
                                  VarType variable,
                                  string name,
                                  Bool plot = false,
                                  Bool restart = false) {
      // Allocate new memory entry
      MemoryEntry* entry = m_memory->newEntry(number,
                                              value,
                                              variable,
                                              name,
                                              plot,
                                              restart);
      // Store new MeshSet entry
      pair<MeshSet::iterator,Bool> n = m_entry.insert(entry);
      if (!n.second)
        throw MemoryException(ExceptType::FATAL, MemExceptType::BAD_INSERT);
      // Get pointer to new entry 
      return m_memory->getPtr<V>(entry);
    }

    //! Reserve element capacity
    void reserveElem(vector<vector<Int>>::size_type n);

    //! Add new element
    void addElem(vector<Int>& nodes);

    //! Add new element tags
    void addElemTags(vector<Int>& tags);

    //! Echo element tags and connectivity in all element sets
    void echoElemSets();

    //! Increase number of node sets
    Int addNodeSet() { return ++m_nodesets; }

    //! Increase number of element sets
    Int addElemSet() { return ++m_elemsets; }

  private:
    //! Don't permit copy constructor
    UnsMesh(const UnsMesh&) = delete;
    //! Don't permit assigment constructor
    UnsMesh& operator=(const UnsMesh&) = delete;
    //! Don't permit move constructor
    UnsMesh(UnsMesh&&) = delete;
    //! Don't permit move assignment
    UnsMesh& operator=(UnsMesh&&) = delete;

    //! Memory object pointer
    Memory* m_memory;

    //! Number of node sets
    Int m_nodesets;

    //! Number of element sets
    Int m_elemsets;

    //! Memory entry keys for Mesh data
    MeshSet m_entry;

    //! Elements
    vector<vector<Int>> m_elem;

    //! Element tags
    vector<vector<Int>> m_tag;
};

} // namespace Quinoa

#endif // UnsMesh_h
