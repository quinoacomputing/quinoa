//******************************************************************************
/*!
  \file      src/Mesh/Mesh.h
  \author    J. Bakosi
  \date      Fri Sep 14 17:48:09 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh base class declaration
  \details   Mesh base class declaration
*/
//******************************************************************************
#ifndef Mesh_h
#define Mesh_h

#include <vector>
#include <unordered_set>

#include <MemoryEntry.h>
#include <Memory.h>

namespace Quinoa {

//! Base names for various mesh memory entries
const string   NODEID_NAME = "nodeID_";
const string    COORD_NAME = "coord_";
const string   ELEMID_NAME = "elmID_";
const string ELEMTYPE_NAME = "elmType_";

//! Mesh dimension
enum class MeshDim { TWOD=0,
                     THREED,
                     NUM_MESH_DIM
};
//! Number of mesh dimensions
const Int NUM_MESH_DIM = static_cast<Int>(MeshDim::NUM_MESH_DIM);

//! Mesh base class
class Mesh {

  //! Memory entry keys holding Mesh data
  typedef unordered_set<MemoryEntry*> MeshSet;

  public:
    //! Constructor
    Mesh(Memory* memory) : m_memory(memory), m_nodesets(0), m_elemsets(0) {}

    //! Destructor
    virtual ~Mesh();

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
    void reserveElem(vector< vector<Int> >::size_type n);

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

    //! Number of node sets
    Int m_nodesets;

    //! Number of element sets
    Int m_elemsets;

    //! Memory entry keys for Mesh data
    MeshSet m_entry;

    //! Elements
    vector< vector<Int> > m_elem;

    //! Element tags
    vector< vector<Int> > m_tag;
};

} // namespace Quinoa

#endif // Mesh_h
