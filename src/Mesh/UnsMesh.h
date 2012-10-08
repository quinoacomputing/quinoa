//******************************************************************************
/*!
  \file      src/Base/UnsMesh.h
  \author    J. Bakosi
  \date      Mon 08 Oct 2012 12:11:11 AM EDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Unstructured mesh class declaration
  \details   Unstructured mesh class declaration
*/
//******************************************************************************
#ifndef UnsMesh_h
#define UnsMesh_h

#include <vector>

#include <Memory.h>
#include <Mesh.h>
#include <MeshException.h>

namespace Quinoa {

//! Base names for various mesh memory entries
const string     NODES_NAME = "nodes";
const string    COORDS_NAME = "coords";
const string     LINES_NAME = "lines";
const string TRIANGLES_NAME = "triangles";

//! UnsMesh : Mesh
class UnsMesh : Mesh {

  public:
    //! Constructor
    UnsMesh(Memory* memory) : m_memory(memory) {}

    //! Destructor: free memory entries held
    ~UnsMesh();

    //! Set mesh version
    void setVersion(const Real version) { m_version = version; }

    //! Set mesh type
    void setType(const Int type) { m_type = type; }

    //! Set mesh data size
    void setDatasize(const Int datasize) { m_datasize = datasize; }

    //! Get mesh version
    Real getVersion() const { return m_version; }

    //! Get mesh type
    Int getType() const { return m_type; }

    //! Get mesh data size
    Int getDatasize() const { return m_datasize; }

    //! Allocate memory for mesh
    void alloc(const Int nnodes, const Int nlins, const Int ntris);

    //! Reserve capacity to store element connectivities and tags
    void reserveElem(const Int nlines, const Int ntriangles);

    //! Add a line element
    void addLine(const vector<Int>& nodes) { m_linpoel.push_back(nodes); }

    //! Add a triangle element
    void addTriangle(const vector<Int>& nodes) { m_tinpoel.push_back(nodes); }

    //! Add line element tags
    void addLineTags(const vector<Int>& tags) { m_lintag.push_back(tags); }

    //! Add triangle element tags
    void addTriangleTags(const vector<Int>& tags) { m_tritag.push_back(tags); }

    //! Coords accessor
    Real* getCoord() const { return m_memory->getPtr<Real>(m_COORD); }

    //! NodeId accessor
    Int* getNodeId() const { return m_memory->getPtr<Int>(m_NODEID); }

    //! Line element id accessor
    Int* getLineId() const { return m_memory->getPtr<Int>(m_LINEID); }

    //! Triangle element id accessor
    Int* getTriangleId() const { return m_memory->getPtr<Int>(m_TRIANGLEID); }

    //! Number of nodes accessor
    Int getNnodes() const { return m_memory->getNumber(m_NODEID); }

    //! Line elements connectivity accessor
    const vector<vector<Int>> getLinpoel() { return m_linpoel; }

    //! Line element tags accessor
    const vector<vector<Int>> getLintag() { return m_lintag; }

    //! Triangles elements connectivity accessor
    const vector<vector<Int>> getTinpoel() { return m_tinpoel; }

    //! Triangle element tags accessor
    const vector<vector<Int>> getTritag() { return m_tritag; }

    //! Echo element tags and connectivity in all element sets
    void echoElemSets() const;

  private:
    //! Don't permit copy constructor
    UnsMesh(const UnsMesh&) = delete;
    //! Don't permit assigment constructor
    UnsMesh& operator=(const UnsMesh&) = delete;
    //! Don't permit move constructor
    UnsMesh(UnsMesh&&) = delete;
    //! Don't permit move assignment
    UnsMesh& operator=(UnsMesh&&) = delete;

    Memory* m_memory;                      //!< Memory object pointer

    Real m_version;                        //!< Mesh version in mesh file
    Int m_type;                            //!< File type in mesh file
    Int m_datasize;                        //!< Data size in mesh file

    MemoryEntry* m_COORD = nullptr;        //!< Node coordinates
    MemoryEntry* m_NODEID = nullptr;       //!< Node Ids
    MemoryEntry* m_LINEID = nullptr;       //!< Line element Ids
    MemoryEntry* m_TRIANGLEID = nullptr;   //!< Triangle element Ids

    vector<vector<Int>> m_linpoel;         //!< Line elements connectivity
    vector<vector<Int>> m_tinpoel;         //!< Triangle elements connectivity
    vector<vector<Int>> m_lintag;          //!< Line element tags
    vector<vector<Int>> m_tritag;          //!< Triangle element tags
};

} // namespace Quinoa

#endif // UnsMesh_h
