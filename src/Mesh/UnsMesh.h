//******************************************************************************
/*!
  \file      src/Mesh/UnsMesh.h
  \author    J. Bakosi
  \date      Fri Apr 26 15:46:26 2013
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

namespace Quinoa {

//! Base names for various mesh memory entries
const string     NODES_NAME = "nodes";
const string    COORDS_NAME = "coords";
const string     LINES_NAME = "lines";
const string TRIANGLES_NAME = "triangles";

//! UnsMesh : Mesh
class UnsMesh : private Mesh {

  public:
    //! Constructor: zero memory entry pointers held
    explicit UnsMesh(Memory* const memory);

    //! Destructor: free memory entries held
    virtual ~UnsMesh() noexcept;

    //! Set mesh version
    void setVersion(const real version) { m_version = version; }

    //! Set mesh type
    void setType(const int type) { m_type = type; }

    //! Set mesh data size
    void setDatasize(const int datasize) { m_datasize = datasize; }

    //! Get mesh version
    real getVersion() const { return m_version; }

    //! Get mesh type
    int getType() const { return m_type; }

    //! Get mesh data size
    int getDatasize() const { return m_datasize; }

    //! Allocate memory for mesh
    void alloc(const int nnodes, const int nlins, const int ntris);

    //! Reserve capacity to store element connectivities and tags
    void reserveElem(const int nlines, const int ntriangles);

    //! Add a line element
    void addLine(const vector<int>& nodes) { m_linpoel.push_back(nodes); }

    //! Add a triangle element
    void addTriangle(const vector<int>& nodes) { m_tinpoel.push_back(nodes); }

    //! Add line element tags
    void addLineTags(const vector<int>& tags) { m_lintag.push_back(tags); }

    //! Add triangle element tags
    void addTriangleTags(const vector<int>& tags) { m_tritag.push_back(tags); }

    //! Coords accessor
    real* getCoord() const { return m_coord.ptr; }

    //! NodeId accessor
    int* getNodeId() const { return m_nodeId.ptr; }

    //! Line element id accessor
    int* getLineId() const { return m_lineId.ptr; }

    //! Triangle element id accessor
    int* getTriangleId() const { return m_triangleId.ptr; }

    //! Number of nodes accessor
    int getNnodes() const { return m_memory->getNumber(m_nodeId.id); }

    //! Line elements connectivity accessor
    const vector<vector<int>> getLinpoel() { return m_linpoel; }

    //! Line element tags accessor
    const vector<vector<int>> getLintag() { return m_lintag; }

    //! Triangles elements connectivity accessor
    const vector<vector<int>> getTinpoel() { return m_tinpoel; }

    //! Triangle element tags accessor
    const vector<vector<int>> getTritag() { return m_tritag; }

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

    Memory* const m_memory;                //!< Memory object pointer

    real m_version;                        //!< Mesh version in mesh file
    int m_type;                            //!< File type in mesh file
    int m_datasize;                        //!< Data size in mesh file

    Data<real> m_coord;                     //!< Node coordinates
    Data<int> m_nodeId;                    //!< Node Ids
    Data<int> m_lineId;                    //!< Line element Ids
    Data<int> m_triangleId;                //!< Triangle element Ids

    vector<vector<int>> m_linpoel;         //!< Line elements connectivity
    vector<vector<int>> m_tinpoel;         //!< Triangle elements connectivity
    vector<vector<int>> m_lintag;          //!< Line element tags
    vector<vector<int>> m_tritag;          //!< Triangle element tags
};

} // namespace Quinoa

#endif // UnsMesh_h
