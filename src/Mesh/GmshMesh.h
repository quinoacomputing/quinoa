//******************************************************************************
/*!
  \file      src/Mesh/GmshMesh.h
  \author    J. Bakosi
  \date      Mon Mar 24 13:57:36 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh class declaration
  \details   Gmsh mesh class declaration
*/
//******************************************************************************
#ifndef GmshMesh_h
#define GmshMesh_h

#include <vector>
#include <memory>

#include <Types.h>

namespace quinoa {

//! Gmsh mesh class
class GmshMesh {

  public:
    //! Constructor: zero memory entry pointers held
    explicit GmshMesh() noexcept :
       m_type(0),
       m_datasize(0) {}

    //! Set mesh version
    void setVersion(const tk::real version) { m_version = version; }

    //! Set mesh type
    void setType(const int type) { m_type = type; }

    //! Set mesh data size
    void setDatasize(const int datasize) { m_datasize = datasize; }

    //! Get mesh version
    tk::real getVersion() const { return m_version; }

    //! Get mesh type
    int getType() const { return m_type; }

    //! Get mesh data size
    int getDatasize() const { return m_datasize; }

    //! Coords accessor
    std::vector< tk::point >& coord() { return m_coord; }

    //! NodeId accessor
    std::vector< int >& nodeId() { return m_nodeId; }

    //! Line element id accessor
    std::vector< int >& lineId() { return m_lineId; }

    //! Triangle element id accessor
    std::vector< int >& triangleId() { return m_triangleId; }

    //! Number of nodes accessor
    int nnodes() const { return m_nnodes; }

    //! Line elements connectivity accessor
    std::vector<std::vector<int>>& linpoel() { return m_linpoel; }

    //! Line element tags accessor
    std::vector<std::vector<int>>& lintag() { return m_lintag; }

    //! Triangles elements connectivity accessor
    std::vector<std::vector<int>>& tinpoel() { return m_tinpoel; }

    //! Triangle element tags accessor
    std::vector<std::vector<int>>& tritag() { return m_tritag; }

    //! Echo element tags and connectivity in all element sets
    void echoElemSets() const;

  private:
    //! Don't permit copy constructor
    GmshMesh(const GmshMesh&) = delete;
    //! Don't permit assigment constructor
    GmshMesh& operator=(const GmshMesh&) = delete;
    //! Don't permit move constructor
    GmshMesh(GmshMesh&&) = delete;
    //! Don't permit move assignment
    GmshMesh& operator=(GmshMesh&&) = delete;

    tk::real m_version;                      //!< Mesh version in mesh file
    int m_type;                              //!< File type in mesh file
    int m_datasize;                          //!< Data size in mesh file
    int m_nnodes;                            //!< Number of nodes in mesh

    std::vector< tk::point > m_coord;        //!< Node coordinates
    std::vector< int > m_nodeId;             //!< Node Ids
    std::vector< int > m_lineId;             //!< Line element Ids
    std::vector< int > m_triangleId;         //!< Triangle element Ids

    std::vector<std::vector<int>> m_linpoel; //!< Line elements connectivity
    std::vector<std::vector<int>> m_tinpoel; //!< Triangle elements connectivity
    std::vector<std::vector<int>> m_lintag;  //!< Line element tags
    std::vector<std::vector<int>> m_tritag;  //!< Triangle element tags
};

} // quinoa::

#endif // GmshMesh_h
