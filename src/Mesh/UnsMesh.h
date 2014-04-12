//******************************************************************************
/*!
  \file      src/Mesh/UnsMesh.h
  \author    J. Bakosi
  \date      Sat 12 Apr 2014 07:29:11 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     3D unstructured mesh class declaration
  \details   3D unstructured mesh class declaration
*/
//******************************************************************************
#ifndef UnshMesh_h
#define UnshMesh_h

#include <vector>
#include <memory>

#include <Types.h>
#include <Exception.h>

namespace quinoa {

//! 3D unstructured mesh class
class UnsMesh {

  public:
    //! Constructor: zero memory entry pointers held
    explicit UnsMesh() = default;

    //! Destructor, default compiler generated
    ~UnsMesh() noexcept = default;

    //! Coords accessor
    std::vector< tk::point >& coord() { return m_coord; }

    //! NodeId accessor
    std::vector< int >& nodeId() { return m_nodeId; }

    //! Line element id accessor
    std::vector< int >& lineId() { return m_lineId; }

    //! Triangle element id accessor
    std::vector< int >& triangleId() { return m_triangleId; }

    //! Tetrahedron element id accessor
    std::vector< int >& tetrahedronId() { return m_tetrahedronId; }

    //! Number of nodes accessor
    std::size_t nnode() { return m_nodeId.size(); }

    //! Line elements connectivity accessor
    std::vector< std::vector< int > >& lininpoel() { return m_lininpoel; }

    //! Line element tags accessor
    std::vector< std::vector< int > >& lintag() { return m_lintag; }

    //! Triangles elements connectivity accessor
    std::vector< std::vector< int > >& triinpoel() { return m_triinpoel; }

    //! Triangle element tags accessor
    std::vector< std::vector< int > >& tritag() { return m_tritag; }

    //! Tetrahedra elements connectivity accessor
    std::vector< std::vector< int > >& tetinpoel() { return m_tetinpoel; }

    //! Tetrahedra element tags accessor
    std::vector< std::vector< int > >& tettag() { return m_tettag; }

    //! Boundary conditions accessor
    std::vector< std::vector< int > >& bc() { return m_bc; }

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

    std::vector< tk::point > m_coord;             //!< Node coordinates
    std::vector< int > m_nodeId;                  //!< Node Ids
    std::vector< int > m_lineId;                  //!< Line element Ids
    std::vector< int > m_triangleId;              //!< Triangle element Ids
    std::vector< int > m_tetrahedronId;           //!< Tetrahedron element Ids

    std::vector< std::vector< int > > m_lininpoel;//!< Line elements conn.
    std::vector< std::vector< int > > m_lintag;   //!< Line element tags

    std::vector< std::vector< int > > m_triinpoel;//!< Triangle elements conn.
    std::vector< std::vector< int > > m_tritag;   //!< Triangle element tags

    std::vector< std::vector< int > > m_tetinpoel;//!< Tetrahedron elements conn.
    std::vector< std::vector< int > > m_tettag;   //!< Tetrahedron element tags

    std::vector< std::vector< int > > m_bc;       //!< Boundary conditions
};

} // quinoa::

#endif // UnsMesh_h
