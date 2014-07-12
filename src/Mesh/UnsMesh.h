//******************************************************************************
/*!
  \file      src/Mesh/UnsMesh.h
  \author    J. Bakosi
  \date      Thu 29 May 2014 05:57:42 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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
#include <Print.h>

namespace quinoa {

//! 3D unstructured mesh class
class UnsMesh {

  public:
    //! Constructor: zero memory entry pointers held
    explicit UnsMesh() = default;

    //! Destructor, default compiler generated
    ~UnsMesh() noexcept = default;

    //! Coords accessors
    std::vector< tk::real >& x() { return m_x; }
    std::vector< tk::real >& y() { return m_y; }
    std::vector< tk::real >& z() { return m_z; }

    //! NodeId accessor
    std::vector< int >& nodeId() { return m_nodeId; }

    //! Line element id accessor
    std::vector< int >& linId() { return m_linId; }

    //! Triangle element id accessor
    std::vector< int >& triId() { return m_triId; }

    //! Tetrahedron element id accessor
    std::vector< int >& tetId() { return m_tetId; }

    //! Number of nodes accessor
    std::size_t nnode() { return m_nodeId.size(); }

    //! Total number of elements accessor
    std::size_t nelem() {
      return m_lininpoel.size() + m_triinpoel.size() + m_tetinpoel.size();
    }

    //! Number of element blocks accessor
    int neblk() {
      return !m_lininpoel.empty() + !m_triinpoel.empty() + !m_tetinpoel.empty();
    }

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

    //! Echo element tags and connectivity in all element sets
    void echoElemSets( const tk::Print& print ) const;

  private:
    //! Don't permit copy constructor
    UnsMesh(const UnsMesh&) = delete;
    //! Don't permit assigment constructor
    UnsMesh& operator=(const UnsMesh&) = delete;
    //! Don't permit move constructor
    UnsMesh(UnsMesh&&) = delete;
    //! Don't permit move assignment
    UnsMesh& operator=(UnsMesh&&) = delete;

    std::vector< int > m_nodeId;                  //!< Node Ids

    std::vector< tk::real > m_x;                  //!< Node coordinates
    std::vector< tk::real > m_y;
    std::vector< tk::real > m_z;

    //! Element ids
    std::vector< int > m_linId;                   //!< Line
    std::vector< int > m_triId;                   //!< Triangle
    std::vector< int > m_tetId;                   //!< Tetrahedron

    //! Element connectivity
    std::vector< std::vector< int > > m_lininpoel;//!< Line
    std::vector< std::vector< int > > m_triinpoel;//!< Triangle
    std::vector< std::vector< int > > m_tetinpoel;//!< Tetrahedron

    //! Element tags
    std::vector< std::vector< int > > m_lintag;   //!< Line
    std::vector< std::vector< int > > m_tritag;   //!< Triangle
    std::vector< std::vector< int > > m_tettag;   //!< Tetrahedron
};

} // quinoa::

#endif // UnsMesh_h
