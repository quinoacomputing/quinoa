//******************************************************************************
/*!
  \file      src/Mesh/UnsMesh.h
  \author    J. Bakosi
  \date      Fri 27 Feb 2015 03:59:46 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     3D unstructured mesh class declaration
  \details   3D unstructured mesh class declaration. This mesh class currently
    supports line, triangle, and tetrahedron elements.
*/
//******************************************************************************
#ifndef UnshMesh_h
#define UnshMesh_h

#include <vector>
#include <memory>

//#include <pup_stl.h>

#include <Types.h>
#include <Print.h>

namespace tk {

//! 3D unstructured mesh class
class UnsMesh {

  public:
    /** @name Point coordinates accessors */
    ///@{
    const std::vector< tk::real >& x() const { return m_x; }
    const std::vector< tk::real >& y() const { return m_y; }
    const std::vector< tk::real >& z() const { return m_z; }
    std::vector< tk::real >& x() { return m_x; }
    std::vector< tk::real >& y() { return m_y; }
    std::vector< tk::real >& z() { return m_z; }
    ///@}

    /** @name Node ID accessors */
    ///@{
    const std::vector< int >& nodeId() const { return m_nodeId; }
    std::vector< int >& nodeId() { return m_nodeId; }
    ///@}

    /** Line element id accessors */
    ///@{
    const std::vector< int >& linId() const { return m_linId; }
    std::vector< int >& linId() { return m_linId; }
    ///@}

    /** Triangle element id accessors */
    ///@{
    const std::vector< int >& triId() const { return m_triId; }
    std::vector< int >& triId() { return m_triId; }
    ///@}

    /** Tetrahedron element id accessors */
    ///@{
    const std::vector< int >& tetId() const { return m_tetId; }
    std::vector< int >& tetId() { return m_tetId; }
    ///@}

    /** Number of nodes accessors */
    ///@{
    std::size_t nnode() const { return m_nodeId.size(); }
    std::size_t nnode() { return m_nodeId.size(); }
    ///@}

    //! Total number of elements accessor
    std::size_t nelem() const
    { return m_lininpoel.size() + m_triinpoel.size() + m_tetinpoel.size(); }

    //! Number of element blocks accessor
    int neblk() const {
      return !m_lininpoel.empty() + !m_triinpoel.empty() + !m_tetinpoel.empty();
    }

    /** Line elements connectivity accessors */
    ///@{
    const std::vector< std::vector< int > >& lininpoel() const
    { return m_lininpoel; }
    std::vector< std::vector< int > >& lininpoel()
    { return m_lininpoel; }
    ///@}

    /** Line element tags accessors */
    ///@{
    const std::vector< std::vector< int > >& lintag() const { return m_lintag; }
    std::vector< std::vector< int > >& lintag() { return m_lintag; }
    ///@}

    /** Triangles elements connectivity accessors */
    ///@{
    const std::vector< std::vector< int > >& triinpoel() const
    { return m_triinpoel; }
    std::vector< std::vector< int > >& triinpoel()
    { return m_triinpoel; }
    ///@}

    /** Triangle element tags accessors */
    ///@{
    const std::vector< std::vector< int > >& tritag() const { return m_tritag; }
    std::vector< std::vector< int > >& tritag() { return m_tritag; }
    ///@}

    /** Tetrahedra elements connectivity accessors */
    ///@{
    const std::vector< std::vector< int > >& tetinpoel() const
    { return m_tetinpoel; }
    std::vector< std::vector< int > >& tetinpoel()
    { return m_tetinpoel; }
    ///@}

    /** Tetrahedra element tags accessors */
    ///@{
    const std::vector< std::vector< int > >& tettag() const { return m_tettag; }
    std::vector< std::vector< int > >& tettag() { return m_tettag; }
    ///@}

    //! Echo element tags and connectivity in all element sets
    void echoElemSets( const tk::Print& print ) const;

//     /** @name Pack/Unpack: Serialize unstructured mesh object for Charm++ */
//     ///@{
//     //! \brief Pack/Unpack serialize member function
//     //! \param[inout] p Charm++'s PUP::er serializer object reference
//     //! \author J. Bakosi
//     void pup( PUP::er& p ) {
//       p | m_nodeId;
//       p | m_x;
//       p | m_y;
//       p | m_z;
//       p | m_linId;
//       p | m_triId;
//       p | m_tetId;
//       p | m_lininpoel;
//       p | m_triinpoel;
//       p | m_tetinpoel;
//       p | m_lintag;
//       p | m_tritag;
//       p | m_tettag;
//     }
//     //! \brief Pack/Unpack serialize operator|
//     //! \param[inout] p Charm++'s PUP::er serializer object reference
//     //! \param[inout] m Unstructured mesh object reference
//     //! \author J. Bakosi
//     friend void operator|( PUP::er& p, UnsMesh& m ) { m.pup(p); }
//     ///@}

  private:
    //!< Node Ids
    std::vector< int > m_nodeId;

    //! Node coordinates
    std::vector< tk::real > m_x;
    std::vector< tk::real > m_y;
    std::vector< tk::real > m_z;

    //! Element ids
    std::vector< int > m_linId;                   //!< Line element ids
    std::vector< int > m_triId;                   //!< Triangle element ids
    std::vector< int > m_tetId;                   //!< Tetrahedron element ids

    //! Element connectivity
    std::vector< std::vector< int > > m_lininpoel;//!< Line connectivity
    std::vector< std::vector< int > > m_triinpoel;//!< Triangle connectivity
    std::vector< std::vector< int > > m_tetinpoel;//!< Tetrahedron connectivity

    //! Element tags
    std::vector< std::vector< int > > m_lintag;   //!< Line tags
    std::vector< std::vector< int > > m_tritag;   //!< Triangle tags
    std::vector< std::vector< int > > m_tettag;   //!< Tetrahedron tags
};

} // tk::

#endif // UnsMesh_h
