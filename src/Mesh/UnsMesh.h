//******************************************************************************
/*!
  \file      src/Mesh/UnsMesh.h
  \author    J. Bakosi
  \date      Tue 17 Mar 2015 07:17:32 AM MDT
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
    const std::vector< tk::real >& x() const noexcept { return m_x; }
    const std::vector< tk::real >& y() const noexcept { return m_y; }
    const std::vector< tk::real >& z() const noexcept { return m_z; }
    std::vector< tk::real >& x() noexcept { return m_x; }
    std::vector< tk::real >& y() noexcept { return m_y; }
    std::vector< tk::real >& z() noexcept { return m_z; }
    ///@}

    /** @name Line element id accessors */
    ///@{
    const std::vector< int >& linId() const noexcept { return m_linId; }
    std::vector< int >& linId() noexcept { return m_linId; }
    ///@}

    /** @name Triangle element id accessors */
    ///@{
    const std::vector< int >& triId() const noexcept { return m_triId; }
    std::vector< int >& triId() noexcept { return m_triId; }
    ///@}

    /** @name Tetrahedron element id accessors */
    ///@{
    const std::vector< int >& tetId() const noexcept { return m_tetId; }
    std::vector< int >& tetId() noexcept { return m_tetId; }
    ///@}

    /** @name Number of nodes accessors */
    ///@{
    std::size_t nnode() const noexcept { return m_x.size(); }
    std::size_t nnode() noexcept { return m_x.size(); }
    ///@}

    //! Total number of elements accessor
    std::size_t nelem() const noexcept {
      return m_lininpoel.size()/2 + m_triinpoel.size()/3 + m_tetinpoel.size()/4;
    }

    //! Number of element blocks accessor
    std::size_t neblk() const noexcept {
      return !m_lininpoel.empty() + !m_triinpoel.empty() + !m_tetinpoel.empty();
    }

    /** @name Line elements connectivity accessors */
    ///@{
    const std::vector< int >& lininpoel() const noexcept { return m_lininpoel; }
    std::vector< int >& lininpoel() noexcept { return m_lininpoel; }
    ///@}

    /** @name Line element tags accessors */
    ///@{
    const std::vector< std::vector< int > >& lintag() const noexcept
    { return m_lintag; }
    std::vector< std::vector< int > >& lintag() noexcept { return m_lintag; }
    ///@}

    /** @name Triangles elements connectivity accessors */
    ///@{
    const std::vector< int >& triinpoel() const noexcept { return m_triinpoel; }
    std::vector< int >& triinpoel() noexcept { return m_triinpoel; }
    ///@}

    /** @name Triangle element tags accessors */
    ///@{
    const std::vector< std::vector< int > >& tritag() const noexcept
    { return m_tritag; }
    std::vector< std::vector< int > >& tritag() noexcept { return m_tritag; }
    ///@}

    /** @name Tetrahedra elements connectivity accessors */
    ///@{
    const std::vector< int >& tetinpoel() const noexcept { return m_tetinpoel; }
    std::vector< int >& tetinpoel() noexcept { return m_tetinpoel; }
    ///@}

    /** @name Tetrahedra element tags accessors */
    ///@{
    const std::vector< std::vector< int > >& tettag() const noexcept
    { return m_tettag; }
    std::vector< std::vector< int > >& tettag() noexcept { return m_tettag; }
    ///@}

  private:
    //! Node coordinates
    std::vector< tk::real > m_x;
    std::vector< tk::real > m_y;
    std::vector< tk::real > m_z;

    //! Element ids
    std::vector< int > m_linId;                 //!< Line element ids
    std::vector< int > m_triId;                 //!< Triangle element ids
    std::vector< int > m_tetId;                 //!< Tetrahedron element ids

    //! Element connectivity
    std::vector< int > m_lininpoel;             //!< Line connectivity
    std::vector< int > m_triinpoel;             //!< Triangle connectivity
    std::vector< int > m_tetinpoel;             //!< Tetrahedron connectivity

    //! Element tags
    std::vector< std::vector< int > > m_lintag; //!< Line tags
    std::vector< std::vector< int > > m_tritag; //!< Triangle tags
    std::vector< std::vector< int > > m_tettag; //!< Tetrahedron tags
};

} // tk::

#endif // UnsMesh_h
