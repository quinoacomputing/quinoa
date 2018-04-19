// *****************************************************************************
/*!
  \file      src/Mesh/UnsMesh.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     3D unstructured mesh class declaration
  \details   3D unstructured mesh class declaration. This mesh class currently
    supports line, triangle, and tetrahedron elements.
*/
// *****************************************************************************
#ifndef UnsMesh_h
#define UnsMesh_h

#include <vector>
#include <array>
#include <memory>
#include <unordered_set>
#include <unordered_map>

#include "Types.h"
#include "ContainerUtil.h"

namespace tk {

//! 3D unstructured mesh class
class UnsMesh {

  public:
    using Coords = std::array< std::vector< real >, 3 >;

    //! Edge: IDs of two end-points
    using Edge = std::array< std::size_t, 2 >;
    //! Hash functor for Edge (node end-point order does not matter)
    struct EdgeHash {
      std::size_t operator()( const Edge& key ) const {
        return std::hash< std::size_t >()( key[0] ) ^
               std::hash< std::size_t >()( key[1] );
      }
    };
    //! Key-equal function for Edge (node end-point order does not matter)
    struct EdgeEq {
      bool operator()( const Edge& left, const Edge& right ) const {
        return (left[0] == right[0] && left[1] == right[1]) ||
               (left[0] == right[1] && left[1] == right[0]);
      }
    };
    //! Map associating IDs of mesh nodes to edges
    using EdgeNodes = std::unordered_map< Edge, std::size_t, EdgeHash, EdgeEq >;
    using EdgeChares = std::unordered_map< Edge,
                                           std::vector< int >,
                                           EdgeHash,
                                           EdgeEq >;
    //! Unique set of edges
    using Edges = std::unordered_set< Edge, EdgeHash, EdgeEq >;

    //! Node ID triplet denoting a tetrahedron face
    using Face = std::array< std::size_t, 3 >;
    // Hash functor for node triplet (Face)
    struct FaceHasher {
      std::size_t operator()( const Face& key ) const {
        return std::hash< std::size_t >()( key[0] ) ^
               std::hash< std::size_t >()( key[1] ) ^
               std::hash< std::size_t >()( key[2] );
      }
    };
    //! Key-equal function for node triplet (Face)
    //! \note This requires all 3 node IDs be present in both the left and right
    //!   hand side but ignores ordering. For example, for A,B,C it will match
    //!   B,A,C, C,B,A, etc.
    struct FaceEq {
      bool operator()( const Face& l, const Face& r ) const {
        return (l[0] == r[0] || l[0] == r[1] || l[0] == r[2]) &&
               (l[1] == r[0] || l[1] == r[1] || l[1] == r[2]) &&
               (l[2] == r[0] || l[2] == r[1] || l[2] == r[2]);
      }
    };
    //! Unique set of node ID triplets representing unique tetrahedron faces
    using FaceSet = std::unordered_set< Face, FaceHasher, FaceEq >;

    /** @name Constructors */
    ///@{
    //! Constructor without initializing anything
    explicit UnsMesh() : m_graphsize(0), m_lininpoel(), m_triinpoel(),
      m_tetinpoel(), m_x(), m_y(), m_z() {}

    //! Constructor copying over element connectivity
    explicit UnsMesh( const std::vector< std::size_t >& tetinp ) :
      m_graphsize( graphsize( tetinp ) ),
      m_lininpoel(), m_triinpoel(),
      m_tetinpoel( tetinp ),
      m_x(), m_y(), m_z()
    {
      Assert( m_tetinpoel.size()%4 == 0,
              "Size of tetinpoel must be divisible by 4" );
    }

    //! Constructor swallowing element connectivity
    explicit UnsMesh( std::vector< std::size_t >&& tetinp ) :
      m_graphsize( graphsize( tetinp ) ),
      m_lininpoel(), m_triinpoel(),
      m_tetinpoel( std::move(tetinp) ),
      m_x(), m_y(), m_z()
    {
      Assert( m_tetinpoel.size()%4 == 0,
              "Size of tetinpoel must be divisible by 4" );
    }

    //! Constructor copying over element connectivity and point coordinates
    explicit UnsMesh( const std::vector< std::size_t >& tetinp,
                      const std::vector< real >& X,
                      const std::vector< real >& Y,
                      const std::vector< real >& Z ) :
      m_graphsize( graphsize( tetinp ) ),
      m_lininpoel(), m_triinpoel(),
      m_tetinpoel( tetinp ),
      m_x( X ),
      m_y( Y ),
      m_z( Z )
    {
      Assert( m_tetinpoel.size()%4 == 0,
              "Size of tetinpoel must be divisible by 4" );
    }

    //! \brief Constructor copying over element connectivity and array of point
    //!   coordinates
    explicit UnsMesh( const std::vector< std::size_t >& tetinp,
                      const std::array< std::vector< real >, 3 >& coord ) :
      m_graphsize( graphsize( tetinp ) ),
      m_lininpoel(), m_triinpoel(),
      m_tetinpoel( tetinp ),
      m_x( coord[0] ),
      m_y( coord[1] ),
      m_z( coord[2] )
    {
      Assert( m_tetinpoel.size()%4 == 0,
              "Size of tetinpoel must be divisible by 4" );
    }

    //! Constructor swallowing element connectivity and point coordinates
    explicit UnsMesh( std::vector< std::size_t >&& tetinp,
                      std::vector< real >&& X,
                      std::vector< real >&& Y,
                      std::vector< real >&& Z ) :
      m_graphsize( graphsize( tetinp ) ),
      m_lininpoel(), m_triinpoel(),
      m_tetinpoel( std::move(tetinp) ),
      m_x( std::move(X) ),
      m_y( std::move(Y) ),
      m_z( std::move(Z) )
    {
      Assert( m_tetinpoel.size()%4 == 0,
              "Size of tetinpoel must be divisible by 4" );
    }


    //! \brief Constructor swallowing element connectivity and array of point
    //!   coordinates
    explicit UnsMesh( std::vector< std::size_t >&& tetinp,
                      std::array< std::vector< real >, 3 >&& coord ) :
      m_graphsize( graphsize( tetinp ) ),
      m_lininpoel(), m_triinpoel(),
      m_tetinpoel( std::move(tetinp) ),
      m_x( std::move(coord[0]) ),
      m_y( std::move(coord[1]) ),
      m_z( std::move(coord[2]) )
    {
      Assert( m_tetinpoel.size()%4 == 0,
              "Size of tetinpoel must be divisible by 4" );
    }
    ///@}

    /** @name Point coordinates accessors */
    ///@{
    const std::vector< real >& x() const noexcept { return m_x; }
    const std::vector< real >& y() const noexcept { return m_y; }
    const std::vector< real >& z() const noexcept { return m_z; }
    std::vector< real >& x() noexcept { return m_x; }
    std::vector< real >& y() noexcept { return m_y; }
    std::vector< real >& z() noexcept { return m_z; }
    ///@}

    /** @name Number of nodes accessors */
    ///@{
    std::size_t nnode() const noexcept { return m_x.size(); }
    std::size_t nnode() noexcept { return m_x.size(); }
    ///@}

    /** @name Graph size accessors */
    ///@{
    const std::size_t& size() const noexcept { return m_graphsize; }
    std::size_t& size() noexcept { return m_graphsize; }
    ///@}

    //! Total number of elements accessor
    std::size_t nelem() const noexcept {
      return m_lininpoel.size()/2 + m_triinpoel.size()/3 + m_tetinpoel.size()/4;
    }

    //! Number of element blocks accessor
    std::size_t neblk() const noexcept {
      return !m_triinpoel.empty() + !m_tetinpoel.empty();
    }

    /** @name Line elements connectivity accessors */
    ///@{
    const std::vector< std::size_t >& lininpoel() const noexcept
    { return m_lininpoel; }
    std::vector< std::size_t >& lininpoel() noexcept { return m_lininpoel; }
    ///@}

    /** @name Triangles elements connectivity accessors */
    ///@{
    const std::vector< std::size_t >& triinpoel() const noexcept
    { return m_triinpoel; }
    std::vector< std::size_t >& triinpoel() noexcept { return m_triinpoel; }
    ///@}

    /** @name Tetrahedra elements connectivity accessors */
    ///@{
    const std::vector< std::size_t >& tetinpoel() const noexcept
    { return m_tetinpoel; }
    std::vector< std::size_t >& tetinpoel() noexcept { return m_tetinpoel; }
    ///@}

  private:
    //! Number of nodes
    //! \details Stores the size (number of nodes) of the mesh graph.
    //!   Used if only the graph is needed but not the node coordinates, e.g.,
    //!   for graph partitioning, in which case only the connectivity is
    //!   required. If the coordinates are also loaded, the member functions
    //!   nnode() and size() return the same.
    std::size_t m_graphsize;

    //! Element connectivity
    std::vector< std::size_t > m_lininpoel;     //!< Line connectivity
    std::vector< std::size_t > m_triinpoel;     //!< Triangle connectivity
    std::vector< std::size_t > m_tetinpoel;     //!< Tetrahedron connectivity

    //! Node coordinates
    std::vector< real > m_x;
    std::vector< real > m_y;
    std::vector< real > m_z;

    //! Compute and return number of unique nodes in element connectivity
    //! \param[in] inpoel Element connectivity
    //! \return Number of unique node ids in connectivity, i.e., the graphsize
    std::size_t
    graphsize( const std::vector< std::size_t >& inpoel ) {
      auto conn = inpoel;
      unique( conn );
      return conn.size();
   }
};

} // tk::

#endif // UnsMesh_h
