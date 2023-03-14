// *****************************************************************************
/*!
  \file      src/Mesh/UnsMesh.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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
#include <tuple>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include "NoWarning/sip_hash.hpp"

#include "Types.hpp"
#include "ContainerUtil.hpp"

namespace tk {

//! Highway hash "secret" key
//! \note No reason for these particular numbers, taken from highwayhash tests.
static constexpr highwayhash::HH_U64 hh_key[2] =
  { 0x0706050403020100ULL, 0x0F0E0D0C0B0A0908ULL };

//! 3D unstructured mesh class
class UnsMesh {

  private:
    //! Union to access an C-style array of std::size_t as an array of bytes
    //! \tparam N Number of entries to hold
    //! \see UnsMesh::Hash
    template< std::size_t N >
    union Shaper {
      char bytes[ N*sizeof(std::size_t) ];
      std::size_t sizets[ N ];
    };

  public:
    using Coords = std::array< std::vector< real >, 3 >;
    using Coord = std::array< real, 3 >;
    using CoordMap = std::unordered_map< std::size_t, Coord >;

    //! Alias for storing a mesh chunk
    //! \details The first vector is the element connectivity (local mesh node
    //!   IDs), the second vector is the global node IDs of owned elements,
    //!   while the third one is a map of global(key)->local(value) node IDs.
    using Chunk = std::tuple< std::vector< std::size_t >,
                              std::vector< std::size_t >,
                              std::unordered_map< std::size_t, std::size_t > >;

    /** @name Aliases for element primitives */
    ///@{
    //! Edge: node IDs of two end-points
    using Edge = std::array< std::size_t, 2 >;
    //! Face: node IDs of a triangle (tetrahedron face)
    using Face = std::array< std::size_t, 3 >;
    //! Tet: node IDs of a tetrahedron
    using Tet = std::array< std::size_t, 4 >;
    ///@}

    //! Hash function class for element primitives, given by node IDs
    //! \tparam N Number of nodes describing element primitive. E.g., Edge:2,
    //!    Face:3, Tet:4.
    template< std::size_t N >
    struct Hash {
      //! Function call operator computing hash of node IDs
      //! \param[in] p Array of node IDs of element primitive
      //! \return Unique hash value for the same array of node IDs
      //! \note The order of the nodes does not matter: the IDs are sorted
      //!   before the hash is computed.
      std::size_t operator()( const std::array< std::size_t, N >& p ) const {
        using highwayhash::SipHash;
        Shaper< N > shaper;
        for (std::size_t i=0; i<N; ++i) shaper.sizets[i] = p[i];
        std::sort( std::begin(shaper.sizets), std::end(shaper.sizets) );
        return SipHash( hh_key, shaper.bytes, N*sizeof(std::size_t) );
      }
    };

    //! Comparitor function class for element primitives, given by node IDs
    //! \tparam N Number of nodes describing element primitive. E.g., Edge:2,
    //!    Face:3, Tet:4.
    template< std::size_t N >
    struct Eq {
      //! Function call operator computing equality of element primitives
      //! \param[in] l Left element primitive given by array of node IDs
      //! \param[in] r Right element primitive given by array of node IDs
      //! \return True if l = r, false otherwise
      //! \note The order of the nodes does not matter: the IDs are sorted
      //!   before equality is determined.
      bool operator()( const std::array< std::size_t, N >& l,
                       const std::array< std::size_t, N >& r ) const
      {
        std::array< std::size_t, N > s = l, p = r;
        std::sort( begin(s), end(s) );
        std::sort( begin(p), end(p) );
        return s == p;
      }
    };

    //! Unique set of edges
    using EdgeSet = std::unordered_set< Edge, Hash<2>, Eq<2> >;

    //! Unique set of faces
    using FaceSet = std::unordered_set< Face, Hash<3>, Eq<3> >;

    //! Unique set of tets
    using TetSet = std::unordered_set< Tet, Hash<4>, Eq<4> >;

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
                      const Coords& coord ) :
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

    //! \brief Constructor copying over triangle element connectivity and array
    //!    of point coordinates
    explicit UnsMesh(
      const Coords& coord,
      const std::vector< std::size_t >& triinp,
      const std::vector< std::string >& nodevarnames = {},
      const std::vector< std::string >& elemvarnames = {},
      const std::vector< tk::real >& vartimes = {},
      const std::vector< std::vector< std::vector< tk::real > > >&
        nodevars = {},
      const std::vector< std::vector< std::vector< tk::real > > >&
        elemvars = {} ) :
      m_graphsize( graphsize( triinp ) ),
      m_lininpoel(),
      m_triinpoel( triinp ),
      m_tetinpoel(),
      m_x( coord[0] ),
      m_y( coord[1] ),
      m_z( coord[2] ),
      m_nodevarnames( nodevarnames ),
      m_elemvarnames( elemvarnames ),
      m_vartimes( vartimes ),
      m_nodevars( nodevars ),
      m_elemvars( elemvars )
    {
      Assert( m_triinpoel.size()%3 == 0,
              "Size of triinpoel must be divisible by 3" );
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
    explicit UnsMesh( std::vector< std::size_t >&& tetinp, Coords&& coord ) :
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

    //! Constructor with connectivities and side set faces
    explicit UnsMesh(
      const std::vector< std::size_t >& tetinp,
      const Coords& coord,
      const std::map< int, std::vector< std::size_t > >& bf,
      const std::vector< std::size_t >& triinp,
      const std::map< int, std::vector< std::size_t > >& fid ) :
      m_graphsize( graphsize( tetinp ) ),
      m_lininpoel(),
      m_triinpoel( triinp ),
      m_tetinpoel( tetinp ),
      m_x( coord[0] ),
      m_y( coord[1] ),
      m_z( coord[2] ),
      m_bface( bf ),
      m_faceid( fid )
    {
      Assert( m_tetinpoel.size() % 4 == 0,
              "Size of tetinpoel must be divisible by 4" );
      Assert( m_triinpoel.size() % 3 == 0,
              "Size of triinpoel must be divisible by 3" );
    }

    //! Constructor with connectivities and side set nodes
    explicit UnsMesh(
      const std::vector< std::size_t >& tetinp,
      const Coords& coord,
      const std::map< int, std::vector< std::size_t > >& bn ) :
      m_graphsize( graphsize( tetinp ) ),
      m_lininpoel(),
      m_triinpoel(),
      m_tetinpoel( tetinp ),
      m_x( coord[0] ),
      m_y( coord[1] ),
      m_z( coord[2] ),
      m_bnode( bn )
    {
      Assert( m_tetinpoel.size() % 4 == 0,
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

    /** @name Triangle elements connectivity accessors */
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

    /** @name Side set face list accessors */
    ///@{
    const std::map< int, std::vector< std::size_t > >& bface() const noexcept
    { return m_bface; }
    std::map< int, std::vector< std::size_t > >& bface() noexcept
    { return m_bface; }
    ///@}

    /** @name Side set face id accessors */
    ///@{
    const std::map< int, std::vector< std::size_t > >& faceid() const noexcept
    { return m_faceid; }
    std::map< int, std::vector< std::size_t > >& faceid() noexcept
    { return m_faceid; }
    ///@}

    /** @name Side set node list accessors */
    ///@{
    const std::map< int, std::vector< std::size_t > >& bnode() const noexcept
    { return m_bnode; }
    std::map< int, std::vector< std::size_t > >& bnode() noexcept
    { return m_bnode; }
    ///@}

    /** @name Node variable names accessors */
    ///@{
    const std::vector< std::string >& nodevarnames() const noexcept
    { return m_nodevarnames; }
    std::vector< std::string >& nodevarnames() noexcept
    { return m_nodevarnames; }
    ///@}

    /** @name Element variable names accessors */
    ///@{
    const std::vector< std::string >& elemvarnames() const noexcept
    { return m_elemvarnames; }
    std::vector< std::string >& elemvarnames() noexcept
    { return m_elemvarnames; }
    ///@}

    /** @name Variable times accessors */
    ///@{
    const std::vector< tk::real >& vartimes() const noexcept
    { return m_vartimes; }
    std::vector< tk::real >& vartimes() noexcept { return m_vartimes; }
    ///@}

    /** @name Node variables accessors */
    ///@{
    const std::vector< std::vector< std::vector< tk::real > > >& nodevars()
    const noexcept { return m_nodevars; }
    std::vector< std::vector< std::vector< tk::real > > >& nodevars() noexcept
    { return m_nodevars; }
    ///@}

    /** @name Element variables accessors */
    ///@{
    const std::vector< std::vector< std::vector< tk::real > > >& elemvars()
    const noexcept { return m_elemvars; }
    std::vector< std::vector< std::vector< tk::real > > >& elemvars() noexcept
    { return m_elemvars; }
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

    //! Side sets storing face ids adjacent to side sets
    //! \details This stores lists of element IDs adjacent to faces associated
    //!   to side set IDs.
    //! \note This is what ExodusII calls side set elem list.
    std::map< int, std::vector< std::size_t > > m_bface;

    //! Side sets storing node ids adjacent to side sets
    //! \details This stores lists of node IDs adjacent to faces associated
    //!   to side set IDs.
    std::map< int, std::vector< std::size_t > > m_bnode;

    //! \brief Sides of faces used to define which face of an element is
    //!   adjacent to side set associated to side set id.
    //! \note This is what ExodusII calls side set side list.
    std::map< int, std::vector< std::size_t > > m_faceid;

    //! Node field data names
    std::vector< std::string > m_nodevarnames;

    //! Element field data names
    std::vector< std::string > m_elemvarnames;

    //! Time values for node field data
    std::vector< tk::real > m_vartimes;

    //! Node field data
    std::vector< std::vector< std::vector< tk::real > > > m_nodevars;

    //! Element field data
    std::vector< std::vector< std::vector< tk::real > > > m_elemvars;

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
