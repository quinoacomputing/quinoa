// *****************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     ExodusII mesh reader
  \details   ExodusII mesh reader class declaration.
*/
// *****************************************************************************
#ifndef ExodusIIMeshReader_h
#define ExodusIIMeshReader_h

#include <cstddef>
#include <iosfwd>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>

#include "NoWarning/exodusII.hpp"

#include "Types.hpp"
#include "Exception.hpp"
#include "UnsMesh.hpp"

namespace tk {

//! Supported ExodusII mesh cell types
//! \see ExodusIIMeshReader::readElemBlockIDs()
enum class ExoElemType : int { TET = 0, TRI = 1 };

//! ExodusII mesh cell number of nodes
//! \details List of number of nodes per element for different element types
//!   supported in the order of tk::ExoElemType
const std::array< std::size_t, 2 > ExoNnpe {{ 4, 3 }};

//! ExodusII face-node numbering for tetrahedron side sets
//! \see ExodusII manual figure on "Sideset side Numbering"
const std::array< std::array< std::size_t, 3 >, 4 >
  expofa{{ {{0,1,3}}, {{1,2,3}}, {{0,3,2}}, {{0,2,1}} }};

//! ExodusII mesh-based data reader
//! \details Mesh reader class facilitating reading from mesh-based field data
//!   a file in ExodusII format.
//! \see https://github.com/trilinos/Trilinos/tree/master/packages/seacas
class ExodusIIMeshReader {

  public:
    //! Constructor
    explicit ExodusIIMeshReader( const std::string& filename,
                                 int cpuwordsize = sizeof(double),
                                 int iowordsize = sizeof(double) );

    //! Destructor
    ~ExodusIIMeshReader() noexcept;

    //! Read ExodusII mesh from file
    void readMesh( UnsMesh& mesh );

    //! Read only connectivity graph from file
    void readGraph( UnsMesh& mesh );

    //! Return total number of mesh points in mesh file
    std::size_t npoin() { return readHeader(); }

    //! Read part of the mesh (graph and coords) from file
    //! \details Total number of PEs defaults to 1 for a single-CPU read, this
    //!    PE defaults to 0 for a single-CPU read.
    void readMeshPart( std::vector< std::size_t >& ginpoel,
                       std::vector< std::size_t >& inpoel,
                       std::vector< std::size_t >& triinp,
                       std::unordered_map< std::size_t, std::size_t >& lid,
                       tk::UnsMesh::Coords& coord,
                       std::unordered_map< std::size_t, std::set< std::size_t > >&
                         elemBlockId,
                       int numpes=1, int mype=0 );

    //! Read coordinates of a number of mesh nodes from ExodusII file
    std::array< std::vector< tk::real >, 3 >
    readCoords( const std::vector< std::size_t >& gid ) const;

    //! Read face list of all side sets from ExodusII file
    void
    readSidesetFaces( std::map< int, std::vector< std::size_t > >& bface,
                      std::map< int, std::vector< std::size_t > >& faces );

    //! Read face connectivity of a number boundary faces from file
    void readFaces( std::vector< std::size_t >& conn );

    //! Read node list of all side sets from ExodusII file
    std::map< int, std::vector< std::size_t > > readSidesetNodes();

    //! Read coordinates of a single mesh node from ExodusII file
    void readNode( std::size_t fid,
                   std::size_t mid,
                   std::vector< tk::real >& x,
                   std::vector< tk::real >& y,
                   std::vector< tk::real >& z ) const;

    //! Read coordinates of a single mesh node from ExodusII file
    void readNode( std::size_t id, std::array< tk::real, 3 >& coord ) const;

    //! Read coordinates of a number of mesh nodes from ExodusII file
    std::array< std::vector< tk::real >, 3 >
    readNodes( const std::vector< std::size_t >& gid ) const;

    //! Read element block IDs from file
    std::size_t readElemBlockIDs();

    //! Read element connectivity of a number of mesh cells from file
    void readElements( const std::array< std::size_t, 2 >& ext,
                       tk::ExoElemType elemtype,
                       std::vector< std::size_t >& conn );

    //! Read local to global node-ID map
    std::vector< std::size_t > readNodemap();

    //! Generate triangle face connectivity for side sets
    std::vector< std::size_t > triinpoel(
      std::map< int, std::vector< std::size_t > >& belem,
      const std::map< int, std::vector< std::size_t > >& faces,
      const std::vector< std::size_t >& ginpoel,
      const std::vector< std::size_t >& triinp ) const;

    //! Read the names of nodal output variables from ExodusII file
    void readNodeVarNames( std::vector< std::string >& nv ) const;

    //! Read time values from ExodusII file
    void readTimeValues( std::vector< tk::real >& tv ) const;

    //! Read node scalar fields from ExodusII file
    void readNodeScalars(
      std::size_t ntime,
      std::size_t nvar,
      std::vector< std::vector< std::vector< tk::real > > >& var ) const;

    //!  Return number of elements in a mesh block in the ExodusII file
    std::size_t nelem( tk::ExoElemType elemtype ) const;

    //! Copy assignment
    // cppcheck-suppress operatorEqVarError
    // cppcheck-suppress operatorEqMissingReturnStatement
    ExodusIIMeshReader& operator=( const ExodusIIMeshReader& x ) {
      m_filename = x.m_filename;
      m_cpuwordsize = x.m_cpuwordsize;
      m_iowordsize = x.m_iowordsize;
      float version;
      m_inFile = ex_open( m_filename.c_str(), EX_READ, &m_cpuwordsize,
                          &m_iowordsize, &version );
      ErrChk( m_inFile > 0, "Failed to open ExodusII file: " + m_filename );
      m_nnode = x.m_nnode;
      m_neblk = x.m_neblk;
      m_neset = x.m_neset;
      m_from = x.m_from;
      m_till = x.m_till;
      m_blockid = x.m_blockid;
      m_blockid_by_type = x.m_blockid_by_type;
      m_nel = x.m_nel;
      m_elemblocks = x.m_elemblocks;
      m_tri = x.m_tri;
      m_elemInBlockId = x.m_elemInBlockId;
      return *this;
    }

    //! Copy constructor: in terms of copy assignment
    // cppcheck-suppress uninitMemberVar
    ExodusIIMeshReader( const ExodusIIMeshReader& x ) { operator=(x); }

    //! Move assignment
    // cppcheck-suppress operatorEqMissingReturnStatement
    // cppcheck-suppress operatorEqVarError
    ExodusIIMeshReader& operator=( ExodusIIMeshReader&& x ) {
      m_filename = x.m_filename;
      m_cpuwordsize = x.m_cpuwordsize;
      m_iowordsize = x.m_iowordsize;
      float version;
      m_inFile = ex_open( m_filename.c_str(), EX_READ, &m_cpuwordsize,
                          &m_iowordsize, &version );
      ErrChk( m_inFile > 0, "Failed to open ExodusII file: " + m_filename );
      m_nnode = x.m_nnode;
      m_neblk = x.m_neblk;
      m_neset = x.m_neset;
      m_from = x.m_from;
      m_till = x.m_till;
      m_blockid = x.m_blockid;
      m_blockid_by_type = x.m_blockid_by_type;
      m_nel = x.m_nel;
      m_elemblocks = x.m_elemblocks;
      m_tri = x.m_tri;
      m_elemInBlockId = x.m_elemInBlockId;
      x.m_cpuwordsize = sizeof(double);
      x.m_iowordsize = sizeof(double);
      x.m_inFile = ex_open( m_filename.c_str(), EX_READ, &x.m_cpuwordsize,
                            &x.m_iowordsize, &version );
      ErrChk( x.m_inFile > 0, "Failed to open ExodusII file: " + m_filename );
      x.m_nnode = 0;
      x.m_neblk = 0;
      x.m_neset = 0;
      x.m_from = 0;
      x.m_till = 0;
      x.m_blockid.clear();
      x.m_blockid_by_type.resize( ExoNnpe.size() );
      x.m_nel.resize( ExoNnpe.size() );
      x.m_elemblocks.clear();
      x.m_tri.clear();
      x.m_elemInBlockId.clear();
      return *this;
    }

    //! Move constructor: in terms of move assignment
    ExodusIIMeshReader( ExodusIIMeshReader&& x ) :
      m_filename(),
      m_cpuwordsize( 0 ),
      m_iowordsize( 0 ),
      m_inFile( 0 ),
      m_nnode( 0 ),
      m_neblk( 0 ),
      m_neset( 0 ),
      m_from( 0 ),
      m_till( 0 ),
      m_blockid(),
      m_blockid_by_type(),
      m_nel(),
      m_elemblocks(),
      m_tri(),
      m_elemInBlockId()
    { *this = std::move(x); }

  private:
    std::string m_filename;             //!< Input file name
    int m_cpuwordsize;                  //!< CPU word size for ExodusII
    int m_iowordsize;                   //!< I/O word size for ExodusII
    int m_inFile;                       //!< ExodusII file handle
    std::size_t m_nnode;                //!< Number of nodes in file
    std::size_t m_neblk;                //!< Number of element blocks in file
    std::size_t m_neset;                //!< Number of element sets in file
    std::size_t m_from;                 //!< Lower bound of tet ids on this PE
    std::size_t m_till;                 //!< Upper bound of tet ids on this PE
    //! Element block IDs in the order as in the file
    std::vector< int > m_blockid;
    //! Element block IDs for each elem type
    std::vector< std::vector< int > > m_blockid_by_type;
    //! Number of elements in blocks for each elem type
    std::vector< std::vector< std::size_t > > m_nel;
    //! Cell type and number of elements in blocks in the order as in the file
    std::vector< std::pair< ExoElemType, std::size_t > > m_elemblocks;
    //! Global->local triangle element ids on this PE
    std::unordered_map< std::size_t, std::size_t > m_tri;
    //! \brief List of elements for each block-id.
    // key: block id; value: set of elements in corresponding block
    std::unordered_map< std::size_t, std::set< std::size_t > > m_elemInBlockId;

    //! Read ExodusII header without setting mesh size
    std::size_t readHeader();

    //! Read ExodusII header with setting mesh size
    void readHeader( UnsMesh& mesh );

    //! Read coordinates of a single mesh node from file
    void readNode( std::size_t id, tk::real& x, tk::real& y, tk::real& z )
    const;

    //! Read all node coordinates from ExodusII file
    void readAllNodes( UnsMesh& mesh ) const;

    //! Read all element blocks and mesh connectivity from ExodusII file
    void readAllElements( UnsMesh& mesh );

    //! Compute element-block-relative element id and element type
    std::pair< tk::ExoElemType, std::size_t >
    blkRelElemId( std::size_t id ) const;
};

} // tk::

#endif // ExodusIIMeshReader_h
