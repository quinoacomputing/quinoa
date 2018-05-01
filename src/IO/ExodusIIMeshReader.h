// *****************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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

#include "NoWarning/exodusII.h"

#include "Types.h"
#include "Exception.h"

namespace tk {

class UnsMesh;

//! \brief Supported ExodusII mesh cell types
//! \details This the order in which ExodusIIMeshReader::m_eid stores the
//!   element block IDs.
//! \see ExodusIIMeshReader::readElemBlockIDs()
enum class ExoElemType : int { TET = 0, TRI = 1 };

//! ExodusII mesh-based data reader
//! \details Mesh reader class facilitating reading from mesh-based field data
//!   a file in ExodusII format.
//! \see also http://sourceforge.net/projects/exodusii
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

    //! Read coordinates of a single mesh node from ExodusII file
    //! \param[in] fid Node id in file whose coordinates to read
    //! \param[in] mid Node id in memory to which to put new cordinates
    //! \param[in,out] x Vector of x coordinates to push to
    //! \param[in,out] y Vector of y coordinates to push to
    //! \param[in,out] z Vector of z coordinates to push to
    void readNode( std::size_t fid,
                   std::size_t mid,
                   std::vector< tk::real >& x,
                   std::vector< tk::real >& y,
                   std::vector< tk::real >& z ) const
    {
      Assert( x.size() == y.size() && x.size() == z.size(), "Size mismatch" );
      Assert( mid < x.size() && mid < y.size() && mid < z.size(),
              "Indexing out of bounds" );
      readNode( fid, x[mid], y[mid], z[mid] );
    }

    //! Read coordinates of a single mesh node from ExodusII file
    //! \param[in] id Node id whose coordinates to read
    //! \param[in,out] coord Array of x, y, and z coordinates
    void readNode( std::size_t id, std::array< tk::real, 3 >& coord ) const
    { readNode( id, coord[0], coord[1], coord[2] ); }

    //! Read coordinates of a number of mesh nodes from ExodusII file
    std::array< std::vector< tk::real >, 3 >
    readNodes( const std::vector< std::size_t >& gid ) const;

    //! Read element block IDs from file
    std::size_t readElemBlockIDs();

    //! Read element connectivity of a number of mesh cells from file
    void readElements( const std::array< std::size_t, 2 >& extent,
                       tk::ExoElemType elemtype,
                       std::vector< std::size_t >& conn ) const;

    //! Read face connectivity of a number boundary faces from file
    void readFaces( std::size_t nbfac,
                    std::vector< std::size_t >& conn );

    //! Read local to global node-ID map
    std::vector< std::size_t > readNodemap();

    //! Read node list of all side sets from ExodusII file
    std::map< int, std::vector< std::size_t > > readSidesets();

    //! Read face list of all side sets from ExodusII file
    std::size_t
    readSidesetFaces( std::map< int, std::vector< std::size_t > >& belem );

    //!  Return number of elements in a mesh block in the ExodusII file
    std::size_t nelem( tk::ExoElemType elemtype ) const;

    //! Read ExodusII header without setting mesh size
    std::size_t readHeader();

  private:
    //! Read ExodusII header with setting mesh size
    void readHeader( UnsMesh& mesh );

    //! Read coordinates of a single mesh node from file
    void readNode( std::size_t id, tk::real& x, tk::real& y, tk::real& z ) const
    {
      ErrChk(
        ex_get_partial_coord( m_inFile, static_cast<int64_t>(id)+1, 1,
                              &x, &y, &z ) == 0,
        "Failed to read coordinates of node " + std::to_string(id) +
        " from ExodusII file: " + m_filename );
    }

    //! Read all node coordinates from ExodusII file
    void readAllNodes( UnsMesh& mesh ) const;

    //! Read all element blocks and mesh connectivity from ExodusII file
    void readAllElements( UnsMesh& mesh );

    const std::string m_filename;          //!< File name
    //! \brief List of number of nodes per element for different element types
    //!   supported in the order of tk::ExoElemType
    const std::array< std::size_t, 2 > m_nnpe {{ 4, 3 }};
    int m_inFile;                       //!< ExodusII file handle
    std::size_t m_nnode;                //!< Number of nodes in file
    std::size_t m_neblk;                //!< Number of element blocks in file
    std::size_t m_neset;                //!< Number of element sets in file
    std::vector< int > m_eid;           //!< Element block IDs
    //! List of element block IDs for each elem type enum
    std::vector< std::vector< int > > m_eidt;
    //! Number of elements in blocks for each elem type enum
    std::vector< std::vector< std::size_t > > m_nel;
};

} // tk::

#endif // ExodusIIMeshReader_h
