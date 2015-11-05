//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.h
  \author    J. Bakosi
  \date      Thu 29 Oct 2015 03:47:46 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     ExodusII mesh reader
  \details   ExodusII mesh reader class declaration.
*/
//******************************************************************************
#ifndef ExodusIIMeshReader_h
#define ExodusIIMeshReader_h

#include <cstddef>
#include <iosfwd>
#include <vector>
#include <array>

#include "Types.h"

namespace tk {

class UnsMesh;

//! \brief Supported ExodusII mesh cell types
//! \details This the order in which ExodusIIMeshReader::m_eid stores the
//!   element block IDs.
//! \see ExodusIIMeshReader::readElemBlockIDs() and
//!    ExodusIIMeshReader::readElement
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

    //! Read coordinates of a single mesh node from file
    void readNode( std::size_t id,
                   std::vector< tk::real >& x,
                   std::vector< tk::real >& y,
                   std::vector< tk::real >& z );

    //! Read element block IDs from file
    void readElemBlockIDs();

    //! Read element connectivity of a single mesh cell from file
    void readElement( std::size_t id,
                      tk::ExoElemType elemtype,
                      std::vector< std::size_t >& conn );

  private:
    //! Read ExodusII header without setting mesh size
    int readHeader();

    //! Read ExodusII header with setting mesh size
    void readHeader( UnsMesh& mesh );

    //! Read all node coordinates from ExodusII file
    void readNodes( UnsMesh& mesh );

    //! Read all element blocks and connectivity from ExodusII file
    void readElements( UnsMesh& mesh );

    const std::string m_filename;          //!< File name

    //! \brief List of number of nodes per element for different element types
    //!   supported in the order of tk::ExoElemType
    const std::array< std::size_t, 2 > m_nnpe {{ 4, 3 }};

    int m_inFile;                       //!< ExodusII file handle
    std::size_t m_nnode;                //!< Number of nodes
    std::size_t m_neblk;                //!< Number of element blocks
    std::vector< int > m_eid;           //!< Element block IDs
    std::vector< int > m_eidt;          //!< Element block IDs mapped to enum
};

} // tk::

#endif // ExodusIIMeshReader_h
