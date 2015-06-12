//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:24:40 PM MDT
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

#include "Types.h"

namespace tk {

class UnsMesh;

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

    //  Read coordinates of a single mesh node from ExodusII file
    void readNode( std::size_t id,
                   std::vector< tk::real >& x,
                   std::vector< tk::real >& y,
                   std::vector< tk::real >& z );

  private:
    //! Read ExodusII header
    void readHeader( UnsMesh& mesh );

    //! Read all node coordinates from ExodusII file
    void readNodes( UnsMesh& mesh );

    //! Read all element blocks and connectivity from ExodusII file
    void readElements( UnsMesh& mesh );

    const std::string m_filename;          //!< File name

    int m_inFile;               //!< ExodusII file handle
    std::size_t m_neblk;        //!< Number of element blocks
    std::size_t m_nnode;        //!< Number of nodes
};

} // tk::

#endif // ExodusIIMeshReader_h
