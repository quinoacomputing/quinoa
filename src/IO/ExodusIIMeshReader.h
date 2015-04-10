//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.h
  \author    J. Bakosi
  \date      Fri 13 Mar 2015 05:08:19 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     ExodusII mesh reader
  \details   ExodusII mesh reader class declaration. Currently, this is a bare
     minimum functionality to interface with the ExodusII reader. It only reads
     3D meshes and only triangle and tetrahedron elements.
*/
//******************************************************************************
#ifndef ExodusIIMeshReader_h
#define ExodusIIMeshReader_h

#include <string>

#include <UnsMesh.h>
#include <Reader.h>

namespace tk {

//! \brief ExodusIIMeshReader : tk::Reader
//! \details Mesh reader class facilitating reading a mesh from a file in
//!   ExodusII format. See also http://sourceforge.net/projects/exodusii.
class ExodusIIMeshReader : public Reader {

  public:
    //! Constructor
    explicit ExodusIIMeshReader( const std::string& filename,
                                 UnsMesh& mesh,
                                 int cpuwordsize = sizeof(double),
                                 int iowordsize = sizeof(double) );

    //! Destructor
    ~ExodusIIMeshReader() noexcept override;

    //! Read ExodusII mesh from file
    void read() override;

    //! Read only connectivity graph from file
    void readGraph();

    //  Read the coordinates of a single mesh node from ExodusII file
    void readNode( std::size_t id,
                   std::vector< tk::real >& x,
                   std::vector< tk::real >& y,
                   std::vector< tk::real >& z );

  private:
    //! Read ExodusII header
    void readHeader();

    //! Read node coordinates from ExodusII file
    void readNodes();

    //! Read element blocks and connectivity from ExodusII file
    void readElements();

    UnsMesh& m_mesh;            //!< Mesh object
    int m_inFile;               //!< ExodusII file handle
    int m_neblk;                //!< Number of element blocks
    int m_nnode;                //!< Number of nodes
};

} // tk::

#endif // ExodusIIMeshReader_h
