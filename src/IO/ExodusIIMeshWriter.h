//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshWriter.h
  \author    J. Bakosi
  \date      Sun 12 Apr 2015 08:20:36 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     ExodusII mesh writer
  \details   ExodusII mesh writer class declaration. Currently, this is a bare
     minimum functionality to interface with the ExodusII writer. It only writes
     3D meshes and only triangle and tetrahedron elements.
*/
//******************************************************************************
#ifndef ExodusIIMeshWriter_h
#define ExodusIIMeshWriter_h

#include <string>

#include <UnsMesh.h>
#include <Writer.h>

namespace tk {

//! \brief ExodusIIMeshWriter : tk::Writer
//! \details Mesh reader class facilitating writing a mesh to a file in
//!   ExodusII format. See also http://sourceforge.net/projects/exodusii.
class ExodusIIMeshWriter : public Writer {

  public:
    //! Constructor
    explicit ExodusIIMeshWriter( const std::string& filename,
                                 const UnsMesh& mesh,
                                 int cpuwordsize = sizeof(double),
                                 int iowordsize = sizeof(double) );

    //! Destructor
    ~ExodusIIMeshWriter() noexcept override;

    //! Write ExodusII mesh to file
    void write() override;

    //!  Write time stamp to ExodusII file
    void writeTimeStamp( int it, tk::real time );

    //! Write the number and names of output variables to ExodusII file
    void writeVarNames( const std::vector< std::string >& nvar );

    //!  Write node scalar field to ExodusII file
    void writeNodeScalar( int it,
                          int varid,
                          const std::vector< tk::real >& var );

  private:
    //! Write ExodusII header
    void writeHeader();

    //! Write nodes coordinates to ExodusII file
    void writeNodes();

    //! Write element conectivity to ExodusII file
    void writeElements();

    //! Write element block to ExodusII file
    void writeElemBlock( int& elclass,
                         int nnpe,
                         const std::string& eltype,
                         std::vector< int >& inpoel );

    const std::string m_filename;          //!< File name
    const UnsMesh& m_mesh;                 //!< Mesh object
    int m_outFile;                         //!< ExodusII file handle
};

} // tk::

#endif // ExodusIIMeshWriter_h
