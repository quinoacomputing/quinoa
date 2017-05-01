// *****************************************************************************
/*!
  \file      src/IO/ExodusIIMeshWriter.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     ExodusII mesh-based data writer
  \details   ExodusII mesh-based data writer class declaration.
*/
// *****************************************************************************
#ifndef ExodusIIMeshWriter_h
#define ExodusIIMeshWriter_h

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "Types.h"

namespace tk {

class UnsMesh;

//! ExodusII writer constructor modes
enum class ExoWriter { CREATE, OPEN };

//! ExodusII mesh-based data writer
//! \details Mesh writer class facilitating writing a mesh and associated
//!   mesh-based field data to a file in ExodusII format.
//! \see http://sourceforge.net/projects/exodusii
class ExodusIIMeshWriter {

  public:
    //! Constructor: create/open ExodusII file
    explicit ExodusIIMeshWriter( const std::string& filename,
                                 ExoWriter mode,
                                 int cpuwordsize = sizeof(double),
                                 int iowordsize = sizeof(double) );

    //! Destructor
    ~ExodusIIMeshWriter() noexcept;

    //! Write ExodusII mesh to file
    void writeMesh( const UnsMesh& mesh ) const;

    //!  Write time stamp to ExodusII file
    void writeTimeStamp( uint64_t it, tk::real time ) const;

    //! Write the names of nodal output variables to ExodusII file
    void writeNodeVarNames( const std::vector< std::string >& nv ) const;

    //! Write the names of element output variables to ExodusII file
    void writeElemVarNames( const std::vector< std::string >& ev ) const;

    //!  Write node scalar field to ExodusII file
    void writeNodeScalar( uint64_t it,
                          int varid,
                          const std::vector< tk::real >& var ) const;

    //!  Write elem scalar field to ExodusII file
    void writeElemScalar( uint64_t it,
                          int varid,
                          const std::vector< tk::real >& var ) const;

  private:
    //! Write ExodusII header
    void writeHeader( const UnsMesh& mesh ) const;

    //! Write nodes coordinates to ExodusII file
    void writeNodes( const UnsMesh& mesh ) const;

    //! Write element conectivity to ExodusII file
    void writeElements( const UnsMesh& mesh ) const;

    //! Write element block to ExodusII file
    void writeElemBlock( int& elclass,
                         int64_t nnpe,
                         const std::string& eltype,
                         const std::vector< std::size_t >& inpoel ) const;

    const std::string m_filename;          //!< File name
    int m_outFile;                         //!< ExodusII file handle
};

} // tk::

#endif // ExodusIIMeshWriter_h
