// *****************************************************************************
/*!
  \file      src/IO/RootMeshWriter.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     ROOT mesh-based data writer
  \details   ROOT mesh-based data writer class declaration.
*/
// *****************************************************************************
#ifndef RootMeshWriter_h
#define RootMeshWriter_h

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "Types.h"

namespace tk {

class UnsMesh;

class RootMeshWriter {

  public:
    //! Constructor: create/open Root file
    explicit RootMeshWriter( const std::string& filename,
                                 RootWriter mode,
                                 int cpuwordsize = sizeof(double),
                                 int iowordsize = sizeof(double) );
    //! Destructor
    ~RootMeshWriter() noexcept;
    //! Write Root mesh to file
    void writeMesh( const UnsMesh& mesh ) const;

  private:
    //! Write Root header
    void writeHeader( const UnsMesh& mesh ) const;
    const std::string m_filename;          //!< File name
    int m_outFile;                         //!< Root file handle
};

} // tk::

#endif // RootMeshWriter_h
