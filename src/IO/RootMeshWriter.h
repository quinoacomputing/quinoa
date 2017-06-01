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

#include <TFile.h>
#include "Types.h"

namespace tk {

class UnsMesh;

class RootMeshWriter {

  public:
    //! Constructor: create/open Root file
    explicit RootMeshWriter( const std::string filename, int option);

    //! Destructor
    ~RootMeshWriter() noexcept;

    //! Write Root mesh to file
    void writeMesh( const UnsMesh& mesh ) const;

    //! Write the names of nodal output variables to Root file
    void writeNodeVarNames( const std::vector< std::string >& nv ) const;

    //!  Write node scalar field to Root file
    void writeNodeScalar( uint64_t it,
                          int varid,
                          const std::vector< tk::real >& var ) const;

    //!  Write time stamp to Root file
    void writeTimeStamp( uint64_t it, tk::real time ) const;

  private:
    //! Write Root header
    void writeHeader( const UnsMesh& mesh ) const;

    TFile *f = 0;
    const std::string m_filename;          //!< File name
    int m_outFile;                         //!< Root file handle
};

} // tk::

#endif // RootMeshWriter_h
