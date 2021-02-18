// *****************************************************************************
/*!
  \file      src/IO/ASCMeshReader.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     ASC mesh reader class declaration
  \details   ASC mesh reader class declaration. Mesh reader facilitating reading
             a mesh from a simple text file used by Jacob Waltz's Chicoma code.
*/
// *****************************************************************************
#ifndef ASCMeshReader_h
#define ASCMeshReader_h

#include <iosfwd>

#include "Reader.hpp"

namespace tk {

class UnsMesh;

//! \brief ASCMeshReader : tk::Reader
//! \details Mesh reader class facilitating reading a mesh from a simple test
//!   file used by Jacob Waltz's Chicoma code.
class ASCMeshReader : public Reader {

  public:
    //! Constructor
    explicit ASCMeshReader( const std::string& filename ) :
      Reader( filename ) {}

    //! Read ASC mesh
    void readMesh( UnsMesh& mesh );

  private:
    //! Read header
    void readHeader();

    //! Read nodes
    void readNodes( UnsMesh& mesh );

    //! Read element connectivity
    void readElements( UnsMesh& mesh );
};

} // tk::

#endif // ASCMeshReader_h
