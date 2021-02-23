// *****************************************************************************
/*!
  \file      src/IO/HyperMeshReader.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Hyper mesh reader class declaration
  \details   Hyper mesh reader class declaration. Only supports tetrahedra.
*/
// *****************************************************************************
#ifndef HyperMeshReader_h
#define HyperMeshReader_h

#include <iosfwd>
#include <utility>

#include "Reader.hpp"

namespace tk {

class UnsMesh;

//! \brief HyperMeshReader : tk::Reader
//! \details Mesh reader class facilitating reading a mesh from a file saved by
//!   the HyperMesh mesh generator by Altair:
//!    http://www.altairhyperworks.com/product/HyperMesh
class HyperMeshReader : public Reader {

  public:
    //! Constructor
    explicit HyperMeshReader( const std::string& filename ) :
      Reader( filename ) {}

    //! Read Hyper mesh
    void readMesh( UnsMesh& mesh );

  private:
    //! Read Hyper mesh metadata and extract filenames we need to read
    std::pair< std::string, std::string > getFileNames() const;

    //! Read nodes
    void readNodes( const std::string& filename, UnsMesh& mesh ) const;

    //! Read element connectivity
    void readElements( const std::string& filename, UnsMesh& mesh ) const;
};

} // tk::

#endif // HyperMeshReader_h
