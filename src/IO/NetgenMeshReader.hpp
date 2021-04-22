// *****************************************************************************
/*!
  \file      src/IO/NetgenMeshReader.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Netgen mesh reader class declaration
  \details   Netgen mesh reader class declaration. Only supports tetrahedra.
*/
// *****************************************************************************
#ifndef NetgenMeshReader_h
#define NetgenMeshReader_h

#include <iosfwd>

#include "Reader.hpp"

namespace tk {

class UnsMesh;

//! \brief NetgenMeshReader : tk::Reader
//! \details Mesh reader class facilitating reading a mesh from a file saved by
//!   the Netgen mesh generator:
//!   http://sourceforge.net/apps/mediawiki/netgen-mesher.
class NetgenMeshReader : public Reader {

  public:
    //! Constructor
    explicit NetgenMeshReader( const std::string& filename ) :
      Reader( filename ) {}

    //! Read Netgen mesh
    void readMesh( UnsMesh& mesh );

  private:
    //! Read nodes
    void readNodes( UnsMesh& mesh );

    //! Read element connectivity
    void readElements( UnsMesh& mesh );
};

} // tk::

#endif // NetgenMeshReader_h
