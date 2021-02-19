// *****************************************************************************
/*!
  \file      src/IO/NetgenMeshWriter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Netgen mesh writer class declaration
  \details   Netgen mesh writer class declaration. Only supports tetrahedra.
*/
// *****************************************************************************
#ifndef NetgenMeshWriter_h
#define NetgenMeshWriter_h

#include <iosfwd>

#include "Writer.hpp"

namespace tk {

class UnsMesh;

//! Netgen mesh-based writer
//! \details Mesh reader class facilitating reading a mesh from a file saved by
//!   the Netgen mesh generator
//! \see http://sourceforge.net/apps/mediawiki/netgen-mesher
class NetgenMeshWriter : public Writer {

  public:
    //! Constructor
    explicit NetgenMeshWriter( const std::string& filename )
      : Writer( filename ) {}

    //! Write Netgen mesh
    void writeMesh( const UnsMesh& mesh );

  private:
    //! Write nodes
    void writeNodes( const UnsMesh& mesh );

    //! Write elements, i.e., connectivity
    void writeElements( const UnsMesh& mesh );
};

} // tk::

#endif // NetgenMeshWriter_h
