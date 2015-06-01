//******************************************************************************
/*!
  \file      src/IO/NetgenMeshWriter.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:24:22 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Netgen mesh writer class declaration
  \details   Netgen mesh writer class declaration. Only supports tetrahedra.
*/
//******************************************************************************
#ifndef NetgenMeshWriter_h
#define NetgenMeshWriter_h

#include <string>

#include "Writer.h"
#include "UnsMesh.h"

namespace tk {

//! \brief NetgenMeshWriter : tk::Writer
//! Mesh reader class facilitating reading a mesh from a file saved by
//!   the Netgen mesh generator:
//!   http://sourceforge.net/apps/mediawiki/netgen-mesher.
class NetgenMeshWriter : public Writer {

  public:
    //! Constructor
    explicit NetgenMeshWriter( const std::string filename,
                               const UnsMesh& mesh )
      : Writer( filename ), m_mesh( mesh ) {}

    //! Write Netgen mesh
    void write() override;

  private:
    //! Write nodes
    void writeNodes();

    //! Write elements, i.e., connectivity
    void writeElements();

    const UnsMesh& m_mesh;                   //!< Mesh object
};

} // tk::

#endif // NetgenMeshWriter_h
