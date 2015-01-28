//******************************************************************************
/*!
  \file      src/IO/NetgenMeshWriter.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 09:00:42 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Netgen mesh writer class declaration
  \details   Netgen mesh writer class declaration. Only supports tetrahedra.
*/
//******************************************************************************
#ifndef NetgenMeshWriter_h
#define NetgenMeshWriter_h

#include <string>

#include <Writer.h>
#include <UnsMesh.h>

namespace quinoa {

//! \brief NetgenMeshWriter : tk::Writer
//! Mesh reader class facilitating reading a mesh from a file saved by
//!   the Netgen mesh generator:
//!   http://sourceforge.net/apps/mediawiki/netgen-mesher.
class NetgenMeshWriter : public tk::Writer {

  public:
    //! Constructor
    explicit NetgenMeshWriter( const std::string filename,
                               const UnsMesh& mesh ) : Writer( filename ),
                                                       m_mesh( mesh ) {}

    //! Write Netgen mesh
    void write() override;

  private:
    //! Write nodes
    void writeNodes();

    //! Write elements, i.e., connectivity
    void writeElements();

    const UnsMesh& m_mesh;                   //!< Mesh object
};

} // quinoa::

#endif // NetgenMeshWriter_h
