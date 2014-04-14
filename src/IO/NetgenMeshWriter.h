//******************************************************************************
/*!
  \file      src/IO/NetgenMeshWriter.h
  \author    J. Bakosi
  \date      Mon Apr 14 11:05:05 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Netgen writer
  \details   Netgen writer
*/
//******************************************************************************
#ifndef NetgenMeshWriter_h
#define NetgenMeshWriter_h

#include <string>

#include <Writer.h>
#include <UnsMesh.h>

namespace quinoa {

//! NetgenMeshWriter : Writer
class NetgenMeshWriter : public tk::Writer {

  public:
    //! Constructor
    explicit NetgenMeshWriter( const std::string filename, UnsMesh& mesh ) :
      Writer( filename ),
      m_mesh( mesh ) {}

    //! Destructor, default compiler generated
    ~NetgenMeshWriter() noexcept override = default;

    //! Write Netgen mesh
    void write();

  private:
    //! Don't permit copy constructor
    NetgenMeshWriter(const NetgenMeshWriter&) = delete;
    //! Don't permit copy assigment
    NetgenMeshWriter& operator=(const NetgenMeshWriter&) = delete;
    //! Don't permit move constructor
    NetgenMeshWriter(NetgenMeshWriter&&) = delete;
    //! Don't permit move assigment
    NetgenMeshWriter& operator=(NetgenMeshWriter&&) = delete;

    //! Write nodes
    void writeNodes();

    //! Write elements, i.e., connectivity
    void writeElements();

    //! Write boundary conditions
    void writeBCs();

    UnsMesh& m_mesh;                   //!< Mesh object
};

} // quinoa::

#endif // NetgenMeshWriter_h
