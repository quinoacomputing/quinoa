//******************************************************************************
/*!
  \file      src/IO/NetgenMeshReader.h
  \author    J. Bakosi
  \date      Thu 10 Apr 2014 09:46:54 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Netgen reader class declaration
  \details   Netgen reader class declaration
*/
//******************************************************************************
#ifndef NetgenMeshReader_h
#define NetgenMeshReader_h

#include <vector>
#include <map>

#include <Reader.h>
#include <UnsMesh.h>
#include <Exception.h>

namespace quinoa {

//! NetgenMeshReader : Reader
class NetgenMeshReader : public tk::Reader {

  public:
    //! Constructor
    explicit NetgenMeshReader( const std::string filename, UnsMesh& mesh ) :
      Reader( filename ),
      m_mesh( mesh ) {}

    //! Destructor, default compiler generated
    ~NetgenMeshReader() noexcept override = default;

    //! Read Netgen mesh
    void read() override;

  private:
    //! Don't permit copy constructor
    NetgenMeshReader(const NetgenMeshReader&) = delete;
    //! Don't permit copy assigment
    NetgenMeshReader& operator=(const NetgenMeshReader&) = delete;
    //! Don't permit move constructor
    NetgenMeshReader(NetgenMeshReader&&) = delete;
    //! Don't permit move assigment
    NetgenMeshReader& operator=(NetgenMeshReader&&) = delete;

    //! Read section
    void readNodes();

    //! Read section
    void readElements();

    UnsMesh& m_mesh;                   //!< Mesh object
};

} // quinoa::

#endif // NetgenMeshReader_h
