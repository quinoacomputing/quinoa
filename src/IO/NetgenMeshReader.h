//******************************************************************************
/*!
  \file      src/IO/NetgenMeshReader.h
  \author    J. Bakosi
  \date      Thu 17 Apr 2014 07:47:12 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Netgen reader class declaration
  \details   Netgen reader class declaration
*/
//******************************************************************************
#ifndef NetgenMeshReader_h
#define NetgenMeshReader_h

#include <Reader.h>
#include <UnsMesh.h>

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

    //! Read nodes
    void readNodes();

    //! Read elements, i.e., connectivity
    void readElements();

    UnsMesh& m_mesh;                   //!< Mesh object
};

} // quinoa::

#endif // NetgenMeshReader_h
