//******************************************************************************
/*!
  \file      src/Base/UnsMesh.h
  \author    J. Bakosi
  \date      Fri 21 Sep 2012 07:10:03 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Unstructured mesh class declaration
  \details   Unstructured mesh class declaration
*/
//******************************************************************************
#ifndef UnsMesh_h
#define UnsMesh_h

#include <vector>
#include <unordered_set>

#include <MemoryEntry.h>
#include <Memory.h>
#include <Mesh.h>

namespace Quinoa {

//! Base names for various mesh memory entries
const string LINECONN_NAME = "ConnLine";
const string  TRICONN_NAME = "ConnTri";

//! UnsMesh : Mesh
class UnsMesh : Mesh {

  public:
    //! Constructor
    UnsMesh(Memory* memory) : m_memory(memory) {}

    //! Destructor: graceful free memory entries held
    ~UnsMesh();

    //! Set mesh version
    void setVersion(Real& version) { m_version = version; }

    //! Set mesh type
    void setType(Int& type) { m_type = type; }

    //! Set mesh data size
    void setDatasize(Int& datasize) { m_datasize = datasize; }

    //! Get mesh version
    Real getVersion() { return m_version; }

    //! Get mesh type
    Int getType() { return m_type; }

    //! Get mesh data size
    Int getDatasize() { return m_datasize; }

    MemoryEntry* m_connLine = nullptr;     //!< Line connectivity
    MemoryEntry* m_connTri  = nullptr;     //!< Triangle connectivity

  private:
    //! Don't permit copy constructor
    UnsMesh(const UnsMesh&) = delete;
    //! Don't permit assigment constructor
    UnsMesh& operator=(const UnsMesh&) = delete;
    //! Don't permit move constructor
    UnsMesh(UnsMesh&&) = delete;
    //! Don't permit move assignment
    UnsMesh& operator=(UnsMesh&&) = delete;

    Real m_version;              //!< Mesh version in mesh file
    Int m_type;                  //!< File type in mesh file
    Int m_datasize;              //!< Data size in mesh file

    Memory* m_memory;            //!< Memory object pointer
};

} // namespace Quinoa

#endif // UnsMesh_h
