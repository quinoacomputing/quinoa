//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.h
  \author    J. Bakosi
  \date      Wed Apr 23 10:50:16 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshConvDriver that drives MeshConv
  \details   MeshConvDriver that drives MeshConv
*/
//******************************************************************************
#ifndef MeshConvDriver_h
#define MeshConvDriver_h

#include <Driver.h>
#include <Print.h>
#include <Reader.h>
#include <Writer.h>
#include <UnsMesh.h>
#include <MeshConv/CmdLine/CmdLine.h>

namespace meshconv {

//! Mesh readers
enum class MeshReaderType : uint8_t { GMSH=0,
                                      NETGEN,
                                      EXODUSII };

//! Mesh reader factory type
using MeshReaderFactory =
  std::map< MeshReaderType, std::function<tk::Reader*()> >;

//! Mesh writers
enum class MeshWriterType : uint8_t { GMSH=0,
                                      NETGEN,
                                      EXODUSII };

//! Mesh writer factory type
using MeshWriterFactory =
  std::map< MeshWriterType, std::function<tk::Writer*()> >;

//! MeshConvDriver : Driver
class MeshConvDriver : public tk::Driver {

  public:
    //! Constructor
    explicit MeshConvDriver(int argc, char** argv, const tk::Print& print);

    //! Destructor
    ~MeshConvDriver() noexcept override = default;

    //! Execute
    void execute() const override;

  private:
    //! Don't permit copy constructor
    MeshConvDriver(const MeshConvDriver&) = delete;
    //! Don't permit assigment constructor
    MeshConvDriver& operator=(const MeshConvDriver&) = delete;
    //! Don't permit move constructor
    MeshConvDriver(MeshConvDriver&&) = delete;
    //! Don't permit move assignment
    MeshConvDriver& operator=(MeshConvDriver&&) = delete;

    //! Detect input mesh file type
    MeshReaderType detectInput() const;

    //! Determine output mesh file type
    MeshWriterType pickOutput() const;

    std::unique_ptr< ctr::CmdLine > m_cmdline;
    quinoa::UnsMesh m_mesh;             //!< Unstructured mesh to store mesh
    MeshReaderFactory m_readers;        //!< Mesh readers factory
    MeshWriterFactory m_writers;        //!< Mesh writers factory
};

} // meshconv::

#endif // MeshConvDriver_h
