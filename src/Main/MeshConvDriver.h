//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 11:36:44 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Mesh converter driver
  \details   Mesh converter driver.
*/
//******************************************************************************
#ifndef MeshConvDriver_h
#define MeshConvDriver_h

#include <Reader.h>
#include <Writer.h>
#include <UnsMesh.h>
#include <MeshConv/CmdLine/CmdLine.h>

//! Mesh converter declarations and definitions
namespace meshconv {

//! Mesh converter driver used polymorphically with tk::Driver
class MeshConvDriver {

  public:
    //! Constructor
    explicit MeshConvDriver( const tk::Print& print,
                             const ctr::CmdLine& cmdline );

    //! Execute
    void execute() const;

  private:
    //! Mesh readers
    enum class MeshReaderType : uint8_t { GMSH=0,
                                          NETGEN,
                                          EXODUSII };

    //! Mesh writers
    enum class MeshWriterType : uint8_t { GMSH=0,
                                          NETGEN,
                                          EXODUSII };

    //! Detect input mesh file type
    MeshReaderType detectInput() const;

    //! Determine output mesh file type
    MeshWriterType pickOutput() const;

    std::string m_input;                //!< Input file name
    std::string m_output;               //!< Output file name
};

} // meshconv::

#endif // MeshConvDriver_h
