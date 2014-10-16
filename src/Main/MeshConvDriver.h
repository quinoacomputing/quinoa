//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.h
  \author    J. Bakosi
  \date      Mon 14 Jul 2014 09:19:44 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Mesh converter driver
  \details   Mesh converter driver
*/
//******************************************************************************
#ifndef MeshConvDriver_h
#define MeshConvDriver_h

#include <Reader.h>
#include <Writer.h>
#include <UnsMesh.h>
#include <MeshConv/CmdLine/CmdLine.h>

namespace meshconv {

//! Mesh converter driver used polymorphically with Driver
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
