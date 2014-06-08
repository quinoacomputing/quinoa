//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.h
  \author    J. Bakosi
  \date      Sun 08 Jun 2014 03:37:20 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh converter driver
  \details   Mesh converter driver
*/
//******************************************************************************
#ifndef MeshConvDriver_h
#define MeshConvDriver_h

#include <Reader.h>
#include <Writer.h>
#include <UnsMesh.h>

namespace meshconv {

//! Mesh converter driver used polymorphically with Driver
class MeshConvDriver {

  public:
    //! Constructor
    explicit MeshConvDriver( int argc, char** argv );

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
