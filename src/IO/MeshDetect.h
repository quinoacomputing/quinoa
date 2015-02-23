//******************************************************************************
/*!
  \file      src/IO/MeshDetect.h
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 03:01:11 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unstructured mesh file type detector
  \details   Unstructured mesh file type detector
*/
//******************************************************************************
#ifndef MeshDetect_h
#define MeshDetect_h

#include <cstdint>
#include <string>

namespace tk {

//! Supported mesh readers
enum class MeshReaderType : uint8_t { GMSH=0,
                                      NETGEN,
                                      EXODUSII };

//! Supported mesh writers
enum class MeshWriterType : uint8_t { GMSH=0,
                                      NETGEN,
                                      EXODUSII };

//! Detect input mesh file type
MeshReaderType detectInput( const std::string& filename );

//! Determine output mesh file type
MeshWriterType pickOutput( const std::string& filename );

} // tk::

#endif // MeshDetect_h
