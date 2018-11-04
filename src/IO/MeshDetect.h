// *****************************************************************************
/*!
  \file      src/IO/MeshDetect.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Unstructured mesh file format detection
  \details   Unstructured mesh file format detection functions.
*/
// *****************************************************************************
#ifndef MeshDetect_h
#define MeshDetect_h

#include <string>

namespace tk {

//! Supported mesh readers
enum class MeshReaderType : uint8_t { GMSH=0,
                                      NETGEN,
                                      EXODUSII,
                                      HYPER,
                                      ASC,
                                      OMEGA_H };

//! Supported mesh writers
enum class MeshWriterType : uint8_t { GMSH=0,
                                      NETGEN,
                                      EXODUSII };

//! Detect input mesh file type
MeshReaderType
detectInput( const std::string& filename );

//! Determine output mesh file type
MeshWriterType
pickOutput( const std::string& filename );

} // tk::

#endif // MeshDetect_h
