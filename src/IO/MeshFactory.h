//******************************************************************************
/*!
  \file      src/IO/MeshFactory.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:25:18 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unstructured mesh reader and writer factory
  \details   Unstructured mesh reader and writer factory.
*/
//******************************************************************************
#ifndef MeshFactory_h
#define MeshFactory_h

#include <cstdint>
#include <string>

#include "UnsMesh.h"

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

//! Read unstructured mesh from file
UnsMesh readUnsMesh( const std::string& filename,
                     std::pair< std::string, tk::real >& timestamp );

//! Write unstructured mesh to file
void writeUnsMesh( const std::string& filename,
                   const UnsMesh& mesh,
                   std::pair< std::string, tk::real >& timestamp );

} // tk::

#endif // MeshFactory_h
