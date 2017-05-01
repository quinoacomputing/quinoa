// *****************************************************************************
/*!
  \file      src/IO/MeshFactory.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unstructured mesh reader and writer factory
  \details   Unstructured mesh reader and writer factory.
*/
// *****************************************************************************
#ifndef MeshFactory_h
#define MeshFactory_h

#include <iosfwd>
#include <utility>

#include "Types.h"
#include "UnsMesh.h"
#include "Print.h"

namespace tk {

//! Supported mesh readers
enum class MeshReader : uint8_t { GMSH=0,
                                  NETGEN,
                                  EXODUSII,
                                  HYPERMESH };

//! Supported mesh writers
enum class MeshWriter : uint8_t { GMSH=0,
                                  NETGEN,
                                  EXODUSII };

//! Detect input mesh file type
MeshReader
detectInput( const std::string& filename );

//! Determine output mesh file type
MeshWriter
pickOutput( const std::string& filename );

//! Read unstructured mesh from file
UnsMesh
readUnsMesh( const tk::Print& print,
             const std::string& filename,
             std::pair< std::string, tk::real >& timestamp );

//! Write unstructured mesh to file
std::vector< std::pair< std::string, tk::real > >
writeUnsMesh( const tk::Print& print,
              const std::string& filename,
              UnsMesh& mesh,
              bool reorder );

} // tk::

#endif // MeshFactory_h
