// *****************************************************************************
/*!
  \file      src/IO/MeshFactory.h
  \copyright 2012-2015, J. Bakosi, 2016-2019, Los Alamos National Security, LLC.
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
