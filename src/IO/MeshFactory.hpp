// *****************************************************************************
/*!
  \file      src/IO/MeshFactory.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unstructured mesh reader and writer factory
  \details   Unstructured mesh reader and writer factory.
*/
// *****************************************************************************
#ifndef MeshFactory_h
#define MeshFactory_h

#include <iosfwd>
#include <utility>

#include "Types.hpp"
#include "UnsMesh.hpp"
#include "Print.hpp"

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
