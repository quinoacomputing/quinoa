//******************************************************************************
/*!
  \file      src/Mesh/DerivedData.h
  \author    J. Bakosi
  \date      Thu 12 Mar 2015 06:34:38 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Generate data structures derived from unstructured mesh
  \details   Generate data structures derived from the connectivity information
     of an unstructured mesh.
*/
//******************************************************************************
#ifndef DerivedData_h
#define DerivedData_h

#include <UnsMesh.h>

namespace tk {

//! Generate derived data structure, elements surrounding points
std::pair< std::unique_ptr< std::size_t[] >, std::unique_ptr< std::size_t[] > >
genEsup( const std::vector< std::size_t >& inpoel, std::size_t npoin );

} // tk::

#endif // DerivedData_h
