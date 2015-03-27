//******************************************************************************
/*!
  \file      src/Mesh/DerivedData.h
  \author    J. Bakosi
  \date      Fri 27 Mar 2015 03:18:12 PM MDT
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

//! Shift node IDs to start with zero in element connectivity
void shiftToZero( std::vector< int >& inpoel );

//! Generate derived data structure, elements surrounding points
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsup( const std::vector< int >& inpoel, std::size_t nnpe );

//! Generate derived data structure, points surrounding points
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genPsup( const std::vector< int >& inpoel,
         std::size_t nnpe,
         const std::pair< std::vector< std::size_t >,
                          std::vector< std::size_t > >& esup );

//! Generate derived data structure, edge connectivity
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genInpoed( const std::vector< int >& inpoel,
           std::size_t nnpe,
           const std::pair< std::vector< std::size_t >,
                            std::vector< std::size_t > >& esup );

//! Generate derived data structure, elements surrounding points of elements
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsupel( const std::vector< int >& inpoel,
           std::size_t nnpe,
           const std::pair< std::vector< std::size_t >,
                            std::vector< std::size_t > >& esup );

//! Generate derived data structure, elements surrounding elements
std::vector< long int >
genEsuel( const std::vector< int >& inpoel,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esupel );

} // tk::

#endif // DerivedData_h
