// *****************************************************************************
/*!
  \file      src/Mesh/DerivedData.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Generate data structures derived from unstructured mesh
  \details   Generate data structures derived from the connectivity information
     of an unstructured mesh.
*/
// *****************************************************************************
#ifndef DerivedData_h
#define DerivedData_h

#include <vector>
#include <map>
#include <utility>
#include <cstddef>

namespace tk {

//! Generate derived data structure, elements surrounding points
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsup( const std::vector< std::size_t >& inpoel, std::size_t nnpe );

//! Generate derived data structure, points surrounding points
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genPsup( const std::vector< std::size_t >& inpoel,
         std::size_t nnpe,
         const std::pair< std::vector< std::size_t >,
                          std::vector< std::size_t > >& esup );

//! Generate derived data structure, edges surrounding points
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEdsup( const std::vector< std::size_t >& inpoel,
          std::size_t nnpe,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup );

//! Generate derived data structure, edge connectivity
std::vector< std::size_t >
genInpoed( const std::vector< std::size_t >& inpoel,
           std::size_t nnpe,
           const std::pair< std::vector< std::size_t >,
                            std::vector< std::size_t > >& esup );

//! Generate derived data structure, elements surrounding points of elements
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsupel( const std::vector< std::size_t >& inpoel,
           std::size_t nnpe,
           const std::pair< std::vector< std::size_t >,
                            std::vector< std::size_t > >& esup );

//! Generate derived data structure, elements surrounding elements
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsuel( const std::vector< std::size_t >& inpoel,
          std::size_t nnpe,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup );

//! Generate derived data structure, elements surrounding elements
//  for tetrahedra in the form of a full vector with voundary elements
std::vector< int >
genEsuelTet( const std::vector< std::size_t >& inpoel,
             const std::pair< std::vector< std::size_t >,
                              std::vector< std::size_t > >& esup );

//! Generate derived data structure, edges of elements
std::vector< std::size_t >
genInedel( const std::vector< std::size_t >& inpoel,
           std::size_t nnpe,
           const std::vector< std::size_t >& inpoed );

//! Generate derived data structure, elements surrounding edges
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsued( const std::vector< std::size_t >& inpoel,
          std::size_t nnpe,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup );

//! Generate derived data structure, total number of faces
std::size_t 
genNtfac( std::size_t nfpe,
          const std::size_t& nbfac,
          const std::vector< int >& esuelTet );

//! Generate derived data structure, elements surrounding faces
std::vector< int >
genEsuf( std::size_t nfpe,
         const std::size_t& ntfac,
         const std::size_t& nbfac,
         const std::map< int, std::vector< std::size_t > >& belem,
         const std::vector< int >& esuelTet );

//! Generate derived data structure, node-face connectivity
std::vector< std::size_t >
genInpofaTet( const std::size_t& ntfac,
              const std::size_t& nbfac,
              const std::vector< std::size_t >& inpoel,
              const std::vector< int >& esuelTet );
} // tk::

#endif // DerivedData_h
