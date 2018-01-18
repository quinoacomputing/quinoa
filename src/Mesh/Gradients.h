// *****************************************************************************
/*!
  \file      src/Mesh/Gradients.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Functions computing gradients on unstructured meshes for tetrahedra
  \details   Functions computing gradients using linear finite element shape
             functions on unstructured meshes for tetrahedra.
*/
// *****************************************************************************
#ifndef Gradients_h
#define Gradients_h

#include <array>
#include <vector>
#include <utility>

#include "Fields.h"
#include "Keywords.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute gradient at a mesh node
std::array< tk::real, 3 >
nodegrad( std::size_t node,
          const std::array< std::vector< tk::real >, 3 >& coord,
          const std::vector< std::size_t >& inpoel,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup,
          const tk::Fields& U,
          ncomp_t c );

//! Compute gradient at a mesh edge
std::array< tk::real, 3 >
edgegrad( std::size_t edge,
          const std::array< std::vector< tk::real >, 3 >& coord,
          const std::vector< std::size_t >& inpoel,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esued,
          const tk::Fields& U,
          ncomp_t c );

} // tk::

#endif // Gradients_h
