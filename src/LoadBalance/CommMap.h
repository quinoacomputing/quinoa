// *****************************************************************************
/*!
  \file      src/LoadBalance/CommMap.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Calculation of communication maps for unstructured meshes
  \details   Calculation of communication maps for unstructured meshes.
*/
// *****************************************************************************
#ifndef CommMap_h
#define CommMap_h

#include <map>
#include <vector>
#include <cstddef>
#include <iosfwd>

namespace tk {

class UnsMesh;

//! Compute point-based communication maps
std::vector< std::map< std::size_t, std::vector< std::size_t > > >
poinCommMaps( std::size_t graphsize,
              const std::vector< std::size_t >& chp,
              const std::vector< std::size_t >& tetinpoel,
              std::size_t nchare,
              std::string&& toofine );

//! Compute element-based communication maps
std::vector< std::map< std::size_t, std::vector< std::size_t > > >
elemCommMaps(
  const std::vector< std::size_t >& chp,
  const std::vector< std::size_t >& tetinpoel,
  const std::vector< std::vector< std::vector< std::size_t > > >& element,
  std::size_t nchare );

} // tk::

#endif // CommMap_h
