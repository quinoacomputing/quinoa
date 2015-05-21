//******************************************************************************
/*!
  \file      src/LoadBalance/CommMap.h
  \author    J. Bakosi
  \date      Thu 21 May 2015 09:38:25 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Calculation of communication maps for unstructured meshes
  \details   Calculation of communication maps for unstructured meshes.
*/
//******************************************************************************
#ifndef CommMap_h
#define CommMap_h

#include <map>

#include <UnsMesh.h>

namespace tk {

//! Compute point-based communication maps
std::vector< std::map< std::size_t, std::vector< std::size_t > > >
poinCommMaps( const tk::UnsMesh& graph,
              const std::vector< std::size_t >& chp,
              const std::vector< std::size_t >& tetinpoel,
              std::size_t nchare,
              std::string&& toofine );

//! Compute element-based communication maps
std::vector< std::map< std::size_t, std::vector< std::size_t > > >
elemCommMaps( const std::vector< std::size_t >& chp,
              const std::vector< std::size_t >& tetinpoel,
              const std::vector< std::vector< std::size_t > >& element,
              std::size_t nchare );

} // tk::

#endif // CommMap_h
