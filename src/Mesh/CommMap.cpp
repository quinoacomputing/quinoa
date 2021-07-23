// *****************************************************************************
/*!
  \file      src/Mesh/CommMap.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions employing communication maps
  \details   Functions employing communication maps.
*/
// *****************************************************************************

#include <algorithm>

#include "CommMap.hpp"

namespace tk {

bool slave( const NodeCommMap& map, std::size_t node, int chare )
// *****************************************************************************
//  Decide if a node is not counted by a chare
//! \param[in] map Node commuinication map to search in
//! \param[in] node Global node id to search for
//! \param[in] chare Caller chare id (but really can be any chare id)
//! \return True if the node is a slave (counted by another chare with a lower
//!   chare id)
//! \details If a node is found in the node communication map and is associated
//! to a lower chare id than the chare id given, it is counted by another chare
//! (and not the caller one), hence a "slave" (for the purpose of this count).
// *****************************************************************************
{
  return
    std::any_of( map.cbegin(), map.cend(),
      [&](const auto& s) {
        return s.second.find(node) != s.second.cend() && s.first > chare; } );
}

tk::real count( const NodeCommMap& map, std::size_t node )
// *****************************************************************************
//  Count the number of contributions to a node
//! \param[in] map Node commuinication map to search in
//! \param[in] node Global node id to search for
//! \return Count of contributions to node
// *****************************************************************************
{
  return 1.0 + std::count_if( map.cbegin(), map.cend(),
    [&](const auto& s) { return s.second.find(node) != s.second.cend(); } );
}

} // tk::
