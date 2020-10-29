// *****************************************************************************
/*!
  \file      src/Inciter/FieldOutputUtil.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field output utility functionality
  \details   Field output utility functions.
*/
// *****************************************************************************
#ifndef FieldOutputUtil_h
#define FieldOutputUtil_h

#include "Fields.hpp"
#include "Vector.hpp"
#include "UnsMesh.hpp"

namespace tk {

//! Evaluate DG solution in nodes
std::tuple< tk::Fields, tk::Fields >
nodeEval( std::size_t offset,
          std::size_t ndof,
          std::size_t rdof,
          const tk::UnsMesh::Coords& coord,
          const std::vector< std::size_t >& inpoel,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup,
          const Fields& U,
          const Fields& P = tk::Fields() );

} // tk::

#endif // FieldOutputUtil_h
