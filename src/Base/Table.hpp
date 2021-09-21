// *****************************************************************************
/*!
  \file      src/Base/Table.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Basic functionality for storing and sampling a discrete
             (y1,y2,...,N) = f(x) function
  \details   Basic functionality for storing and sampling a discrete
             (y1,y2,...,N) = f(x) function.
*/
// *****************************************************************************
#ifndef Table_h
#define Table_h

#include <array>
#include <vector>

#include "Types.hpp"
#include "Exception.hpp"

namespace tk {

//! Type alias for storing a discrete (y1,y2,...,N) = f(x) function
//! \tparam N Number of ordinates in the table
template< std::size_t N >
using Table = std::vector< std::array< real, N+1 > >;

//! Sample a discrete (y1,y2,...,N) = f(x) function at x
//! \tparam N Number of ordinates in the table
template< std::size_t N >
std::array< real, N > sample( real x, const Table< N >& table ) {

  Assert( not table.empty(), "Empty table to sample from" );

  // Shortcut for the type of a single line
  using Line = std::array< tk::real, N+1 >;
  // Shortcut for the type of all ordinates
  using Ord = std::array< tk::real, N >;

  // Lambda to return the abscissa of a Table (return the first value)
  auto abscissa = []( const Line& t ){ return t[0]; };

  // Lambda to return ordinates of a tk::Table
  auto ordinate = []( const Line& t ){
    Ord o;
    for (std::size_t i=0; i<N; ++i) o[i] = t[i+1];
    return o;
  };

  auto eps = std::numeric_limits< real >::epsilon();
  if (x < abscissa(table.front())+eps) return ordinate( table.front() );

  for (std::size_t i=0; i<table.size()-1; ++i) {
    auto t1 = abscissa( table[i] );
    auto t2 = abscissa( table[i+1] );
    if (t1 < x and x < t2) {
      auto d = (t2-t1)/(x-t1);
      auto p = ordinate( table[i] );
      auto q = ordinate( table[i+1] );
      Ord r;
      for (std::size_t j=0; j<N; ++j) r[j] = p[j]+(q[j]-p[j])/d;
      return r;
    }
  }

  return ordinate( table.back() );
}

} // tk::

#endif // Table_h
