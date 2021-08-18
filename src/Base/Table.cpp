// *****************************************************************************
/*!
  \file      src/Base/Table.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Basic functionality for storing and sampling a discrete y = f(x)
             function
  \details   Basic functionality for storing and sampling a discrete y = f(x)
             function.
*/
// *****************************************************************************

#include <utility>
#include <stddef.h>

#include "Table.hpp"
#include "Exception.hpp"

tk::real
tk::sample( tk::real x, const tk::Table& table )
// *****************************************************************************
//  Sample a discrete y = f(x) function at x
//! \param[in] x Value of abscissa at which to sample y = f(x)
//! \param[in] table tk::Table to sample
//! \details If x is lower than the first x value in the function table, the
//!   first function value is returned. If x is larger than the last x value in
//!   the function table, the last function value is returned. In other words,
//!   no extrapolation is performed. If x falls between the first/lowest and the
//!   last/largest value in the table, linear interpolation is used to compute a
//!   sample between the two closest x values of the table around the abscissa
//!   given.
//! \return Sampled value from discrete table
//! \note The x column in the table is assumed to be in increasing order.
//! \see walker::invhts_eq_A005H, walker::prod_A005H for example tables
// *****************************************************************************
{
  Assert( !table.empty(), "Empty table to sample from" );

  if (x < std::get<0>(table.front())) return std::get<1>(table.front());

  for (std::size_t i=0; i<table.size()-1; ++i) {
    if (std::get<0>(table[i]) < x && x < std::get<0>(table[i+1])) {
      auto t1 = std::get<0>( table[i] );
      auto y1 = std::get<1>( table[i] );
      auto t2 = std::get<0>( table[i+1] );
      auto y2 = std::get<1>( table[i+1] );
      return y1 + (y2-y1)/(t2-t1)*(x-t1);
    }
  }

  return std::get<1>(table.back());
}
