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

  auto eps = std::numeric_limits< tk::real >::epsilon();
  if (x < std::get<0>(table.front())+eps) return std::get<1>(table.front());

  for (std::size_t i=0; i<table.size()-1; ++i) {
    if (std::get<0>(table[i]) < x and x < std::get<0>(table[i+1])) {
      auto t1 = std::get<0>( table[i] );
      auto y1 = std::get<1>( table[i] );
      auto t2 = std::get<0>( table[i+1] );
      auto y2 = std::get<1>( table[i+1] );
      return y1 + (y2-y1)/(t2-t1)*(x-t1);
    }
  }

  return std::get<1>(table.back());
}

std::array< tk::real, 3 >
tk::sample( tk::real x, const tk::Table3& table )
// *****************************************************************************
//  Sample a discrete (y1,y2,y3) = f(x) function at x
//! \param[in] x Value of abscissa at which to sample (y1,y2,y3) = f(x)
//! \param[in] table tk::Table to sample
//! \details If x is lower than the first x value in the function table, the
//!   first 3 function values, corresponding to the first x value, are returned.
//!   If x is larger than the last x value in the function table, the last 3
//!   function values, corresponding to the last x value, are returned. In other
//!   words, no extrapolation is performed. If x falls between the first/lowest
//!   and the last/largest value in the table, linear interpolation is used to
//!   compute a sample between the two closest x values of the table around the
//!   abscissa given.
//! \return 3 sampled values from discrete table
//! \note The x column in the table is assumed to be in increasing order.
// *****************************************************************************
{
  Assert( !table.empty(), "Empty table to sample from" );

  // Lambda to return the abscissa of a 4-tuple (return the first value)
  auto abs = []( const tk::Table3::value_type& t ){ return std::get<0>(t); };

  // Lambda to return ordinates of a 4-tuple as an array of 3 values
  auto ord = []( const tk::Table3::value_type& t ){
    return std::array< tk::real, 3 >{
             std::get<1>(t), std::get<2>(t), std::get<3>(t) };
  };

  auto eps = std::numeric_limits< tk::real >::epsilon();
  if (x < abs(table.front())+eps) return ord( table.front() );

  for (std::size_t i=0; i<table.size()-1; ++i) {
    auto t1 = abs( table[i] );
    auto t2 = abs( table[i+1] );
    if (t1 < x and x < t2) {
      auto d = (t2-t1)/(x-t1);
      auto p = ord( table[i] );
      auto q = ord( table[i+1] );
      return { p[0]+(q[0]-p[0])/d, p[1]+(q[1]-p[1])/d, p[2]+(q[2]-p[2])/d };
    }
  }

  return ord( table.back() );
}
