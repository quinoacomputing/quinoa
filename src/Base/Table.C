// *****************************************************************************
/*!
  \file      src/Base/Table.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Basic functionality for storing and sampling a discrete y = f(x)
             function
  \details   Basic functionality for storing and sampling a discrete y = f(x)
             function.
*/
// *****************************************************************************

#include <vector>
#include <utility>

#include "Table.h"

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
//! \author J. Bakosi
// *****************************************************************************
{
  if (x < table.front().first) return table.front().second;

  for (std::size_t i=0; i<table.size()-1; ++i) {
    if (table[i].first < x && x < table[i+1].first) {
      auto t1 = table[i].first;
      auto y1 = table[i].second;
      auto t2 = table[i+1].first;
      auto y2 = table[i+1].second;
      return y1 + (y2-y1)/(t2-t1)*(x-t1);
    }
  }

  return table.back().second;
}
