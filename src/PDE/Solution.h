// *****************************************************************************
/*!
  \file      src/PDE/Solution.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Function prototype for evaluating known solutions or ICs
  \details   This file defines a prototype used to define functions for
     evaluating known (e.g., analytical) solutions or setting initial
     conditions (ICs).
*/
// *****************************************************************************
#ifndef Solution_h
#define Solution_h

#include <vector>
#include <functional>

#include "Types.h"
#include "Keywords.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Function propotype for Problem::solution() functions
using SolutionFn = std::function<
  std::vector< real >( ncomp_t, ncomp_t, real, real, real, real ) >;

} // tk::

#endif // Solution_h
