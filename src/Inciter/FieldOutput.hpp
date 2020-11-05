// *****************************************************************************
/*!
  \file      src/Inciter/FieldOutput.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Extract field output for inciter
  \details   Extract field output for inciter.
*/
// *****************************************************************************
#ifndef FieldOutput_h
#define FieldOutput_h

#include "Types.hpp"
#include "Fields.hpp"

namespace inciter {

//! Collect node field output names based on user input
std::vector< std::string >
nodeFieldNames();

//! Collect node field output from solution based on user input
std::vector< std::vector< tk::real > >
nodeFieldOutput( const tk::Fields& Un );

} // inciter::

#endif // FieldOutput_h
